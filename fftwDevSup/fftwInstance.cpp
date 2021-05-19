/*************************************************************************\
* Copyright (c) 2015 Brookhaven Science Assoc. as operator of
*     Brookhaven National Laboratory.
* Copyright (c) 2021 ITER Organization.
* This module is distributed subject to a Software License Agreement found
* in file LICENSE that is included with this distribution.
\*************************************************************************/

/*
 *  Author: Ralph Lange <ralph.lange@gmx.de>
 *
 *  based on pscdrv/sigApp by Michael Davidsaver <mdavidsaver@ospreydcs.com>
 */

#include <cmath>
#include <cstring>
#include <iostream>
#include <iomanip>

#include <dbScan.h>
#include <epicsThread.h>

#include "fftwConnector.h"
#include "fftwInstance.h"

// Windows implementation of clock_gettime
// see: https://stackoverflow.com/questions/5404277/porting-clock-gettime-to-windows
#ifdef _WIN32
#include <Windows.h>
#include <minwinbase.h>
//struct timespec { long tv_sec; long tv_nsec; };    //header part
int clock_gettime(int, struct timespec *spec)      //C-file part
{  __int64 wintime; GetSystemTimeAsFileTime((FILETIME*)&wintime);
    wintime      -=116444736000000000i64;  //1jan1601 to 1jan1970
    spec->tv_sec  =wintime / 10000000i64;           //seconds
    spec->tv_nsec =wintime % 10000000i64 *100;      //nano-seconds
    return 0;
}
#endif

std::vector<FFTWInstance *> FFTWInstance::instances;
FFTWThreadPool FFTWInstance::workers;

FFTWThreadPool::FFTWThreadPool()
{
    epicsThreadPoolConfigDefaults(&poolConfig);
    pool = epicsThreadPoolCreate(&poolConfig);
    assert(pool != nullptr);
}

FFTWThreadPool::~FFTWThreadPool()
{
    epicsThreadPoolDestroy(pool);
}

FFTWInstance::FFTWInstance(const std::string &name)
    : name(name)
    , valid(false)
    , triggerSrc(nullptr)
{
    scanIoInit(&valueScan);
    scanIoInit(&scaleScan);
    scanIoInit(&windowScan);
    instances.push_back(this);
    job = epicsJobCreate(workers.pool, calcJob, this);
    assert(job != nullptr);
}

void
FFTWInstance::calculate()
{
    PTimer runtime;

    for (auto conn : inputs) {
        switch (conn->sigtype) {
        case FFTWConnector::InputReal:
            fftw.set_input_real(conn->getNextInputValue());
            break;
        case FFTWConnector::SetSampleFreq:
            fftw.set_fsamp(conn->getSampleFreq());
            break;
        case FFTWConnector::SetWindowType:
            fftw.set_wtype(conn->getWindowType());
            break;
        default:
            break;
        }
    }

    bool window_changed = fftw.apply_window();
    runtime.maybeSnap("calculate() prepare", 5e-3);

    bool fscale_changed = fftw.replan();
    runtime.maybeSnap("calculate() replan", 0.1);

    fftw.transform();
    runtime.maybeSnap("calculate() execute", 3e-3);

    valid = true;
    if (fftw.output.size() == 0 || fftw.window.size() == 0 || fftw.fscale.size() == 0)
        valid = false;

    // Trying to do some optimization while letting the compiler still do vectorization
    // - always do the simple transactions
    // - do a second loop with the complex transactions if required

#define creal(C) C[0]
#define cimag(C) C[1]

    if (!outReal)
        outReal = std::unique_ptr<std::vector<double>>(new std::vector<double>(fftw.nfreq));
    double *outr = outReal->data();
    if (!outImag)
        outImag = std::unique_ptr<std::vector<double>>(new std::vector<double>(fftw.nfreq));
    double *outi = outImag->data();

    for (size_t i = 0; i < fftw.nfreq; i++) {
        fftw_complex &out = fftw.output[i];
        outr[i] = creal(out);
        outi[i] = cimag(out);
    }

    if (useMagn || usePhas) {
        outMagn = std::unique_ptr<std::vector<double>>(new std::vector<double>(fftw.nfreq));
        double *outm = outMagn->data();
        outPhas = std::unique_ptr<std::vector<double>>(new std::vector<double>(fftw.nfreq));
        double *outp = outPhas->data();

        for (size_t i = 0; i < fftw.nfreq; i++) {
            fftw_complex &out = fftw.output[i];
            outm[i] = 20. * log(sqrt(creal(out) * creal(out) + cimag(out) * cimag(out)));
            outp[i] = atan(cimag(out) / creal(out));
        }
    }

#undef creal
#undef cimag

    if (useWindow && window_changed) {
        outWindow = std::unique_ptr<std::vector<double>>(new std::vector<double>(fftw.ntime));
        double *outw = outWindow->data();
        double *getw = fftw.window.data();
        for (size_t i = 0; i < fftw.ntime; i++)
            outw[i] = getw[i];
    }

    if (useFscale && fscale_changed) {
        outFscale = std::unique_ptr<std::vector<double>>(new std::vector<double>(fftw.nfreq));
        double *outf = outFscale->data();
        double *getf = fftw.fscale.data();
        for (size_t i = 0; i < fftw.nfreq; i++)
            outf[i] = getf[i];
    }

    runtime.maybeSnap("calculate() post-proc", 1e-3);

    for (auto conn : outputs) {
        switch (conn->sigtype) {
        case FFTWConnector::OutputImag:
            conn->setNextOutputValue(std::move(outImag));
            break;
        case FFTWConnector::OutputReal:
            conn->setNextOutputValue(std::move(outReal));
            break;
        case FFTWConnector::OutputMagn:
            conn->setNextOutputValue(std::move(outMagn));
            break;
        case FFTWConnector::OutputPhas:
            conn->setNextOutputValue(std::move(outPhas));
            break;
        case FFTWConnector::OutputFscale:
            if (fscale_changed)
                conn->setNextOutputValue(std::move(outFscale));
            break;
        case FFTWConnector::OutputWindow:
            if (window_changed)
                conn->setNextOutputValue(std::move(outWindow));
            break;
        default:
            break;
        }
    }

    scanIoRequest(valueScan);
    if (fscale_changed)
        scanIoRequest(scaleScan);
    if (window_changed)
        scanIoRequest(windowScan);
}

void
FFTWInstance::trigger()
{
    if (triggerSrc && triggerSrc->prec->tpro > 5)
        std::cerr << "Queueing calculation job for " << name << std::endl;
    epicsJobQueue(job);
}

void
FFTWInstance::show(const unsigned int verbosity) const
{
    std::cout << "Instance " << name << " using job " << job << " of pool " << workers.pool
              << "\nConnected records:";
    for (auto &conn : inputs) {
        const char *type = FFTWConnector::SignalTypeName(conn->sigtype);
        std::cout << "\n  " << type << std::setw(15 - strlen(type)) << std::left << ": "
                  << std::right << conn->prec->name;
    }
    for (auto &conn : outputs) {
        const char *type = FFTWConnector::SignalTypeName(conn->sigtype);
        std::cout << "\n  " << type << std::setw(15 - strlen(type)) << std::left << ": "
                  << std::right << conn->prec->name;
    }
    if (triggerSrc)
        std::cout << "\nTriggered by: " << triggerSrc->prec->name;
    else
        std::cout << "\nNo trigger set";
    std::cout << "\nInput size: " << fftw.input_sz
              << "\nWindow type: " << FFTWCalc::WindowTypeName(fftw.wintype)
              << "\nSample freq: " << fftw.fsamp
              << "\nExec time: " << lasttime;
    std::cout << std::endl;
}

FFTWInstance *
FFTWInstance::findOrCreate(const std::string &name)
{
    if (FFTWInstance *inst = find(name))
        return inst;
    else
        return new FFTWInstance(name);
}

FFTWInstance *
FFTWInstance::find(const std::string &name)
{
    for (auto it : instances)
        if (it->name == name)
            return it;
    return nullptr;
}

void
FFTWInstance::calcJob(void *arg, epicsJobMode mode)
{
    auto instance = reinterpret_cast<FFTWInstance *>(arg);
    if (mode == epicsJobModeCleanup) {
        epicsJobDestroy(instance->job);
        return;
    }
    if (FFTWDebug)
        std::cerr << "Running calculation for instance " << instance->name << std::endl;
    instance->calculate();
}
