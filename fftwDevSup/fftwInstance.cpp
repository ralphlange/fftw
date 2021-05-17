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
#include <iostream>

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
    instances.push_back(this);
    job = epicsJobCreate(workers.pool, calcJob, this);
    assert(job != nullptr);
}

void
FFTWInstance::calculate()
{
    PTimer runtime;

    for (auto conn : inputs) {
        if (conn->sigtype == FFTWConnector::InputReal) {
            fftw.set_input(conn->getNextInputValue());
        }
    }

    // number of time samples
    size_t ntime = fftw.input_sz;
    // number of frequency samples
    size_t nfreq = ntime / 2 + 1;

    assert(ntime > 0);
    assert(nfreq > 0);

    fftw.apply_window();
    runtime.maybeSnap("calculate() prepare", 5e-3);

    bool fscale_changed = fftw.replan();
    runtime.maybeSnap("calculate() replan", 0.1);

    fftw.transform();
    runtime.maybeSnap("calculate() execute", 3e-3);

    double *outr = nullptr;
    double *outi = nullptr;
    double *outm = nullptr;
    double *outp = nullptr;

    if (useReal) {
        outReal = std::unique_ptr<std::vector<double>>(new std::vector<double>(nfreq));
        outr = outReal->data();
    }
    if (useImag) {
        outImag = std::unique_ptr<std::vector<double>>(new std::vector<double>(nfreq));
        outi = outImag->data();
    }
    if (useMagn) {
        outMagn = std::unique_ptr<std::vector<double>>(new std::vector<double>(nfreq));
        outm = outMagn->data();
    }
    if (usePhas) {
        outPhas = std::unique_ptr<std::vector<double>>(new std::vector<double>(nfreq));
        outp = outPhas->data();
    }
    if (useFscale && fscale_changed) {
        outFscale = std::unique_ptr<std::vector<double>>(new std::vector<double>(fftw.fscale));
    }

    for (size_t i = 1; i < nfreq; i++) {
        fftw_complex &out = fftw.output[i];
#define creal(C) C[0]
#define cimag(C) C[1]
        if (useReal)
            outr[i] = creal(out);
        if (useImag)
            outi[i] = cimag(out);
        if (useMagn)
            outm[i] = 20. * log(sqrt(creal(out) * creal(out) + cimag(out) * cimag(out)));
        if (usePhas)
            outp[i] = atan(cimag(out) / creal(out));
#undef creal
#undef cimag
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
        default:
            break;
        }
    }
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
        std::cout << "\n  " << FFTWConnector::SignalTypeName(conn->sigtype) << ": "
                  << conn->prec->name;
    }
    for (auto &conn : outputs) {
        std::cout << "\n  " << FFTWConnector::SignalTypeName(conn->sigtype) << ": "
                  << conn->prec->name;
    }
    if (triggerSrc)
        std::cout << "\nTriggered by: " << triggerSrc->prec->name;
    else
        std::cout << "\nNo trigger set";
    std::cout << "\nInput size: " << fftw.input_sz;
    std::cout << "\nWindow type: " << FFTWCalc::WindowTypeName(fftw.wintype);
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
