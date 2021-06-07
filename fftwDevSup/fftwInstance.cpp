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
    , lasttime(0.0)
    , valid(false)
    , triggerSrc(nullptr)
    , useReal(false)
    , useImag(false)
    , useMagn(false)
    , usePhas(false)
    , useFscale(false)
    , useWindow(false)
    , sizeReal(0)
    , sizeImag(0)
    , sizeMagn(0)
    , sizePhas(0)
    , sizeFscale(0)
    , sizeWindow(0)
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

    epicsTimeStamp ts = triggerSrc->getTimestamp();

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

    if (useReal || useImag) {
        outReal = std::shared_ptr<std::vector<double>>(new std::vector<double>(fftw.nfreq));
        outReal->reserve(sizeReal);
        double *outr = outReal->data();
        outImag = std::shared_ptr<std::vector<double>>(new std::vector<double>(fftw.nfreq));
        outImag->reserve(sizeImag);
        double *outi = outImag->data();

        for (size_t i = 0; i < fftw.nfreq; i++) {
            fftw_complex &out = fftw.output[i];
            outr[i] = creal(out);
            outi[i] = cimag(out);
        }
    }

    if (useMagn || usePhas) {
        outMagn = std::shared_ptr<std::vector<double>>(new std::vector<double>(fftw.nfreq));
        outMagn->reserve(sizeMagn);
        double *outm = outMagn->data();
        outPhas = std::shared_ptr<std::vector<double>>(new std::vector<double>(fftw.nfreq));
        outPhas->reserve(sizePhas);
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
        outWindow = std::shared_ptr<std::vector<double>>(new std::vector<double>(fftw.ntime));
        outWindow->reserve(sizeWindow);
        double *outw = outWindow->data();
        double *getw = fftw.window.data();
        for (size_t i = 0; i < fftw.ntime; i++)
            outw[i] = getw[i];
    }

    if (useFscale && fscale_changed) {
        outFscale = std::shared_ptr<std::vector<double>>(new std::vector<double>(fftw.nfreq));
        outFscale->reserve(sizeFscale);
        double *outf = outFscale->data();
        double *getf = fftw.fscale.data();
        for (size_t i = 0; i < fftw.nfreq; i++)
            outf[i] = getf[i];
    }

    runtime.maybeSnap("calculate() post-proc", 1e-3);

    lasttime = calctime.snap();
    for (auto conn : outputs) {
        conn->setTimestamp(ts);
        switch (conn->sigtype) {
        case FFTWConnector::OutputImag:
            conn->setNextOutputValue(outImag);
            break;
        case FFTWConnector::OutputReal:
            conn->setNextOutputValue(outReal);
            break;
        case FFTWConnector::OutputMagn:
            conn->setNextOutputValue(outMagn);
            break;
        case FFTWConnector::OutputPhas:
            conn->setNextOutputValue(outPhas);
            break;
        case FFTWConnector::OutputFscale:
            if (fscale_changed)
                conn->setNextOutputValue(outFscale);
            break;
        case FFTWConnector::OutputWindow:
            if (window_changed)
                conn->setNextOutputValue(outWindow);
            break;
        case FFTWConnector::ExecutionTime:
            conn->setRuntime(lasttime);
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
    calctime.start();
    epicsJobQueue(job);
}

void
FFTWInstance::show(const unsigned int verbosity) const
{
    std::cout << "Instance " << name << " using job " << job << " of pool " << workers.pool
              << "\nConnected records:";
    for (auto &conn : inputs)
        conn->show(verbosity, 2);
    for (auto &conn : outputs)
        conn->show(verbosity, 2);
    if (verbosity > 1) {
        std::cout << "\nRequired output vector sizes:\n ";
        if (useReal)
            std::cout << " Real:" << sizeReal;
        if (useImag)
            std::cout << " Imag:" << sizeImag;
        if (useMagn)
            std::cout << " Magn:" << sizeMagn;
        if (usePhas)
            std::cout << " Phas:" << sizePhas;
        if (useFscale)
            std::cout << " Fscale:" << sizeFscale;
        if (useWindow)
            std::cout << " Window:" << sizeWindow;
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

// Careful: not thread safe (ok during record initialization)
void FFTWInstance::setRequiredOutputSize(const FFTWConnector::SignalType type, const epicsUInt32 size)
{
    switch (type) {
    case FFTWConnector::OutputReal:
        if (size > sizeReal)
            sizeReal = size;
        break;
    case FFTWConnector::OutputImag:
        if (size > sizeImag)
            sizeImag = size;
        break;
    case FFTWConnector::OutputMagn:
        if (size > sizeMagn)
            sizeMagn = size;
        break;
    case FFTWConnector::OutputPhas:
        if (size > sizePhas)
            sizePhas = size;
        break;
    case FFTWConnector::OutputFscale:
        if (size > sizeFscale)
            sizeFscale = size;
        break;
    case FFTWConnector::OutputWindow:
        if (size > sizeWindow)
            sizeWindow = size;
        break;
    default:
        break;
    }
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
