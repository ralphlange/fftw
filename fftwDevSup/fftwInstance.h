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

#ifndef FFTWINSTANCE_H
#define FFTWINSTANCE_H

#include <string>
#include <vector>

#include <fftw3.h>

#include <dbScan.h>
#include <epicsGuard.h>
#include <epicsMutex.h>
#include <epicsThreadPool.h>
#include <epicsTime.h>

#include "fftwConnector.h"
#include "fftwCalc.h"

// Windows implementation of clock_gettime
#ifdef _WIN32
#define CLOCK_PROCESS_CPUTIME_ID 0
#    include <Windows.h>
#    include <minwinbase.h>
int clock_gettime(int, struct timespec *spec);
#endif

// Performance timer
struct PTimer {
    timespec tstart;
    PTimer() {start();}
    void start()
    {
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &tstart);
    }
    double snap()
    {
        timespec now;
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &now);
        double ret = now.tv_sec-tstart.tv_sec + 1e-9*(now.tv_nsec-tstart.tv_nsec);
        tstart = now;
        return ret;
    }
    // Print msg if elapsed time is greater than the given threshold
    void maybeSnap(const char *msg, double threshold=0.0)
    {
        double interval = snap();
        if(FFTWDebug && interval>threshold) {
            errlogPrintf("%s over threshold %f > %f\n", msg, interval, threshold);
        }
    }
};

typedef epicsGuard<epicsMutex> Guard;
typedef epicsGuardRelease<epicsMutex> UnGuard;

class FFTWConnector;

//typedef std::vector<double, FFTWAllocator<double>> FFTWvector_d;
//typedef std::vector<fftw_complex, FFTWAllocator<fftw_complex>> FFTWvector_c;

struct FFTWThreadPool
{
    FFTWThreadPool();
    ~FFTWThreadPool();
    epicsThreadPool *pool;
    epicsThreadPoolConfig poolConfig;
};

class FFTWInstance
{
public:
    std::string name;
    epicsMutex lock;

    epicsTimeStamp timein, timeout;
    double lasttime;
    bool valid;

    // Pointers to the trigger, the list of inputs and outputs
    FFTWConnector *triggerSrc;
    std::vector<FFTWConnector *> inputs;
    std::vector<FFTWConnector *> outputs;

    std::unique_ptr<std::vector<double>> outReal, outImag, outMagn, outPhas, outFscale, outWindow;
    bool useReal, useImag, useMagn, usePhas, useFscale, useWindow;

    FFTWCalc fftw;

    IOSCANPVT valueScan, scaleScan, windowScan;

    void trigger();

    // Show method to print the setup
    void show(const unsigned int verbosity) const;

    // Find an instance
    static FFTWInstance *find(const std::string &name);

    // Factory method to create an instance
    static FFTWInstance *findOrCreate(const std::string &name);

    // epicsThreadPool interface
    static void calcJob(void *arg, epicsJobMode mode);

private:
    FFTWInstance(const std::string &name);

    // Transformation routine called from the job
    void calculate();

    static std::vector<FFTWInstance *> instances;
    static FFTWThreadPool workers;
    epicsJob *job;
};

#endif // FFTWINSTANCE_H
