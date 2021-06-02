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

#ifndef FFTWCONNECTOR_H
#define FFTWCONNECTOR_H

#include <memory>
#include <vector>

#include <dbScan.h>
#include <dbCommon.h>
#include <epicsMutex.h>
#include <epicsGuard.h>
#include <epicsTime.h>

#include "fftwCalc.h"

typedef epicsGuard<epicsMutex> Guard;
typedef epicsGuardRelease<epicsMutex> UnGuard;

class FFTWInstance;

// FFTWConnector
// - per-record configuration parameters
// - points to an FFTWInstance and a record

class FFTWConnector
{
public:
    FFTWConnector(dbCommon *prec);

    enum SignalType {
        None = 0,
        SetWindowType,
        SetSampleFreq,
        ExecutionTime,
        InputReal,
        OutputReal,
        OutputImag,
        OutputMagn,
        OutputPhas,
        OutputFscale,
        OutputWindow
    };
    enum TransformType {
        R2c_1d = 0,
    };

    static inline const char *
    SignalTypeName(const SignalType s)
    {
        switch (s) {
        case None:
            return "None";
        case SetWindowType:
            return "SetWindowType";
        case SetSampleFreq:
            return "SetSampleFreq";
        case ExecutionTime:
            return "ExecutionTime";
        case InputReal:
            return "InputReal";
        case OutputReal:
            return "OutputReal";
        case OutputImag:
            return "OutputImag";
        case OutputMagn:
            return "OutputMagn";
        case OutputPhas:
            return "OutputPhas";
        case OutputFscale:
            return "OutputFscale";
        case OutputWindow:
            return "OutputWindow";
        }
        return "<none>";
    }

    epicsMutex lock;
    FFTWInstance *inst;
    dbCommon *prec;

    SignalType sigtype;
    TransformType trftype;

    long get_ioint(int cmd, dbCommon *prec, IOSCANPVT *io);

    // Record side interface
    //     called from record processing
    //     holds record lock

    // Copy value into connector (next)
    void setNextInputValue(void *bptr, epicsUInt32 nelm);

    // Move value from connector into record
    void getNextOutputValue(void **bptr, epicsUInt32 nelm, epicsUInt32 *nord);

    // Create new value and move into record
    void createEmptyOutputValue(void **bptr, epicsUInt32 nelm);

    // Set sampling frequency
    void setSampleFreq(const double f);

    // Set window type
    void setWindowType(const FFTWCalc::WindowType t);

    // FFTW instance side interface

    // Move value from connector into instance
    std::unique_ptr<std::vector<double, FFTWAllocator<double>>> getNextInputValue();

    // Move value from instance into connector (next)
    void setNextOutputValue(std::unique_ptr<std::vector<double>> value);

    // Get the sampling frequency
    double getSampleFreq();

    // Get the window type
    FFTWCalc::WindowType getWindowType();

    // Trigger the next transform
    void trigger();

    // Set timestamp
    void setTimestamp(const epicsTimeStamp &ts);

    // Get timestamp
    epicsTimeStamp getTimestamp();

private:
    std::unique_ptr<std::vector<double>> curr_out, next_out;
    std::unique_ptr<std::vector<double, FFTWAllocator<double>>> next_inp;
    FFTWCalc::WindowType wintype;
    double fsample;
    epicsTimeStamp ts;
};

#endif // FFTWCONNECTOR_H
