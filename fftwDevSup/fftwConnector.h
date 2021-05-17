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

#include "fftwCalc.h"

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

    FFTWInstance *inst;
    dbCommon *prec;

    SignalType sigtype;
    TransformType trftype;

    long get_ioint(int cmd, dbCommon *prec, IOSCANPVT *io);

    // Copy value into connector (next)
    //     must be called under the record lock
    //     called from record processing
    void setNextInputValue(void *bptr, epicsUInt32 nelm);
    // Move value from connector into instance
    //     called from the FFTW Instance
    std::unique_ptr<std::vector<double, FFTWAllocator<double>>> getNextInputValue();

    // Move value into connector (next)
    //     called from the FFTW Instance
    void setNextOutputValue(std::unique_ptr<std::vector<double>> value);
    // Move value from connector into record
    //     must be called under the record lock
    //     called from record processing
    void getNextOutputValue(void **bptr, epicsUInt32 nelm, epicsUInt32 *nord);
    // Create new value and move into record
    //     must be called under the record lock
    //     called from record initialization
    void createEmptyOutputValue(void **bptr, epicsUInt32 nelm);

    // Trigger the next transform
    void trigger();

private:
    std::unique_ptr<std::vector<double>> curr_out, next_out;
    std::unique_ptr<std::vector<double, FFTWAllocator<double>>> next_inp;
};

#endif // FFTWCONNECTOR_H
