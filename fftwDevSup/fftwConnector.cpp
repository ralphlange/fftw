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

#include <memory>
#include <utility>
#include <algorithm>
#include <cstring>
#include <iostream>
#include <iomanip>

#include <dbScan.h>
#include <dbCommon.h>

#include "fftwCalc.h"
#include "fftwInstance.h"
#include "fftwConnector.h"

FFTWConnector::FFTWConnector(dbCommon *prec)
    : inst(nullptr)
    , prec(prec)
    , sigtype(None)
    , offset(0)
{}

long
FFTWConnector::get_ioint(int cmd, dbCommon *prec, IOSCANPVT *io)
{
    switch (sigtype) {
    case OutputReal:
    case OutputImag:
    case OutputMagn:
    case OutputPhas:
    case ExecutionTime:
        *io = inst->valueScan;
        return 0;
    case OutputFscale:
        *io = inst->scaleScan;
        return 0;
    case OutputWindow:
        *io = inst->windowScan;
        return 0;
    default:
        return 1;
    }
}

void
FFTWConnector::show(const unsigned int verbosity, const unsigned char indent) const
{
    const char *type = FFTWConnector::SignalTypeName(sigtype);
    size_t tlen = strlen(type);
    std::cout << "\n"
              << std::setw(indent + tlen) << type << std::setw(15 - tlen) << std::left << ": "
              << std::right << prec->name;
    if (offset)
        std::cout << " offset=" << offset;
}

void
FFTWConnector::setRequiredOutputSize(const epicsUInt32 nelm)
{
    inst->setRequiredOutputSize(sigtype, nelm + offset);
}

void
FFTWConnector::setNextInputValue(void *bptr, epicsUInt32 elements)
{
    std::unique_ptr<std::vector<double, FFTWAllocator<double>>> vec(new std::vector<double, FFTWAllocator<double>>());

    const double *src = static_cast<double *>(bptr);
    vec->insert(vec->end(), src, src + elements);
    next_inp = std::move(vec);
}

std::unique_ptr<std::vector<double, FFTWAllocator<double>>>
FFTWConnector::getNextInputValue()
{
    return std::move(next_inp);
}

void
FFTWConnector::setNextOutputValue(std::shared_ptr<std::vector<double>> value)
{
    next_out = std::move(value);
}

void
FFTWConnector::getNextOutputValue(void **bptr, epicsUInt32 nelm, epicsUInt32 *nord)
{
    std::shared_ptr<std::vector<double>> vec = std::move(next_out);
    if (vec) {
        assert(vec->capacity() >= nelm + offset);
        *nord = static_cast<epicsUInt32>(std::min<size_t>(nelm, vec->size() - offset));
        double *data = vec.get()->data();
        if (offset && vec->size() > offset)
            data += offset;
        *bptr = data;
        curr_out = std::move(vec);
    }
}

void
FFTWConnector::createEmptyOutputValue(void **bptr, epicsUInt32 nelm)
{
    std::shared_ptr<std::vector<double>> vec(new std::vector<double>(nelm));
    *bptr = vec.get()->data();
    curr_out = std::move(vec);
}

void
FFTWConnector::setSampleFreq(const double f)
{
    Guard G(lock);
    fsample = f;
}

double FFTWConnector::getSampleFreq()
{
    Guard G(lock);
    return fsample;
}

void
FFTWConnector::setWindowType(const FFTWCalc::WindowType t)
{
    Guard G(lock);
    wintype = t;
}

double
FFTWConnector::getRuntime()
{
    return runtime;
}

void
FFTWConnector::setOffset(const size_t offset)
{
    this->offset = offset;
}

FFTWCalc::WindowType FFTWConnector::getWindowType()
{
    Guard G(lock);
    return wintype;
}

void
FFTWConnector::trigger()
{
    inst->trigger();
}

void
FFTWConnector::setRuntime(const double time)
{
    Guard G(lock);
    runtime = time;
}

void
FFTWConnector::setTimestamp(const epicsTimeStamp &ts)
{
    Guard G(lock);
    this->ts = ts;
}

epicsTimeStamp FFTWConnector::getTimestamp()
{
    return ts;
}
