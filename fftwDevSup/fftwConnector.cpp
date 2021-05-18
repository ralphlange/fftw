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

#include <dbScan.h>
#include <dbCommon.h>

#include "fftwCalc.h"
#include "fftwInstance.h"
#include "fftwConnector.h"

FFTWConnector::FFTWConnector(dbCommon *prec)
    : inst(nullptr)
    , prec(prec)
    , sigtype(None)
{}

long
FFTWConnector::get_ioint(int cmd, dbCommon *prec, IOSCANPVT *io)
{
    switch (sigtype) {
    case OutputReal:
    case OutputImag:
    case OutputMagn:
    case OutputPhas:
        *io = inst->valueScan;
        return 0;
    case OutputFscale:
        Guard(inst->lock);
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
FFTWConnector::setNextInputValue(void *bptr, epicsUInt32 nelm)
{
    std::unique_ptr<std::vector<double, FFTWAllocator<double>>> vec(new std::vector<double, FFTWAllocator<double>>(nelm));

    const double *src = static_cast<double *>(bptr);
    vec->insert(vec->end(), src, src + nelm);
    next_inp = std::move(vec);
}

std::unique_ptr<std::vector<double, FFTWAllocator<double>>>
FFTWConnector::getNextInputValue()
{
    return std::move(next_inp);
}

void
FFTWConnector::setNextOutputValue(std::unique_ptr<std::vector<double>> value)
{
    next_out = std::move(value);
}

void
FFTWConnector::getNextOutputValue(void **bptr, epicsUInt32 nelm, epicsUInt32 *nord)
{
    if (next_out) {
        epicsUInt32 N = std::min<epicsUInt32>(nelm, (epicsUInt32) next_out->size());
        if (N < nelm) {
            next_out->resize(nelm);
            *nord = N;
        } else {
            *nord = nelm;
        }
        *bptr = next_out.get();
        curr_out = std::move(next_out);
    }
}

void
FFTWConnector::createEmptyOutputValue(void **bptr, epicsUInt32 nelm)
{
    std::unique_ptr<std::vector<double>> vec(new std::vector<double>(nelm));
    curr_out = std::move(vec);
    *bptr = curr_out.get();
}

void
FFTWConnector::trigger()
{
    inst->trigger();
}
