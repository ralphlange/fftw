/*************************************************************************\
* Copyright (c) 2021 ITER Organization.
* This module is distributed subject to a Software License Agreement found
* in file LICENSE that is included with this distribution.
\*************************************************************************/

/*
 *  Author: Ralph Lange <ralph.lange@gmx.de>
 */

#include <iostream>
#include <string>
#include <cstring>

#include <errlog.h>
#include <iocsh.h>

#include "fftwInstance.h"

#include <epicsExport.h> // defines epicsExportSharedSymbols

namespace {

static const iocshArg fftwShowArg0 = {"instance name", iocshArgString};
static const iocshArg fftwShowArg1 = {"verbosity level [0]", iocshArgInt};

static const iocshArg *const fftwShowArg[2] = {&fftwShowArg0, &fftwShowArg1};

static const iocshFuncDef fftwShowFuncDef = {"fftwShow", 2, fftwShowArg};

static void
fftwShowCallFunc(const iocshArgBuf *args)
{
    bool ok = true;
    int verb = 0;

    if (args[0].sval == nullptr) {
        errlogPrintf("missing argument #1 (instance name)\n");
        ok = false;
    } else if (strchr(args[0].sval, ' ')) {
        errlogPrintf("invalid argument #1 (instance name) '%s'\n", args[0].sval);
        ok = false;
    }

    if (!ok)
        return;

    auto instance = FFTWInstance::find(args[0].sval);

    if (instance == nullptr) {
        errlogPrintf("'%s' no such instance\n", args[0].sval);
        ok = false;
    }

    if (args[1].ival < 0) {
        errlogPrintf("invalid argument #2 (verbosity level) '%d'\n", args[1].ival);
    } else {
        verb = args[1].ival;
    }

    if (ok)
        instance->show(verb);
}

static void
fftwIocshRegister()
{
    iocshRegister(&fftwShowFuncDef, fftwShowCallFunc);
}

extern "C" {
epicsExportRegistrar(fftwIocshRegister);
}

} // namespace
