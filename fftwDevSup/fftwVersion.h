/*************************************************************************\
* Copyright (c) 2021 ITER Organization.
* This module is distributed subject to a Software License Agreement found
* in file LICENSE that is included with this distribution.
\*************************************************************************/

/*
 *  Author: Ralph Lange <ralph.lange@gmx.de>
 *
 *  based on code by Michael Davidsaver <mdavidsaver@ospreydcs.com>
 */

#ifndef FFTWVERSION_H
#define FFTWVERSION_H

#include <fftw3.h>
#include <epicsVersion.h>

#ifndef VERSION_INT
#  define VERSION_INT(V,R,M,P) ( ((V)<<24) | ((R)<<16) | ((M)<<8) | (P))
#endif

/* include generated headers with:
 *   EPICS_FFTW_MAJOR_VERSION
 *   EPICS_FFTW_MINOR_VERSION
 *   EPICS_FFTW_MAINTENANCE_VERSION
 *   EPICS_FFTW_DEVELOPMENT_FLAG
 */
#include "fftwVersionNum.h"

#define FFTW_VERSION_INT VERSION_INT(EPICS_FFTW_MAJOR_VERSION, EPICS_FFTW_MINOR_VERSION, EPICS_FFTW_MAINTENANCE_VERSION, 0)

#define EPICS_FFTW_LIBRARY_VERSION fftw_version

#endif // FFTWVERSION_H
