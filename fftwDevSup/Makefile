TOP=..

include $(TOP)/configure/CONFIG
#----------------------------------------
#  ADD MACRO DEFINITIONS AFTER THIS LINE

USR_CPPFLAGS += -DUSE_TYPED_RSET
USR_CPPFLAGS += -I$(FFTW3)

#=============================

EXPANDVARS += EPICS_FFTW_MAJOR_VERSION
EXPANDVARS += EPICS_FFTW_MINOR_VERSION
EXPANDVARS += EPICS_FFTW_MAINTENANCE_VERSION
EXPANDVARS += EPICS_FFTW_DEVELOPMENT_FLAG

EXPANDFLAGS += $(foreach var,$(EXPANDVARS),-D$(var)="$(strip $($(var)))")

LIBRARY_IOC += fftwSup

DBD += fftw.dbd

INC += fftwVersionNum.h
INC += fftwVersion.h

fftwSup_SRCS += fftwSupport.cpp
fftwSup_SRCS += fftwConnector.cpp
fftwSup_SRCS += fftwInstance.cpp
fftwSup_SRCS += fftwCalc.cpp
fftwSup_SRCS += iocshIntegration.cpp

fftwSup_LIBS += $(EPICS_BASE_IOC_LIBS)
fftwSup_SYS_LIBS_Linux += fftw3

fftwSup_LIBS_WIN32 += libfftw3-3
libfftw3-3_DIR = $(FFTW3)

SHRLIB_VERSION ?= $(EPICS_FFTW_MAJOR_VERSION).$(EPICS_FFTW_MINOR_VERSION)

#===========================

include $(TOP)/configure/RULES
#----------------------------------------
#  ADD RULES AFTER THIS LINE

# Can't use EXPAND as generated headers must appear
# in O.Common, but EXPAND emits rules for O.$(T_A)
../O.Common/fftwVersionNum.h: ../fftwVersionNum.h@
	$(EXPAND_TOOL) $(EXPANDFLAGS) $($@_EXPANDFLAGS) $< $@
