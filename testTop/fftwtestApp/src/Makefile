TOP=../..

include $(TOP)/configure/CONFIG
#----------------------------------------
#  ADD MACRO DEFINITIONS AFTER THIS LINE
#=============================

#=============================
# Build the IOC application

PROD_IOC = fftwtest

DBD += fftwtest.dbd
fftwtest_DBD += base.dbd
fftwtest_DBD += fftw.dbd

fftwtest_LIBS += fftwSup

# fftwtest_registerRecordDeviceDriver.cpp derives from fftwtest.dbd
fftwtest_SRCS += fftwtest_registerRecordDeviceDriver.cpp

# Build the main IOC entry point on workstation OSs.
fftwtest_SRCS_DEFAULT += fftwtestMain.cpp
fftwtest_SRCS_vxWorks += -nil-

# Finally link to the EPICS Base libraries
fftwtest_LIBS += $(EPICS_BASE_IOC_LIBS)

#===========================

include $(TOP)/configure/RULES
#----------------------------------------
#  ADD RULES AFTER THIS LINE
