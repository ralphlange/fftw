#!../../bin/linux-x86_64/fftwtest

## You may have to change fftwtest to something else
## everywhere it appears in this file

#< envPaths

## Register all support components
dbLoadDatabase("../../dbd/fftwtest.dbd",0,0)
fftwtest_registerRecordDeviceDriver(pdbbase) 

## Load record instances
dbLoadRecords("../../db/single.db","P=A1,R=:")

iocInit()

## Start any sequence programs
#seq sncfftwtest,"user=ralph"
