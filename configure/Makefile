TOP=..

include $(TOP)/configure/CONFIG

CFG += CONFIG_FFTW_VERSION

TARGETS = $(CONFIG_TARGETS)
CONFIGS += $(subst ../,,$(wildcard $(CONFIG_INSTALLS)))

include $(TOP)/configure/RULES
