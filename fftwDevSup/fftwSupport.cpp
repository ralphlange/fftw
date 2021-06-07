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

#include <string>
#include <cstring>
#include <set>
#include <vector>
#include <iostream>
#include <sstream>

#include <epicsEvent.h>
#include <epicsGuard.h>
#include <epicsThread.h>
#include <errlog.h>
#include <alarm.h>
#include <dbCommon.h>
#include <dbStaticLib.h>
#include <dbAccess.h>
#include <dbScan.h>
#include <devLib.h> // for S_dev_* codes
#include <devSup.h>
#include <recGbl.h>
#include <menuConvert.h>
#include <menuFtype.h>

#include <mbboRecord.h>
#include <aoRecord.h>
#include <aiRecord.h>
#include <aaoRecord.h>
#include <aaiRecord.h>

#include "fftwVersion.h"
#include "fftwConnector.h"
#include "fftwInstance.h"

namespace {

template<typename R>
struct dset6
{
    long N;
    long (*report)(R *);
    long (*init)(int);
    long (*init_record)(R *);
    long (*get_io_intr_info)(int, dbCommon *prec, IOSCANPVT *);
    long (*readwrite)(R *);
    long (*linconv)(R *);
};

struct SB
{
    std::ostringstream strm;
    operator std::string() { return strm.str(); }
    template<typename T>
    SB &
    operator<<(const T &v)
    {
        strm << v;
        return *this;
    }
};

class DBEntry {
    DBENTRY entry;
public:
    DBENTRY *pentry() const { return const_cast<DBENTRY*>(&entry); }
    explicit DBEntry(dbCommon *prec) {
        dbInitEntry(pdbbase, &entry);
        if (dbFindRecord(&entry, prec->name))
            throw std::logic_error(SB() << "can't find record " << prec->name);
    }
    DBEntry(const DBEntry& ent) {
        dbCopyEntryContents(const_cast<DBENTRY*>(&ent.entry), &entry);
    }
    DBEntry& operator=(const DBEntry& ent) {
        dbFinishEntry(&entry);
        dbCopyEntryContents(const_cast<DBENTRY*>(&ent.entry), &entry);
        return *this;
    }
    ~DBEntry() {
        dbFinishEntry(&entry);
    }
    DBLINK *getDevLink() const {
        if (dbFindField(pentry(), "INP") && dbFindField(pentry(), "OUT"))
            throw std::logic_error(SB() << entry.precnode->recordname << " has no INP/OUT?!");
        if (entry.pflddes->field_type != DBF_INLINK &&
            entry.pflddes->field_type != DBF_OUTLINK)
            throw std::logic_error(SB() << entry.precnode->recordname << " INP/OUT is not a link?!");
        return static_cast<DBLINK*>(entry.pfield);
    }
    bool isOutput() const {
        return !dbFindField(pentry(), "OUT");
    }
    const char *info(const char *name, const char *deflt) const
    {
        if (dbFindInfo(pentry(), name))
            return deflt;
        else
            return entry.pinfonode->string;
    }
};

template<typename F, typename I, typename R>
F
analogRaw2EGU(R *prec, I rval)
{
    F v = static_cast<F>(rval) + static_cast<F>(prec->roff);

    if (prec->aslo != 0.0)
        v *= prec->aslo;
    v += prec->aoff;

    switch (prec->linr) {
    case menuConvertLINEAR:
    case menuConvertSLOPE:
        v += prec->eoff;
        v *= prec->eslo;
    }

    return v;
}

template<typename I, typename F, typename R>
I
analogEGU2Raw(R *prec, F v)
{
    switch (prec->linr) {
    case menuConvertLINEAR:
    case menuConvertSLOPE:
        v -= prec->eoff;
        v /= prec->eslo;
    }

    v -= prec->aoff;
    if (prec->aslo != 0.0)
        v /= prec->aslo;

    v -= prec->roff;

    return static_cast<I>(v);
}

inline FFTWConnector::SignalType
signalTypeIndex(const std::string &name)
{
    if (name == "input-real")
        return FFTWConnector::InputReal;
    else if (name == "input-imag")
        return FFTWConnector::InputReal;
    else if (name == "windowtype")
        return FFTWConnector::SetWindowType;
    else if (name == "sample-freq")
        return FFTWConnector::SetSampleFreq;
    else if (name == "exectime")
        return FFTWConnector::ExecutionTime;
    else if (name == "output-real")
        return FFTWConnector::OutputReal;
    else if (name == "output-imag")
        return FFTWConnector::OutputImag;
    else if (name == "output-magn")
        return FFTWConnector::OutputMagn;
    else if (name == "output-phas")
        return FFTWConnector::OutputPhas;
    else if (name == "output-fscale")
        return FFTWConnector::OutputFscale;
    else if (name == "output-window")
        return FFTWConnector::OutputWindow;
    else
        return FFTWConnector::None;
}

// Helper: checks if c is one of the characters meaning "yes" or "no"
bool
isYes (const char c)
{
    if (strchr("YyTt1", c))
        return true;
    else if (strchr("NnFf0", c))
        return false;
    else
        throw std::runtime_error(SB() << "illegal value '" << c << "'");
}

// Helper: split string at delimiter, return list of tokens
std::vector<std::string>
splitString(const std::string &str, const char delim)
{
    std::vector<std::string> tokens;
    size_t prev = 0, sep = 0;
    do {
        sep = str.find_first_of(delim, prev);
        if (sep == std::string::npos)
            sep = str.length();
        std::string token = str.substr(prev, sep - prev);
        // allow escaping delimiters
        while (sep < str.length() && sep > 0 && str[sep - 1] == '\\') {
            prev = sep + 1;
            sep = str.find_first_of(delim, prev);
            if (sep == std::string::npos)
                sep = str.length();
            token.pop_back();
            token.append(str.substr(prev - 1, sep - prev + 1));
        }
        tokens.push_back(token);
        prev = sep + 1;
    } while (sep < str.length() && prev <= str.length());
    return tokens;
}

// Link parameter parser
FFTWConnector *
parseLink(dbCommon *prec, const DBEntry &ent)
{
    DBLINK *link = ent.getDevLink();
    std::unique_ptr<FFTWConnector> conn (new FFTWConnector(prec));

    if (link->type != INST_IO)
        throw std::logic_error("link is not INST_IO");

    if (prec->tpro > 10)
        std::cerr << prec->name << ": parsing link '" << link->value.instio.string << "'"
                  << std::endl;

    std::vector<std::string> tokens = splitString(link->value.instio.string, ' ');

    for (const auto &token : tokens) {
        // First token: instance name
        if (!conn->inst) {
            conn->inst = FFTWInstance::findOrCreate(token);
            continue;
        }
        // Second token: signal type
        if (conn->sigtype == FFTWConnector::None) {
            FFTWConnector::SignalType sig = signalTypeIndex(token);
            if (sig != FFTWConnector::None) {
                conn->sigtype = sig;
                switch (sig) {
                case FFTWConnector::InputReal:
                case FFTWConnector::SetWindowType:
                case FFTWConnector::SetSampleFreq:
                    conn->inst->inputs.push_back(conn.get());
                    break;
                case FFTWConnector::OutputReal:
                    conn->inst->outputs.push_back(conn.get());
                    conn->inst->useReal = true;
                    break;
                case FFTWConnector::OutputImag:
                    conn->inst->outputs.push_back(conn.get());
                    conn->inst->useImag = true;
                    break;
                case FFTWConnector::OutputMagn:
                    conn->inst->outputs.push_back(conn.get());
                    conn->inst->useMagn = true;
                    break;
                case FFTWConnector::OutputPhas:
                    conn->inst->outputs.push_back(conn.get());
                    conn->inst->usePhas = true;
                    break;
                case FFTWConnector::OutputFscale:
                    conn->inst->outputs.push_back(conn.get());
                    conn->inst->useFscale = true;
                    break;
                case FFTWConnector::OutputWindow:
                    conn->inst->outputs.push_back(conn.get());
                    conn->inst->useWindow = true;
                    break;
                case FFTWConnector::ExecutionTime:
                    conn->inst->outputs.push_back(conn.get());
                    break;
                case FFTWConnector::None:
                    break;
                }
                if (prec->tpro > 1)
                    std::cerr << prec->name << ": connected as "
                              << FFTWConnector::SignalTypeName(conn->sigtype)
                              << " to FFTW instance '" << conn->inst->name << "'" << std::endl;
                continue;
            }
        }
        // All others: option=value
        std::vector<std::string> option = splitString(token, '=');
        if (option[0] == "trigger") {
            if (isYes(option[1][0]) && !conn->inst->triggerSrc) {
                conn->inst->triggerSrc = conn.get();
                if (prec->tpro > 1)
                    std::cerr << prec->name << ": this record will trigger FFTW instance '"
                              << conn->inst->name << "'" << std::endl;
            }
            continue;
        }
    }
    return conn.release();
}

// Device Support routines

#define TRY \
    if (!prec->dpvt) \
        return 0; \
    FFTWConnector *conn = static_cast<FFTWConnector *>(prec->dpvt); \
    try
#define CATCH(NAME) \
    catch (std::exception & e) \
    { \
        std::cerr << prec->name << " (" #NAME ") Error : " << e.what() << std::endl; \
        (void) recGblSetSevr(prec, COMM_ALARM, INVALID_ALARM); \
        return S_dev_badRequest; \
    }
#define CHKVALID \
    if (!conn->inst->valid) { \
        (void) recGblSetSevr(prec, READ_ALARM, INVALID_ALARM); \
        return 0; \
    }

long
global_init(int pass)
{
    if (pass)
        return 0;

    static bool blurp = false;
    if (!blurp) {
        errlogPrintf("FFTW Device Support %u.%u.%u%s (%s); using %s\n",
                     EPICS_FFTW_MAJOR_VERSION,
                     EPICS_FFTW_MINOR_VERSION,
                     EPICS_FFTW_MAINTENANCE_VERSION,
                     (EPICS_FFTW_DEVELOPMENT_FLAG ? "-dev" : ""),
                     "TODO:git hash",
                     EPICS_FFTW_LIBRARY_VERSION);
        blurp = true;
    }

    return 0;
}

template<typename REC>
long
init_record(REC *prec)
{
    try {
        dbCommon *pdbc = reinterpret_cast<dbCommon *>(prec);
        prec->dpvt = parseLink(pdbc, DBEntry(pdbc));
    }
    CATCH(__FUNCTION__)
    return 0;
}

template<typename REC>
long
init_record_write(REC *prec)
{
    long status = init_record<REC>(prec);
    if (!status)
        status = 2;
    return status;
}

template<typename REC>
long
init_record_read(REC *prec)
{
    return init_record<REC>(prec);
}

template<typename REC>
long
init_record_mask(REC *prec)
{
    long status = init_record_write<REC>(prec);
    if (prec->nobt == 0)
        prec->mask = 0xffffffff;
    prec->mask <<= prec->shft;
    return status;
}

template<typename REC>
long
init_record_write_arr(REC *prec)
{
    long status = init_record<REC>(prec);

    if (prec->ftvl != menuFtypeDOUBLE)
        throw std::runtime_error("Unsupported FTVL");

    return status;
}

template<typename REC>
long
init_record_read_arr(REC *prec)
{
    long status = init_record<REC>(prec);
    FFTWConnector *conn = static_cast<FFTWConnector *>(prec->dpvt);

    if (prec->ftvl != menuFtypeDOUBLE)
        throw std::runtime_error("Unsupported FTVL");

    if (prec->bptr) {
        free(prec->bptr); // get rid of record support allocated buffer
        conn->createEmptyOutputValue(&prec->bptr, prec->nelm);
    }
    return status;
}

long
get_iointr(int cmd, dbCommon *prec, IOSCANPVT *io)
{
    TRY { return conn->get_ioint(cmd, prec, io); }
    CATCH(__FUNCTION__)
}

template<typename REC>
long
write_enum(REC *prec)
{
    TRY
    {
        bool failed = true;
        if (conn->sigtype == FFTWConnector::SetWindowType) {
            switch (prec->rval) {
            case FFTWCalc::None:
            case FFTWCalc::Hann:
                failed = false;
                conn->setWindowType(static_cast<FFTWCalc::WindowType>(prec->rval));
                if (prec->tpro > 1)
                    std::cerr << prec->name << ": set window type " << conn->inst->fftw.wintype
                              << std::endl;
                break;
            }
        }
        if (!failed && conn->inst->triggerSrc == conn) {
            conn->setTimestamp(prec->time);
            conn->trigger();
        }
        if (failed) {
            (void) recGblSetSevr(prec, WRITE_ALARM, INVALID_ALARM);
            return S_dev_badRequest;
        }
        return 0;
    }
    CATCH(__FUNCTION__)
}

template<typename REC>
long
write_double(REC *prec)
{
    TRY
    {
        bool failed = true;
        if (conn->sigtype == FFTWConnector::SetSampleFreq) {
            failed = false;
            conn->setSampleFreq(analogEGU2Raw<double>(prec, prec->val));
            if (prec->tpro > 1)
                std::cerr << prec->name << ": set sample freq " << conn->getSampleFreq() << std::endl;
        }
        if (!failed && conn->inst->triggerSrc == conn) {
            conn->setTimestamp(prec->time);
            conn->trigger();
        }
        if (failed) {
            (void) recGblSetSevr(prec, WRITE_ALARM, INVALID_ALARM);
            return S_dev_badRequest;
        }
        return 0;
    }
    CATCH(__FUNCTION__)
}

template<typename REC>
long
write_double_arr(REC *prec)
{
    TRY
    {
        bool failed = true;
        if (conn->sigtype == FFTWConnector::InputReal) {
            if (prec->tpro > 1)
                std::cerr << prec->name << ": set input (real)" << std::endl;
            conn->setNextInputValue(prec->bptr, prec->nord);
            failed = false;
        }
        if (!failed && conn->inst->triggerSrc == conn) {
            conn->setTimestamp(prec->time);
            conn->trigger();
        }
        if (failed) {
            (void) recGblSetSevr(prec, WRITE_ALARM, INVALID_ALARM);
            return S_dev_badRequest;
        }
        return 0;
    }
    CATCH(__FUNCTION__)
}

template<typename REC>
long
read_double(REC *prec)
{
    TRY
    {
        long status = 0;
        bool failed = true;

        if (conn->sigtype == FFTWConnector::ExecutionTime) {
            double val = analogRaw2EGU<double>(prec, conn->getRuntime());
            prec->val = val;
            prec->udf = 0;
            prec->time = conn->getTimestamp();
            status = 2;
            failed = false;
        }
        if (failed) {
            (void) recGblSetSevr(prec, READ_ALARM, INVALID_ALARM);
            return S_dev_badRequest;
        }
        return status;
    }
    CATCH(__FUNCTION__)
}

template<typename REC>
long
read_double_arr(REC *prec)
{
    TRY
    {
        CHKVALID

        conn->getNextOutputValue(&prec->bptr, prec->nelm, &prec->nord);
        prec->udf = 0;
        prec->time = conn->getTimestamp();
        return 0;
    }
    CATCH(__FUNCTION__)
}


// Device Support jump tables

#define DEVSUPN(NAME, REC, IOINTR, OP, DIR) \
    static dset6<REC##Record> NAME = {6, \
                                      nullptr, \
                                      global_init, \
                                      init_record_##DIR<REC##Record>, \
                                      IOINTR, \
                                      &DIR##_##OP<REC##Record>, \
                                      nullptr};

#define DEVSUPI(NAME, REC, INIT, IOINTR, OP, DIR) \
    static dset6<REC##Record> NAME = {6, \
                                      nullptr, \
                                      global_init, \
                                      init_record_##INIT<REC##Record>, \
                                      IOINTR, \
                                      &DIR##_##OP<REC##Record>, \
                                      nullptr};

//      devsup name, record,      init, get_iointr, value type, direction
DEVSUPI(devMBBOfftw,   mbbo,      mask,    nullptr,       enum,     write)
DEVSUPN(  devAOfftw,     ao,            get_iointr,     double,     write)
DEVSUPN(  devAIfftw,     ai,            get_iointr,     double,      read)
DEVSUPI( devAAOfftw,    aao, write_arr, get_iointr, double_arr,     write)
DEVSUPI( devAAIfftw,    aai,  read_arr, get_iointr, double_arr,      read)

} // namespace

#include <epicsExport.h>

extern "C" {
    epicsExportAddress(dset, devMBBOfftw);
    epicsExportAddress(dset, devAOfftw);
    epicsExportAddress(dset, devAIfftw);
    epicsExportAddress(dset, devAAOfftw);
    epicsExportAddress(dset, devAAIfftw);
}
