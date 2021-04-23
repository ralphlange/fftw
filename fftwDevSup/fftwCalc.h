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

#ifndef FFTWCALC_H
#define FFTWCALC_H

#include <cstdlib>
#include <ctime>
#include <vector>
#include <memory>

#include <fftw3.h>

#include <errlog.h>

extern int FFTWDebug;

// STL compatible allocator which uses fftw_alloc_*() to ensure aligned arrays
template<typename T>
class FFTWAllocator : public std::allocator<T>
{
public:
#define PTYPE(P) typedef typename std::allocator<T>::P P
    PTYPE(pointer);
    PTYPE(const_pointer);
    PTYPE(reference);
    PTYPE(const_reference);
    PTYPE(size_type);
#undef PTYPE

    inline FFTWAllocator() {}

    template<typename U> struct rebind {typedef FFTWAllocator<U> other;};

    inline pointer address(reference x) const { return &x; }
    inline const_pointer address(const_reference x) const { return &x; }
    inline size_type max_size()const {return ((size_t)-1)/sizeof(T);}

    inline void construct(pointer p, const_reference val);

    inline void destroy(pointer p)
    {p->~T();}

    inline pointer allocate(size_type n, const void * =0)
    {
        return (T*)fftw_malloc(n*sizeof(T));
    }

    inline void deallocate(pointer p, size_type n)
    {
        fftw_free(p);
    }
};

template<>
inline void FFTWAllocator<double>::construct(pointer p, const_reference val)
{::new((void*)p) double(val);}

template<>
inline void FFTWAllocator<fftw_complex>::construct(pointer p, const_reference val)
{
    (*p)[0] = val[0];
    (*p)[1] = val[1];
}

template<>
inline void FFTWAllocator<double>::destroy(pointer p)
{}

template<>
inline void FFTWAllocator<fftw_complex>::destroy(pointer p)
{}

// Helper to ensure that plans are destroyed
class Plan
{
public:
    Plan() :plan(nullptr) {}
    ~Plan() {clear();}
    void clear()
    {
        if(plan)
            fftw_destroy_plan(plan);
        plan = nullptr;
    }
    void reset(fftw_plan p=nullptr)
    {
        if(!p)
            throw std::bad_alloc();
        if(plan)
            fftw_destroy_plan(p);
        plan = p;
    }
    Plan& operator=(fftw_plan p) {reset(p); return *this;}
    fftw_plan get() const {return plan;}
private:
    fftw_plan plan;
};

struct FFTWCalc
{
    enum WindowType {
        None = 0,
        Hann,
    };

    static inline const char *
    WindowTypeName(const WindowType s)
    {
        switch (s) {
        case None:
            return "None";
        case Hann:
            return "Hann";
        default:
            break;
        }
        return "?";
    }

    WindowType wintype;
    std::vector<double> window;

    std::unique_ptr<std::vector<double, FFTWAllocator<double>>> input;
    std::vector<fftw_complex, FFTWAllocator<fftw_complex>> output;

    size_t input_sz;
    size_t ntime, nfreq;

    Plan plan;

    double fsamp;
    std::vector<double> fscale;

    bool redo_plan, newval;

    FFTWCalc();
    ~FFTWCalc();

    bool set_fsamp(double f);
    bool set_wtype(FFTWCalc::WindowType type);

    void set_input(std::unique_ptr<std::vector<double, FFTWAllocator<double>>> inp);
    bool apply_window();
    bool replan();
    void transform();
};

#endif // FFTWCALC_H
