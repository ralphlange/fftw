<a target="_blank" href="http://semver.org">![Version][badge.version]</a>
<a target="_blank" href="https://www.codacy.com/gh/ralphlange/fftw">![Codacy grade][badge.codacy]</a>

# fftw - EPICS Device Support for the FFTW Library

C++ [EPICS](https://epics-controls.org) Device Support module wrapping around
the [FFTW](https://www.fftw.org/) library v3.

The Device Support uses C++ shared pointers and the feature to swap the BPTR
of EPICS array records to minimize copying of arrays. 

## Status

One of the DFT algorithm variants provided by the FFTW library is supported:
One-dimensional, real input data transformed from the time domain (sampled data)
to the frequency domain (amplitude/phase).

The module can be extended to support other variants.

## Prerequisites

*   A C++ compiler that supports the C++11 standard. \
    Microsoft Visual C++ needs to be from Visual Studio 2015 or newer.
    g++ needs to be 4.6 or above.

*   [EPICS Base](https://epics-controls.org/resources-and-support/base/)
    release 3.15 (and up; EPICS 7 is supported).

*   The [pyepics module](https://pyepics.github.io/pyepics) for running
    the tests.

*   The [FFTW](https://www.fftw.org/) library (v3).

## Building the module

This is a standard EPICS module.

Inside the `configure` subdirectory or one level above the TOP location
(TOP is where this README file resides), create a file `RELEASE.local`
that sets `EPICS_BASE` to the absolute path inside your EPICS
installation.

Configure the compiler on Linux to use the C++11 standard by adding
```makefile
USR_CXXFLAGS_Linux += -std=c++11
```
to the `CONFIG_SITE` file (or one of the host/target specific site
configuration files). \
It is preferable to set this option globally in EPICS Base.

## Using the module

IOC applications that use the module need to

*   add an entry to the Device Support module in their `RELEASE.local` file
*   include `fftw.dbd` when building the IOC's DBD file
*   include `fftwSup` in the support libraries for the IOC binary.

## Documentation

Sparse.

## Feedback / Reporting issues

Please use the GitHub project's
[issue tracker](https://github.com/ralphlange/fftw/issues).

## Credits

This module is based on ideas and code snippets from
Michael Davidsaver (Osprey DCS).

## License

This module is distributed subject to a Software License Agreement found
in file LICENSE that is included with this distribution.

<!-- Links -->
[badge.version]: https://img.shields.io/github/v/release/ralphlange/fftw?sort=semver
[badge.codacy]: https://app.codacy.com/project/badge/Grade/5b02a83f024f4ef7ab725d699ac8fb3a
