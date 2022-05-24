# Operation

The central entity that runs the FFT for one set of inputs and outputs
is called an *instance*.

The job of calculating the transformation is handed to a set of
background worker threads within an `epicsThreadPool`.
This allows for smooth and efficient execution of any number of FFT
instances on a single IOC, without blocking EPICS database execution.

FFT instances are identified by a unique name.
The EPICS database connects input and output arrays, configuration
parameters and statistics records to a specific instance, using
the instance name in the INP/OUT link.
A *connector* item is used for each record to link to an FFT instance,
keeping the record specific configuration and data.

One of the (input) records is configured to trigger the
transformation.
When this record processes, the connected FFT instance is triggered,
which adds its `calculate()` method to the EPICS job queue to be picked
up by a worker thread.
At the end of the `calculate()` method, the instance is pushing new
data to the output and stats connectors and registers the connected
records for processing.

C++ 11 `shared_ptr<>` and `unique_ptr<>` are used to avoid unnecessary
copying of arrays. Multiple output records presenting different parts
of the same FFT instance's output data will all use the same shared
array, with their BPTR fields possibly pointing to different sections
of it. (Allowing for low overhead Region-of-Interest records.)

## Code Overview

### fftwCalc

Thin wrapper around the FFTW library functions. It's the only class
that needs access to FFTW headers.

### fftwInstance

Instance of the transformation. Keeps lists of input and output
connectors and contains the main `calculate()` method.

### fftwConnector

Connects one EPICS record to an FFTW instance, keeping the record's
signal type and parameters and handling the `shared_ptr<>` mechanics
and BPTR swapping.

### fftwSupport

C++ implementation of the Device Support functions.

### iocshIntegration

iocShell integration (sic!).

## iocShell interface

### fftwShow - Print Diagnostic and Stats

Called with an instance name and a verbosity level.

Prints a report on the status and configuration of the instance,
lists the connected records and their signals, vector sizes etc.
