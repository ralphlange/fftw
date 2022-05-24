# Database Interface

DTYP = "FFTW".

All record configurations use the instance name (macro `P` in the
test/sample database) as first token in the INP/OUT link.

## Parameters

### windowtype

Window function to apply to the input data.
Used with an mbbo record.
*   0 = no window function
*   1 = Hanning window

### sample-freq

Sampling frequency of the input data [Hz].
Used with an ao record.

## Inputs

One of the defined input records can set a link option
"trigger=y". This record will trigger the tranformation.

### input-real

Real part of the input data.
Used with an aao record of type DOUBLE.

### input-real using aSub

Fetching the real part of the input data from a different array
record.
Using an aSub record.
The record needs to set INAM = "FFTW_init", SNAM = "FFTW_input",
FTA = "DOUBLE", NOA = <size of input array>, and an info item with
the usual configuration (<instance name> input-real). 

## Outputs

The maximum used size of the output arrays is
half of the input size plus one.
(Except for the output-window signal that uses the original size
of the input array.)

Output records can set the link option "skipDC=y" to not include the
first value. They can set "offset=<n>" to start a an arbitrary offset
into the output array. ("skipDC=y" is equivalent to "offset=1".)

The records of outputs are driven by the instance and must be set
to SCAN = "I/O Intr".

### output-real

Real part of the output data.
Used with an aai record of type DOUBLE.

### output-imag

Imaginary part of the output data.
Used with an aai record of type DOUBLE.

### output-magn

Magnitude of the output data.
Used with an aai record of type DOUBLE.

### output-phas

Phase of the output data.
Used with an aai record of type DOUBLE.

### output-fscale

Frequency scale of the output data. Contains the x-axis values
[Hz] for plotting the output arrays.
Used with an aai record of type DOUBLE.
Frequency scales matching output records with skipDC or offset
options can be obtained by using the same link option.

### output-window

Window function used on the input data.
Used with an aai record of type DOUBLE.
The maximum used size of the output-window array is
the size of the input.

### exectime

Execution time of the last transformation [s].
Used with an ai record.
