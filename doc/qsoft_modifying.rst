Adding New Functions
********************

The top-level subroutines are written in Fortran 77.

Subroutines names are prefixed with:

* QRA\_ astronomy/astrophysics
* QRX\_ X-ray astronomy/physics
* QRI\_ image processing
* QR\_FITS for the FITS file interface
* QRT\_ sequential rays tracing interface
* QR\_ general utilities

Internal Fortran routines also have prefixes:

* SRT\_ internal sequential ray tracing
* AX\_ coordinate transformations
* XX\_ X-ray physics
* SYS\_ system utilities
  
Two general Fortran routines are also suppied, SCAN and LEN\_TRIM.

A few internal routines are written in C or C++. These are best left alone!

All the source code is held in directory:

$QSOFT/src/

.. toctree::

   qsoft_source_code
   xsrt_ray_tracing
   xsrt_new_surfaces
