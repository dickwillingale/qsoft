Module Source Code
******************
The modules are independent and have no extenal dependencies.
If a routine is required
by more than one module then a copy of the source code
is included in each of the relevant source directories.

The source directories are:

* $QSOFT/src/qfits
* $QSOFT/src/astro
* $QSOFT/src/images
* $QSOFT/src/xscat
* $QSOFT/src/xsrt

Each directory contains the Fortran interface routines, prefixed by QR, that
are called by Python, R or IDL and any internal Fortran subroutines required.

The source code for common blocks is held in files with uppercase
names without a Fortran extension \*.f, e.g. SRT\_COM, SPX\_COM.
These are dragged into the subroutine source using the Fortran INCLUDE
statement.

In addition the directories contain the module definition scripts, e.g.
images.py, images.R. The IDL function interface is provided by
individual \*.pro files for each command/function and these are held in
a subdirectory qIDL/ e.g. qIDL/qri\_getpos.pro, qIDL/qri\_beam.pro.

Each module directory contains a Makefile used to compile all the routines
and create shareable libraries for Python, R and IDL. Note that IDL
uses the library image created by R so R is required to produce the IDL
shareable. Once compiled you can install the library and definition
scripts using

$ make install

This copies the shareable libraries and definition scripts to the
directories $QSOFT/python\_modules, $QSOFT/R\_libraries and
$QSOFT/qIDL.

