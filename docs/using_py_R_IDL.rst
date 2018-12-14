Python, R and IDL
*****************

All the Fortran functions can be called from Python, R or IDL.
Because of peculiarities in the syntax and structure of the scripting languages
there are minor differences in the way the functions are accessed.

The documentation of all the functions uses the Python implementation. Where
there are significant differences in the R or IDL versions these are
mentioned in the text.

Python
======

The directory $QSOFT/python_modules is included in the PYTHONPATH at set up so
that the python  modules can be imported in the usual way. Here is a 
snippet of a Python script using the astro.cosmo() function

.. code-block:: python

    #!/usr/bin/env python
    # Test of Cosmological parameter calculations etc.
    import numpy as np
    import astro
    import matplotlib.pylab as plt
    #
    zmax=5
    # Einstein de Sitter
    c1=astro.cosmo(70,1,0,zmax)
    # Low density
    c2=astro.cosmo(70,0.05,0,zmax)
    # High Lambda
    c3=astro.cosmo(70,0.2,0.8,zmax)
    ...

RScript
=======

The file .Rprofile in the users home directory is executed by Rscript at start
up to dynamically load the shareable libraries. The QSOFT R function names
are prefixed according to the module library/subject as follows

* utilities: qr\_
* qfits: qr\_fits
* images: qri\_
* astro: qra\_
* xscat: qrx\_
* xsrt: qrt\_

Here is a snippet of a Rscript using the astro.cosmo() function

.. code-block:: R

    #!/usr/bin/env Rscript
    # Test of Cosmological parameter distance calculations
        zmax<-5
    # Einstein de Sitter
        c1<-qra_cosmo(70,1,0,zmax)
    # Low density
        c2<-qra_cosmo(70,0.05,0,zmax)
    # High Lambda
        c3<-qra_cosmo(70,0.2,0.8,zmax)
    ...

IDL
===
