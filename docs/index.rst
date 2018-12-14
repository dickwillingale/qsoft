.. QSOFT documentation master file, created by
   sphinx-quickstart on Thu Nov  1 08:24:59 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to QSOFT
****************

QSOFT is a collection of data analysis and modelling applications for use in X-ray astronomy and related disciplines.

The applications are run as commands/functions within Python, R or IDL.

The core code is written in Fortran and C, compiled to produce shareable object libraries and imported as modules in Python or loaded by R or IDL.

The commands/functions are defined in Python modules, R or IDL scripts.

The QSOFT collection should be built using the gcc and gfortran compilers. The Python module f2py and/or R must be available to create the shareable objects.

.. toctree::
   :maxdepth: 2

   installation
   using_py_R_IDL
   qfits
   images
   astro
   xscat
   xsrt
   qsoft_modifying
   authors
