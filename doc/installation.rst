Build and Installation
**********************
The following software items are required for the build:

* Fortran compiler - GNU gfortran used in development
* C compiler - GNU gcc used in development
* Python module f2py - to build the modules for Python
* R - to build the shareable libraries for R and IDL
* Python module sphinx - to build the documentation

1. Download from GitHub

        $ git clone git://github.com/dickwillingale/qsoft.git

        This will create a directory qsoft/

2) Move into /your/files/top/qsoft/src

        $ cd /your/files/top/qsoft/src

        Edit the compile.config so that the compilers CC and F77 and R, F2PY and IDL
        point to the correct executables on your system.
        If you don't have R or F2PY leave them blank. If your target is IDL you will
        need R to compile the shareable library.

3. Move into /your/files/top/qsoft and execute build

        $ cd /your/files/top/qsoft

        $ ./build

        This will check you have gcc gfortran and R and/or Python with f2py.  It will then create the src/compiler.config file for make and compile the shareable objects.

4. Put the following line into your /home/.profile

        . /your/files/top/qsoft/setup_q

        Qsoft will be available when you launch a login terminal.

5. That's it.

        When you start R (or Rscript) it will automatically load the QSOFT applications
        using the /home/.Rprofile file.

        The environment variable PYTHONPATH will point to the qsoft/python_modules
        directory so you can load the modules into Python.
