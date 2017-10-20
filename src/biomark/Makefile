# Makefile for Python/R module biomark
# RW 26-Jun-2017
# compiler definitions
include ../compiler.config
# Target directory for Python/R modules
PYMODDIR = $(QSOFT)/python_modules
RLIBDIR = $(QSOFT)/R_libraries
LIBDIR = $(QSOFT)/lib
# Python module
MODNAME = biomark
# Internal Fortran object library
LIB = libbiomark.a
# Python Fortran module source files
FSOURCES = qrb_init.f \
qrb_biomapping.f  qrb_bioroc.f      qrb_biosort.f \
qrb_bioback.f     qrb_biomarkers.f  qrb_bioscale.f    qrb_biospread.f \
qrb_biobin.f      qrb_biomatch.f    qrb_bioscat.f     qrb_biotest.f \
qrb_bioclass.f    qrb_bionorm.f     qrb_bioscore.f \
qrb_biohunt.f     qrb_biopeaks.f    qrb_biosigpks.f
# Internal object library source files
LSOURCES = \
pda_dsort.f  pda_xermsg.f sx_chisqd.f sx_pchisq.f sx_gaussd.f \
sx_gamma.f sx_factor.f

# Routines to skip in Python Fortran module
SKIP = skip: qr_init escale arcsinh scatfun diffsec :

OBJECTS = $(LSOURCES:.f=.o)
FOBJECTS= $(FSOURCES:.f=.o)

shareables:  $(FSOURCES) $(LIB) $(FOBJECTS)
ifdef F2PY
	$(F2PY) -c -m $(MODNAME)for $(FSOURCES) $(SKIP) $(LIB)
endif
ifdef R
	$(R) CMD SHLIB -o $(MODNAME)R.so $(FOBJECTS) $(LIB) 
endif

$(LIB): $(OBJECTS)
	ar rsu $(LIB) $(OBJECTS)

.SUFFIXES: .o .f
.f.o:
	$(F77) $(FFLAGS) -c $<


install: $(MODNAME).R $(MODNAME).py
ifdef F2PY
	cp $(MODNAME)for*.so $(PYMODDIR)
	cp $(MODNAME).py $(PYMODDIR)
endif
ifdef R
	cp $(MODNAME)R.so $(RLIBDIR)
	cp $(MODNAME).R $(RLIBDIR)
endif
	cp $(LIB) $(LIBDIR)

clean:
	rm $(OBJECTS)
	rm *.so
	rm *.a
