# Makefile for Python module astro
# RW 05-Jun-2016
# compiler definitions
include ../compiler.config
# Target directory for Python modules
PYMODDIR = $(QSOFT)/python_modules
RLIBDIR = $(QSOFT)/R_libraries
LIBDIR = $(QSOFT)/lib
# Python module
MODNAME = astro
# Internal Fortran object library
LIB = libastro.a
# Python Fortran module source files
FSOURCES = qrx_init.f qr_init.f \
qra_habs.f qra_cosmo.f \
qrx_ismtau.f     qrx_iismtau.f    qrx_igmtau.f     qrx_iigmtau.f \
qrx_setabund.f   qrx_lyftau.f     qrx_lyftauvz.f   qrx_iigmtauvz.f \
qra_kcorrb.f qra_bandint.f   qra_band.f qrx_brems.f

# Internal object library source files
LSOURCES = \
spf_bbody.f spf_brem.f spf_feabs.f spf_hyabs.f spf_gaunt.f \
spf_hytau.f spf_plaw.f spx_bbody.f spx_brems.f spx_embhyabs.f \
spx_feabs.f spx_hyabs.f spx_hyabspat.f spx_plaw.f spx_sa.f \
spf_band.f spx_band.f spf_bplaw.f spx_bplaw.f \
xsphab.f xszphb.f xszvab.f xszvph.f gphoto.f photo.f fgabn.f fgxsct.f \
setxsec.f phfit2.f xsabsori_tau.f xspecinit.f \
len_trim.f sys_getenv.f sys_getlun.f

CSOURCES = ismatten.c ismatten_.c

# Routines to skip in Python Fortran module
SKIP = skip: qr_init escale arcsinh scatfun diffsec :

OBJECTS = $(LSOURCES:.f=.o) $(CSOURCES:.c=.o)
FOBJECTS= $(FSOURCES:.f=.o)

shareables: $(FSOURCES) $(LIB) $(FOBJECTS)
ifdef F2PY
	$(F2PY) -c -m $(MODNAME)for $(FSOURCES) $(SKIP) $(LIB)
endif
ifdef R
	$(R) CMD SHLIB -o $(MODNAME)R.so $(FOBJECTS) $(LIB)
endif
	        
$(LIB): $(OBJECTS)
	ar rsu $(LIB) $(OBJECTS)

.SUFFIXES: .o .f .c
.f.o:
	$(F77) $(FFLAGS) -c $<

.c.o:
	$(CC) $(CFLAGS) -c $<

install: $(MODNAME).R $(MODNAME).py
	cp $(MODNAME)for*.so $(PYMODDIR)
	cp $(MODNAME).py $(PYMODDIR)
	cp $(MODNAME)R.so $(RLIBDIR)
	cp $(MODNAME).R $(RLIBDIR)
	cp $(LIB) $(LIBDIR)

clean:
	rm $(OBJECTS)
	rm *.so
	rm *.a
