# Makefile for Python module xscat
# RW 05-Jun-2016
# compiler definitions
include ../compiler.config
# Target directory for Python modules
PYMODDIR = $(QSOFT)/python_modules
RLIBDIR = $(QSOFT)/R_libraries
LIBDIR = $(QSOFT)/lib
# Python module
MODNAME = xscat
# Internal Fortran object library
LIB = libxscat.a
# Python Fortran module source files
FSOURCES = qrs_init.f miev0_py.f ErrPack.f \
qra_dust.f qrt_xopt.f qrt_xfresnel.f qrt_mlayer.f

# Internal object library source files
LSOURCES = \
srt_cpii.f  srt_frnl.f  srt_mcm.f   srt_mlay.f  srt_mlti.f \
srt_famp.f      srt_acos.f \
xx_cromer.f       xx_get_hagemann.f xx_lgndr.f        xx_sigma1.f \
xx_decode.f       xx_get_henke.f    xx_optical.f      xx_sigma2.f \
xx_dielec.f       xx_hagemann.f     xx_optl.f         xx_zref.f \
xx_gauss.f        xx_henke.f        xx_ref.f \
xx_get_cromer.f   xx_href.f         xx_sigma0.f \
sys_getenv.f sys_getlun.f sys_upcase.f scan.f len_trim.f

# Routines to skip in Python Fortran module
SKIP = skip: LPCOEF LPCO1T LPCO2T BIGA MIPRNT SMALL1 SMALL2 TESTMI ErrMsg :

OBJECTS = $(LSOURCES:.f=.o)
FOBJECTS = $(FSOURCES:.f=.o)

$(MODNAME)for.so:  $(FSOURCES) $(LIB) $(FOBJECTS)
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

install: $(MODNAME).py $(MODNAME).R
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
