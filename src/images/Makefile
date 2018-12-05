# Makefile for Python module images
# RW 06-Jun-2016
# compiler definitions
include ../compiler.config
# Target directories for Python, R and IDL modules
PYMODDIR = $(QSOFT)/python_modules
RLIBDIR = $(QSOFT)/R_libraries
IDLDIR = $(QSOFT)/qIDL
LIBDIR = $(QSOFT)/lib
# Python module
MODNAME = images
# Internal Fortran object library
LIB = libimages.a
# Python Fortran module source files
FSOURCES = qr_init.f qri_init.f qr_xyztrans.f \
qri_binxy.f qri_setfield.f qri_setsky.f qri_setpos.f qri_getpos.f \
qri_beam.f qri_annulus.f qri_blur.f qri_lecbeam.f qri_lecimage.f \
qri_lebin.f qri_lepsf.f qri_sqbeam.f

# Internal object library source files
LSOURCES = \
ax_barcor.f     ax_dammerinv.f  ax_donvcp.f     ax_locpos.f \
ax_barpre.f     ax_dayhms.f     ax_donvrt.f     ax_mjddat.f \
ax_barvel.f     ax_dec2dms.f    ax_dtcor.f      ax_mjdgmst.f \
ax_boxsiz.f     ax_deg2dms.f    ax_dtrans.f     ax_mjdname.f \
ax_cona2v.f     ax_degdms.f     ax_ecliptic.f   ax_mjdutd.f \
ax_congal.f     ax_dgetgal.f    ax_galcon.f     ax_mjdymd.f \
ax_congen.f     ax_dmat.f       ax_gendmat.f    ax_ra2hms.f \
ax_conmin.f     ax_dona2v.f     ax_getgal.f     ax_sunrad.f \
ax_conmxm.f     ax_dongal.f     ax_gorbit.f     ax_tcor.f \
ax_conv2a.f     ax_dongen.f     ax_gtcirc.f     ax_trisiz.f \
ax_convcp.f     ax_donmin.f     ax_iauname.f    ax_utdymd.f \
ax_convrt.f     ax_donmxm.f     ax_lambert.f    ax_ymdmjd.f \
ax_dammer.f     ax_donv2a.f     ax_lambertinv.f

# Routines to skip in Python Fortran module
SKIP = skip: LPCOEF LPCO1T LPCO2T BIGA MIPRNT SMALL1 SMALL2 TESTMI ErrMsg :

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
	cp locator.py $(PYMODDIR)
endif   
ifdef R 
	cp $(MODNAME)R.so $(RLIBDIR)
	cp $(MODNAME).R $(RLIBDIR)
endif   
ifdef IDL
	cp qIDL/*.pro $(IDLDIR)/$(MODNAME)
endif
	cp $(LIB) $(LIBDIR)

clean:
	rm $(OBJECTS)
	rm *.so
	rm *.a
