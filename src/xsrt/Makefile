# Makefile for Python module xsrt
# RW 04-Jun-2016
# compiler definitions
include ../compiler.config
# Target directory for Python modules
PYMODDIR = $(QSOFT)/python_modules
RLIBDIR = $(QSOFT)/R_libraries
LIBDIR = $(QSOFT)/lib
# Python module
MODNAME = xsrt
# Internal Fortran object library
LIB = libxsrt.a
# Python Fortran module source files
FSOURCES = qr_init.f qrt_init.f qr_rseed.f \
qrt_aperture.f  qrt_kbs.f       qrt_moa.f       qrt_source.f    qrt_wolter2.f \
qrt_baffle.f    qrt_lens.f      qrt_opgrat.f    qrt_spider.f    \
qrt_c1nest.f    qrt_list.f      qrt_prism.f     qrt_sqpore.f    \
qrt_defs.f      qrt_mats.f      qrt_rotate.f    qrt_surface.f \
qrt_detector.f  qrt_mirror.f    qrt_shift.f     qrt_trace.f \
qrt_fresnel.f   qrt_sipore.f    qrt_w1nest.f \
qrt_bfield.f    qrt_bdipole.f   qrt_prtathena.f qrt_eltmxt.f \
qrt_elips.f	qrt_sqmpoarr.f
# Internal object library source files
LSOURCES = \
srt_cnic.f      srt_cnin.f      srt_dalr.f      srt_dfrm.f      srt_didi.f \
srt_in2d.f      srt_plne.f      srt_plna.f      srt_qufd.f \
srt_qufu.f      srt_gray.f      srt_qurt.f      srt_rcon.f      srt_refs.f \
srt_relk.f      srt_sphr.f      srt_mray.f      srt_trc.f       srt_vcrs.f \
srt_vdot.f      srt_vnrm.f      srt_wlt1.f      srt_pray.f      srt_nsur.f \
srt_su1.f       srt_su2.f       srt_su3.f       srt_su4.f       srt_su5.f \
srt_su6.f       srt_su7.f       srt_su8.f       srt_su9.f       srt_su10.f \
srt_su11.f      srt_su12.f      srt_su13.f      srt_su14.f      srt_su15.f \
srt_su16.f      srt_su17.f      srt_su18.f      srt_su19.f      srt_su20.f \
srt_su21.f      srt_su22.f      srt_su23.f      srt_su24.f      srt_su25.f \
srt_su26.f      srt_su27.f 	srt_su28.f \
srt_setf.f      srt_sett.f      srt_setd.f      srt_sets.f \
srt_mats.f      srt_cone.f      srt_fleg.f      srt_scrn.f      srt_list.f \
srt_tdef.f      srt_grad.f      srt_mcm.f       srt_mltr.f \
srt_famp.f      srt_acos.f      srt_cpii.f      srt_mlay.f      srt_setr.f \
srt_fpow.f      srt_fpor.f      srt_cpsd.f      srt_mlier.f     srt_liri.f \
srt_wlt2.f      srt_sph1.f      srt_save.f      srt_moa.f \
srt_shift.f     srt_rotate.f    srt_pore.f      srt_getpix.f    srt_nxtshft.f \
srt_sipore.f    srt_kbs.f       srt_kbmakeconstellation.f \
srt_bfield.f    srt_bdipole.f   srt_scat.f srt_frnl.f \
srt_trace.f     srt_kbfindmodule.f  \
srt_idfrm.f     srt_iin2d.f     srt_changeaxes.f  \
srt_setkbscom.f srt_findmod.f   srt_sfindmod.f srt_setconcom.f \
srt_sqmpoarr.f srt_findmpo.f srt_setmpoarr.f \
sys_gauss.f sys_srand.f sys_drand.f sys_dran0.f  sys_getlun.f 

CSOURCES = stray_xray_psf.c
# Routines to skip in Python Fortran module
SKIP = skip: qr_init escale arcsinh scatfun diffsec :

COBJECTS= $(CSOURCES:.c=.o)
OBJECTS = $(LSOURCES:.f=.o)
FOBJECTS= $(FSOURCES:.f=.o)

shareables:  $(FSOURCES) $(LIB) $(FOBJECTS)
ifdef F2PY
	$(F2PY) -c -m $(MODNAME)for $(FSOURCES) $(SKIP) $(LIB)
endif
ifdef R
	$(R) CMD SHLIB -o $(MODNAME)R.so $(FOBJECTS) $(LIB) 
endif

$(LIB): $(OBJECTS) $(COBJECTS)
	ar rsu $(LIB) $(OBJECTS) $(COBJECTS)

.SUFFIXES: .o .f .c
.f.o:
	$(F77) $(FFLAGS) -c $<
.c.o:
	$(CC) $(CFLAGS) $(LARGEFLAGS) -c $<


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
