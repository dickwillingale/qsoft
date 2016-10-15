# Makefile for Python module qfits
# RW 14-Jul-2016
# compiler definitions
include ../compiler.config
# Target directory for Python modules
PYMODDIR = $(QSOFT)/python_modules
RLIBDIR = $(QSOFT)/R_libraries
LIBDIR = $(QSOFT)/lib
# Python modules
TARGET = qfits
# Python Fortran/C object library
LIB = libqfits.a
# Source files for module
FSOURCES = qrf_init.f \
qr_fitsclose.f   qr_fitsgcvm.f    qr_fitsparrd.f   qr_fitspcols.f \
qr_fitscolnam.f  qr_fitsgcx.f     qr_fitsparrj.f   qr_fitspcom.f \
qr_fitsempty.f   qr_fitsgetkey.f  qr_fitspbtab.f   qr_fitsphis.f \
qr_fitsgcvl.f    qr_fitsgpvd.f    qr_fitspcolb.f   qr_fitspkeyd.f \
qr_fitsgcs.f     qr_fitsgpvj.f    qr_fitspcolc.f   qr_fitspkeyj.f \
qr_fitsgcvc.f    qr_fitshdu.f     qr_fitspcold.f   qr_fitspkeyl.f \
qr_fitsgcvd.f    qr_fitsnew.f     qr_fitspcolj.f   qr_fitspkeys.f \
qr_fitsgcvj.f    qr_fitsopen.f    qr_fitspcoll.f   qr_fitstypes.f \
qr_fitsexists.f  qr_fitscolunit.f qr_fitspcolx.f   qr_fitspcolm.f
# C source files for internal library
CSOURCES =buffers.c cfileio.c checksum.c drvrfile.c drvrmem.c \
editcol.c edithdu.c eval_l.c \
eval_y.c eval_f.c fitscore.c getcol.c getcolb.c getcold.c getcole.c \
getcoli.c getcolj.c getcolk.c getcoll.c getcols.c getcolsb.c \
getcoluk.c getcolui.c getcoluj.c getkey.c group.c grparser.c \
histo.c iraffits.c \
modkey.c putcol.c putcolb.c putcold.c putcole.c putcoli.c \
putcolj.c putcolk.c putcoluk.c putcoll.c putcols.c putcolsb.c \
putcolu.c putcolui.c putcoluj.c putkey.c region.c scalnull.c \
swapproc.c wcssub.c wcsutil.c imcompress.c quantize.c ricecomp.c \
pliocomp.c fits_hcompress.c fits_hdecompress.c zuncompress.c zcompress.c \
adler32.c crc32.c inffast.c inftrees.c trees.c zutil.c \
deflate.c infback.c inflate.c uncompr.c \
f77_wrap1.c f77_wrap2.c f77_wrap3.c f77_wrap4.c
# Fortran source files for internal library
LSOURCES = len_trim.f sys_getlun.f

SKIP = skip: qr_fitsexists :

COBJECTS= $(CSOURCES:.c=.o)
OBJECTS = $(FSOURCES:.f=.o)
LOBJECTS = $(LSOURCES:.f=.o)

$(TARGET)for.so: $(FSOURCES) $(LIB) $(OBJECTS)
ifdef F2PY
	$(F2PY) -c -m $(TARGET)for $(FSOURCES) $(SKIP) $(LIB)
endif
ifdef R
	$(R) CMD SHLIB -o $(TARGET)R.so $(OBJECTS) $(LIB)
endif

$(LIB): $(COBJECTS) $(LOBJECTS)
	ar rsu $(LIB) $(COBJECTS) $(LOBJECTS)

.SUFFIXES: .o .f
.f.o:
	$(F77) $(FFLAGS) -c $<
.c.o:
	$(CC) $(CFLAGS) $(LARGEFLAGS) -c $<

install: $(TARGET).py $(TARGET).R
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
