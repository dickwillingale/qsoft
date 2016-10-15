# Image processing routines
import imagesfor 
# Initialisation and reseting
def init():
    imagesfor.qri_init()
def reset():
    imagesfor.qri_init()
init()
def binxy(x,y,iq,w,xleft,xright,ybot,ytop,nx,ny):
    a=imagesfor.qri_binxy(x,y,iq,w,xleft,xright,ybot,ytop,nx,ny)
    return a
class bdata: pass
def beam(arr,rbeam,blev,bvar):
    a=imagesfor.qri_beam(arr,rbeam,blev,bvar)
    b=bdata()
    b.nsam=a[0]
    b.bflux=a[1]
    b.bsigma=a[2]
    b.flux=a[3]
    b.fsigma=a[4]
    b.peak=a[5]
    b.cen=a[6]
    b.tha=a[7]
    b.rmsa=a[8]
    b.rmsb=a[9]
    b.fwhm=a[10]
    b.hew=a[11]
    b.w90=a[12]
    b.fwhmp=a[13]
    b.hewp=a[14]
    b.w90p=a[15]
    b.fwhmc=a[16]
    b.hewc=a[17]
    b.w90c=a[18]
    return b
def setfield(nx,xleft,xright,ny,ybot,ytop):
    imagesfor.qri_setfield(nx,xleft,xright,ny,ybot,ytop)
def setsky(xtodeg,ytodeg,ipr,mjd,ra,dec,roll):
    imagesfor.qri_setsky(xtodeg,ytodeg,ipr,mjd,ra,dec,roll)
def setpos(ipos,p):
    imagesfor.qri_setpos(ipos,p)
class iposition: pass
def getpos():
    a=imagesfor.qri_getpos()
    b=iposition()
    b.pix=a[0]
    b.xyl=a[1]
    b.aes=a[2]
    b.equ=a[3]
    b.ecl=a[4]
    b.gal=a[5]
    return b
