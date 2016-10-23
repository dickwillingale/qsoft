#!/usr/bin/env python3
# Test of the sequential ray tracing routines
from __future__ import print_function
import numpy as np
import xsrt
import xscat
import images
import matplotlib.pylab as plt
import pandas as pd
# Set up SPO modules
from athena_ladapt_modules import *
# Set up ray tracing of finite source distance testing
panlen= 1000000
# useful vectors
sn=np.array([1,0,0])
nn=np.array([0,0,0])
rx=np.array([0,1,0])
# reflecting surface index
isu=1
print("energy",ekev,"keV")
if ist==1:
    print("xopt alpha gamma",alpha,gamma)
else:
    print("reflectivity look-up table")
    na=angs.size
    for i in range(angs.size):
        print(angs[i],refs[i])
print("surface roughness parameters",nrough,rind,brk)
# Set position of aperture centre
a2j=max(lm)
pcen=np.array([flen+a2j,0,0])
dshift= panlen*flen/(panlen-flen)-flen
flenprime= flen+dshift
print("detector shift for Panter finite distance",dshift)
# set source to cover aperture
ssize=rmax*1.05
smin=rmin*0.95
slim=np.array([smin,ssize,0,0,0,0])
spos=np.array([flen+300,0,0])
dsrc=np.array([flen+panlen,0,0])
nray=50000
# Detector
drad=ssize
dpos=np.array([-dshift,0,0])
dlim=np.array([0,drad,0,0])
radet=0.0
# set deformation index
ide=1
# set image parameters - pixel 0.25 arc secs
dpix= (0.25/3600)*np.pi/180*flen
irad=4
nx= np.int(irad/dpix*2)
ny= nx
irad= nx*dpix/2
print("image parameters",dpix,nx,ny,irad)
#
print("arcmin area  hew w90 fwhm dshift")
#
offarcmin=0.0
offrad=(offarcmin/60)*np.pi/180.
ccy=np.sin(offrad)
di=np.array([-np.sqrt(1.0-ccy**2),ccy,0.0])
print("di",di)
xsrt.reset()
xsrt.rseed(55)
xsrt.surface(isu,ist,ekev,nrough,brk,rind,alpha,gamma,angs,refs,0,0,0)
xsrt.deform(ide,1,5,nm,1)
xsrt.defmat(ide,1,xd,yd,dx)
xsrt.defmat(ide,2,xd,yd,dy)
xsrt.defmat(ide,3,xd,yd,dz)
xsrt.defmat(ide,4,xd,yd,db)
xsrt.defmat(ide,5,xd,yd,dc)
xsrt.sipore(pcen,sn,rx,flen,rpitch,apitch,wall,rm,pm,te,wm,hm,lm,cm,gm,wfr,
a2j,ide,isu)
xsrt.source(3,di,dsrc,spos,sn,rx,slim,0,nray,0)
xsrt.detector(1,dpos,sn,rx,dlim,radet)
results=xsrt.trace(0,drad,2)
dshft=results.dshft
#
detpos=pd.read_table("detected.dat",sep="\s+")
rad= flen*max(gr)*offrad/8
xcen= offrad*flen+rad
ycen= 0
xmin= xcen-irad
xmax= xcen+irad
ymin= ycen-irad
ymax= ycen+irad
arr=images.binxy(detpos["YDET"],detpos["ZDET"],0,detpos["AREA"],
xmin,xmax,ymin,ymax,nx,ny)
xsam=(xmax-xmin)/nx
ysam=(ymax-ymin)/ny
print("sum of array",np.sum(arr))
print("xsam",xsam)
# 
plt.figure()
cm=plt.get_cmap("gnuplot2")
cm=plt.get_cmap("hot")
plt.imshow(arr,extent=[xmin,xmax,ymin,ymax],origin="lower",cmap=cm)
images.setfield(nx,xmin,xmax,ny,ymin,ymax)
images.setpos(2,[0,0])
bb=images.beam(arr,nx,0.0,0.0)
print("peak cen",bb.peak,bb.cen)
area=bb.flux/100.
hew=(bb.hew*xsam/flen)*(180/np.pi)*3600
w90=(bb.w90*xsam/flen)*(180/np.pi)*3600
fwhm=(bb.fwhm*xsam/flen)*(180/np.pi)*3600
print(area,hew,w90,fwhm,dshft)
#a$xlab= "mm"
#a$ylab= "mm"
#qri_displayimage(a)
#go_on= locator(1)
#plot(ras,sb$y,log="y",type="l",
#main="CDF aperture 15 mm short",
#xlab="radius arc secs",ylab="surface brightness")
#
#cum= cumsum(sb$y)
#cum= cum/max(cum)
#plot(ras,cum,type="l",xlab="radius arc secs",ylab="flux fraction")
#go_on= locator(1)
#dev.off()
plt.show()
