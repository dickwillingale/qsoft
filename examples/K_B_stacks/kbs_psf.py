#!/usr/bin/env python
# Point spread function of Kirpatrick-Baez stacks
from __future__ import print_function
import numpy as np
import xsrt
import xscat
import images
import matplotlib.pylab as plt
import pandas as pd
# useful vectors
sn=np.array([1,0,0])
nn=np.array([0,0,0])
rx=np.array([0,1,0])
# reflecting surface definition
rough=0
rind=1.4
brk=1.0
ekev=1.0
iss=1
Ir=xscat.xopt("Ir",22.65,ekev,1)
alpha=Ir.alpha
gamma=Ir.gamma
angs=np.array([0])
refs=np.array([0])
# Deformations
idd=0
# K-B stack parameters
fl=5000
#ipack=1
#rmin=244
ipack=3
rmin=0
rmax=1900
pnor=np.array([1,0,0])
raxi=np.array([0,1,0])
pdd=0.610
wall=0.15
pitch=pdd+wall
csize=(rmax-rmin)/10
grdeg=1.0
gr=grdeg*np.pi/180
pl=pdd/gr
print("grazing angle",grdeg,"degrees, axial length",pl,"mm")
pcen1=np.array([fl+pl,0,0])
pcen2=np.array([fl,0,0])
# set off-axis direction
offarcmin=-2.0
offrad=(offarcmin/60)*np.pi/180.
ccy=np.sin(offrad)
di=np.array([-np.sqrt(1.0-ccy**2),-ccy,0.0])
rlim=np.array([rmin-10,rmax+10,0,0,0,0])
spos=np.array([fl+pl+3,0,0])
nray=1000000
# Detector - spherical with radius of curvature equal focal length
#dnorm=-di
dnorm=np.array([1,0,0])
dpos=np.array([-fl,0,0])+dnorm*fl
dlim=np.array([0,rmax,0,0])
radet=fl
halfd=4
#
xsrt.reset()
s1=xsrt.kbs(pcen1,pnor,raxi,ipack,rmin,rmax,fl,csize,pitch,wall,pl,pl,idd,iss)
s2=xsrt.kbs(pcen2,pnor,raxi,ipack,rmin,rmax,-fl,csize,pitch,wall,pl,pl,idd,iss)
xsrt.source(1,di,nn,spos,sn,rx,rlim,0,nray,0)
xsrt.detector(3,dpos,dnorm,rx,dlim,radet)
xsrt.surface(iss,1,ekev,rough,brk,rind,alpha,gamma,angs,refs,0,0,0)
xsrt.srtlist()
results=xsrt.trace(0,100,2)
print("detector axial shift",results.dshft)
# Get detected postions
detpos=pd.read_table("detected.dat",sep="\s+")
# create image
xmin= -halfd
xmax= halfd
ymin= -halfd
ymax= halfd
nx=100
ny=100
arr=images.binxy(detpos["YDET"],detpos["ZDET"],0,detpos["AREA"],
xmin,xmax,ymin,ymax,nx,ny)
xsam=(xmax-xmin)/nx
#
plt.figure()
#cm=plt.get_cmap("gnuplot2")
cm=plt.get_cmap("hot")
plt.imshow(arr,extent=[xmin,xmax,ymin,ymax],origin="lower",cmap=cm)
images.setfield(nx,xmin,xmax,ny,ymin,ymax)
rbeam=pdd*5
ibeam=rbeam/xsam
images.setpos(2,np.array([offrad*fl,0]))
b=images.beam(arr,ibeam,0,0)
hew=(b.hew*xsam/fl)*(180/np.pi)*3600
print("HEW arc sec",hew)
fwhm=(b.fwhm*xsam/fl)*(180/np.pi)*3600
print("FWHM arc sec",fwhm)
area=b.flux/100
print("area",area,"cm2")
psfarea=(b.hew*xsam/10)**2
fgain=area/psfarea
print("focusing gain",fgain)
#
plt.show()
