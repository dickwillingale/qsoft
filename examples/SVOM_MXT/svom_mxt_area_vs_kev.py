#!/usr/bin/env python
# Ray tracing of SVOM MXT optic V4 - area vs. energy
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
# Set up plates
from svom_mxt_plates_v4 import *
# Global loss factor
floss= 0.78
# reflecting surface index
isu=1
if ist==1:
    print("reflectivity using optical constant")
else:
    print("reflectivity using look-up table")
    na=angs.size
    for i in range(angs.size):
        print(angs[i],refs[i])
print("surface roughness parameters",nrough,rind,brk)
#
ppp= np.array([rm,pm,-pm,wm,hm,lm])
pp=ppp.transpose()
nplts= xm.size
pitch= ps+wall
plm= 0
plx= 10.0
ipack=7
idef= 1
lekev=np.arange(np.log10(0.19),np.log10(11),0.01)
ekev=np.power(10.0,lekev)
xop= xscat.xopt(mspec,rho,ekev,1)
alpha= xop.alpha
gamma= xop.gamma
ne=ekev.size
iss=1
# set up source
ccy=np.sin(0.0*np.pi/180)
ccz=np.sin(0*np.pi/180)
di=np.array([-np.sqrt(1.0-ccy**2-ccz**2),ccy,ccz])
hlen=1.05*hmcp
slim=np.array([-hlen,-hlen,hlen,hlen,0,0])
nray=500000
psrc=np.array([flen*1.05,0,0])
# set position of square pore MCP
pmcp=np.array([flen,0,0])
# Plane detector with cartesian limits
pix= 0.075
ni= 512
hdet= ni*pix/2
dlim=np.array([-hdet,-hdet,hdet,hdet])
# rotate detector axes wrt the plate reference axes
thetadet=-0*np.pi/180
rxdet=np.array([0,np.cos(thetadet),np.sin(thetadet)])
farea=np.empty(ne)
tarea=np.empty(ne)
fwhm=np.empty(ne)
for i in range(ne):
    xsrt.reset()
#   xsrt.rseed(55)
    xsrt.surface(iss,ist,ekev[i],nrough,brk,rind,alpha[i],gamma[i],
    angs,refs,0,0,0)
    xsrt.deform(idef,1,5,nd,1)
    xsrt.defmat(idef,1,xd,yd,zx)
    xsrt.defmat(idef,2,xd,yd,zy)
    xsrt.defmat(idef,3,xd,yd,za)
    xsrt.defmat(idef,4,xd,yd,zd)
    xsrt.defmat(idef,5,xd,yd,fd)
    xsrt.source(2,di,nn,psrc,sn,rx,slim,0,nray,0)
    xsrt.sqpore(pmcp,sn,rx,2*flen,ipack,hmcp,pitch,wall,plx,idef,
    iss,plm,plx,fibre,pp)
    xsrt.detector(2,nn,sn,rxdet,dlim,0)
    results=xsrt.trace(0,hdet,-1)
# get detected positions
    detpos=pd.read_table("detected.dat",sep="\s+")
# create image
    a=images.binxy(detpos["YDET"],detpos["ZDET"],0,detpos["AREA"],
    -hdet,hdet,-hdet,hdet,ni,ni)
#
    images.setfield(ni,-hdet,hdet,ni,-hdet,hdet)
    xp= 0/60*np.pi/180*flen
    images.setpos(2,np.array([-xp,0]))
# Beam 6 arc minute radius (about FWHM)
    rbm=(6/60)*(np.pi/180)*flen
    rpix=rbm/pix
    b=images.beam(a.data,rpix,0,-1)
#   fwhm[i]=b.fwhm*pix/flen*180*60/np.pi
    fwhm[i]= b.fit.x[3]*2.355*pix/flen*180*60/np.pi
    farea[i]=b.flux*0.01*floss
# Get total area
    b=images.beam(a.data,ni*2,0,0)
    tarea[i]=b.flux*0.01*floss
    print(ekev[i],farea[i],tarea[i],fwhm[i])
# Write tabulation file
atab=np.array([ekev,farea,tarea,fwhm])
np.savetxt("svom_mxt_area_vs_kev_v4.dat",atab.transpose(),
header="ekev,farea,tarea,fwhm",fmt=["%f","%f","%f","%f"])
