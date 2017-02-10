#!/usr/bin/env python
# Set up for SVOM MXT plates version 4
from __future__ import print_function
import numpy as np
# Focal length
flen= 1150
# Plate dimension (aperture)
wplate= 38
# Gap size
gap= 4
# Pore size
ps= 0.04
openf= 0.67
wall= ps/np.sqrt(openf)-ps
fibre= (ps+wall)*25
# Scale for L/d
xsi= 1.25
# number of plates on a side
nplate= 5
#
pitch= wplate+gap
#
xmin= -(nplate-1)*pitch/2
ymin= xmin
#
rmax= np.sqrt(pitch**2+(nplate*pitch/2)**2)
rmin= -1
hmcp= 110
print("rmin rmax",rmin,rmax)
#
nd=0
im=[]
ix=[]
iy=[]
xm=[]
ym=[]
wm=[]
hm=[]
lm=[]
tm=[]
am=[]
rm=[]
pm=[]
for i in range(nplate):
    xp= xmin+i*pitch
    for j in range(nplate):
        yp= ymin+j*pitch
        rp= np.sqrt(xp**2+yp**2)	
        if rp<rmax and rp>rmin:
            nd=nd+1
            im.append(nd)
            ix.append(i)
            iy.append(j)
            xm.append(xp)
            ym.append(yp)
            wm.append(wplate)
            hm.append(wplate)
            if rp<1:
                rr=lm[nd-2]
                lm.append(rr)
            else:
                lm.append(2.0*xsi*ps*flen/rp)
            tm.append(0)
            rm.append(rp)
            pm.append(np.arctan2(yp,xp))
print(nd,"plates")
im=np.array(im)
ix=np.array(ix)
iy=np.array(iy)
xm=np.array(xm)
ym=np.array(ym)
wm=np.array(wm)
hm=np.array(hm)
lm=np.array(lm)
tm=np.array(tm)
rm=np.array(rm)
pm=np.array(pm)
# set up deformations etc. for the plates
xd= np.array(im)
yd= np.array([1])
# ds is scaling factor for intrinsic slumping errors
ds=1
zx=np.empty(nd)
zx.fill(ds)
# de is the maximum thermoelelastic pore axial pointing error radians
# If -ve then is the rms pointing error radians on either axis chosen
# randomly
de=  -1.7*np.pi/180.0/60.0
zy=np.empty(nd)
zy.fill(de)
# trms is the rms pore rotation error about pore axis in radians
trms= 0.5*np.pi/180
za=np.empty(nd)
za.fill(trms)
# dshr is the shear distance for each pore mm used in multifibre model
# of the shear error
dshr= ps*0.07
zd=np.empty(nd)
zd.fill(dshr)
# qrms error of figure in arc mins to radians
qrms= 0.9*np.pi/180.0/60.0
fd=np.empty(nd)
fd.fill(qrms)
# reflecting surface definition
rough= 11
rind=1.4
brk= 10.0
gam= rind-1.0
nrough= rough**2/((1+1/gam)*brk**(-gam))
# Optical constants for bare MCP glass
#mspec<- "Si40 O98 K8 Na7 Pb6 Bi2"
#rho<- 3.3
# Optical constants for Iridium
mspec= "Ir"
rho= 22.65
# Set reflectivity type to be optical constants
ist=1
# Dummy arrays for angle and reflectivity
refs=np.array([1,1])
angs=np.array([0,90])
# ascii tabulation of plate data
atab=np.array([im,ix,iy,xm,ym,wm,hm,lm,rm,pm,tm])
np.savetxt("svom_mxt_plates_v4.dat",atab.transpose(),
header="im,ix,iy,xm,ym,wm,hm,lm,rm,pm,tm",
fmt=["%d","%d","%d","%f","%f","%f","%f","%f","%f","%f","%f"])
