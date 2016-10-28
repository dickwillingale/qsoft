#!/usr/bin/env python
# Set up Athena SPO modules with large adaptor
from __future__ import print_function
import numpy as np
from ir_b4cover import *
# Aperture radial limits
rmin=258.0
rmax= 1450.0
# rib spacing, radial and azimuthal pore pitch
wall=0.17
pr=0.605
rpitch=pr+wall
apitch=2.3
# iax controls the axial curvature of 1st and 2nd surfaces and principle plane
#
# Conical approximation
#iax=0
# Wolter I
#iax= 1
# Wolter I with spherical principal surface
#iax= 4
# Conical approximation with spherical principal surface
#iax= 5
# Wolter I with spherical principal surface and grid
#iax= 6
# Spherical principal plane with curved-plane curvature
#iax= 7
# Spherical principal plane with plane-curved curvature
iax= 8
# Spherical principal plane with -1/2, 5/2 curvature
#iax= 9
# Spherical principal plane with 1/2, 3/2 curvature
#iax= 10
# Spherical principal plane with -1, 3 curvature
#iax= 11
# Spherical principal plane with -1.5, 3.5 curvature
#iax= 12
# Spherical principal plane with -2, 4 curvature
#iax= 13
# reflecting surface roughness definition
rough= 5.0
rind= 1.4
brk= 1.
gam= rind-1.0
nrough= rough**2/((1+1/gam)*brk**(-gam))
# Reflectivity
ekev= np.array([1.0])
# Tabulated reflectivity for Ir + B4C overcoat
angs,refs=ir_b4c(ekev)
nag=angs.size
alpha=0
gamma=0
ist=2
# Optical constants for Iridium
#mspec= "Ir"
#rho= 22.65
#xop.xscat(mspec,rho,ekev,1)
#alpha=xop.alpha
#gamma=xop.gamma
#ist=1
#refs=np.array([1,1])
#angs=np.array([0,90])
#nag=2
#
# radial gap and module aperture height
gapr= 8.0
mr=54.1
# step size in radius
dring=mr+gapr
# azimuthal gap
gapw= 15.0
# width of spider arms
gaps= 20.0
# module aperture width outer
dang= 100.0
# module aperture width inner
dangi= 60.0
# module frame width (around module aperture)
wfr= 1.0
# focal length
flen=12000.0
# ratio of grazing angles
gra= 1.0
# Set up modules in rings
nring=np.int((rmax-rmin)/dring)
rring= np.empty(nring)
hring= np.empty(nring)
rsplit= np.empty(nring)
graz= np.empty(nring)
graz1= np.empty(nring)
graz2= np.empty(nring)
imod=0
nmod=1
im=[]
ir=[]
ii=[]
xm=[]
ym=[]
wm=[]
hm=[]
lm=[]
tm=[]
cm=[]
gm=[]
for i in range(nring):
    ring=rmin+(i+0.5)*dring
    rring[i]= ring
# Work in 60 degree sectors
    if(ring<500):
        na=np.int((ring*np.pi/3-gaps)/dangi)+1
    else:
        na=np.int((ring*np.pi/3-gaps)/dang)+1
    dth=(np.pi/3-gaps/ring)/na
    mw= dth*ring-gapw
    graz[i]=np.arctan(ring/flen)/2.0/(1.0+gra)
    ml=pr/graz[i]
    graz1[i]=np.arctan((ring-mr/2)/flen)/2.0/(1.0+gra)
    graz2[i]=np.arctan((ring+mr/2)/flen)/2.0/(1.0+gra)
    hring[i]= ml
    for k in range(na):
        for l in range(6):
            theta=k*dth+l*np.pi/3+gaps/ring/2+mw/ring/2
            imod=imod+1
            im.append(imod)
            ir.append(i)
            ii.append(k)
            xm.append(ring*np.cos(theta))
            ym.append(ring*np.sin(theta))
            wm.append(mr)
            hm.append(mw)
            lm.append(ml)
            tm.append(theta)
            cm.append(iax)
            gm.append(gra)
    rsplit[i]= rring[i]+dring/2
#    print(i,rsplit[i],graz[i]*180/np.pi*60,graz1[i]*180/np.pi*60,
#    graz2[i]*180/np.pi*60)
im=np.array(im)
ir=np.array(ir)
ii=np.array(ii)
# Radius of module aperture centres
xm=np.array(xm)
ym=np.array(ym)
rm=np.sqrt(xm**2+ym**2)
# azimuthal position of module centres
pm=np.arctan2(ym,xm)
# Add frame width
wm=np.array(wm)
hm=np.array(hm)
wm=wm+2*wfr
hm=hm+2*wfr
#
lm=np.array(lm)
tm=np.array(tm)
cm=np.array(cm)
gm=np.array(gm)
# grazing angle at module centres
gr=np.arctan(rm/flen)/4
#
nm=imod
print(nm,"modules")
print("minimum radius",np.min(rm)-mr/2)
print("maximum radius",np.max(rm)+mr/2)
# set up deformation vectors that give module alignment and pore figure errors
xd=np.array(im)+1
yd=np.array([1])
# module rotation about optical axis radians (3 arc secs)
rrms=(3.0/3600.0)*np.pi/180
# shift of module in aperture plane rms mm (applied to x and y in aperture)
xyrms=0.03
# axial module shift rms mm (includes error in kink angle)
zrms=1.0
# in-plane figure errors rms radians (1.2 arc secs)
brms=(1.2/3600)*np.pi/180
# out-of-plane figure errors rms radians (1.5 arc secs)
crms=(1.5/3600)*np.pi/180
# set up deformation vectors for modules
te=np.random.normal(0.0,rrms,nm)
dx=np.random.normal(0.0,xyrms,nm)
dy=np.random.normal(0.0,xyrms,nm)
dz=np.random.normal(0.0,zrms,nm)
db=np.empty(nm)
db.fill(brms)
dc=np.empty(nm)
dc.fill(crms)
# ASCII tabulation
atab=np.array([im,ir,ii,xm*0.001,ym*0.001,wm,hm,lm,rm])
np.savetxt("athena_ladapt_modules.dat",atab.transpose(),
header="im ir ii xm ym wm hm lm rm",
fmt=["%d","%d","%d","%f","%f","%f","%f","%f","%f"])
# estimate volume
vol=sum(wm*hm*lm)/1000
print("volume cm3",vol)
print("min max width",min(wm),max(wm))
print("min max height",min(hm),max(hm))
print("min max length",min(lm),max(lm))
#
