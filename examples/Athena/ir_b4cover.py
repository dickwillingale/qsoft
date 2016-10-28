#!/usr/bin/env python
# Function to calculate reflectivity of Ir with B4C overcoat
import numpy as np
import xscat
def ir_b4c(ekev):
# Ir refractive index
    irdensity= 22.65
    fir= xscat.xopt("Ir",irdensity,ekev,1)
    nri= -fir.alpha/2.0+1.0
    nii= fir.gamma/2.0
# B4C refractive index
    b4cdensity= 2.52
    fbc= xscat.xopt("B4 C",b4cdensity,ekev,1)
    nrc= -fbc.alpha/2.0+1.0
    nic= fbc.gamma/2.0
# Set up Ir base with B4C overcoat
    nrlay= np.empty(3)
    nilay= np.empty(3)
    dlay= np.empty(3)
    nrlay[0]= 1
    nilay[0]= 0.0
    dlay[0]= 0.0
    dlay[1]= 80.0
    dlay[2]= 0.0
#
    amin= 87.0
    amax= 89.9
    asam= 0.04
    iangle= np.arange(amin,amax,asam)
    gangle= -iangle+90.0
#
    na=iangle.size
    ne= ekev.size
    refs= np.empty([na,ne])
#
    for k in range(ne):
        nrlay[1]=nrc[k]
        nilay[1]=nic[k]
        nrlay[2]=nri[k]
        nilay[2]=nii[k]
        mdat= xscat.mlayer(iangle,ekev[k],nrlay,nilay,dlay,1)
        refs[:,k]= (mdat.rsig+mdat.rpi)*0.5
#
    return iangle,refs
