#!/usr/bin/env python
# Functions for dust halo X-ray scattering
import numpy as np
#
def thetascat(td,ds,zd):
    # Calculate scattering angle radians
    # td delay time seconds after source flare (assumed delta function)
    # ds distance to source parsecs (convert to m - 3.086e16/parsec)
    # zd fractional distance to dust 
    return np.sqrt(td*2.0*3.0e8/(1.0-zd)/(zd*ds*3.086e16))
def scatfun(x):
    return ((np.sin(x)/x**2-np.cos(x)/x)/x)**2
def diffsec(ekv,amin,amax,qa,sig1,ts):
    # differential scattering cross section of dust grains
    # Based on Mauche and Gorenstein 1986 and Smith and Dwek 1998
    # ekv photon energy keV
    # amin, amax grain size range microns
    # qa grain size distribution index
    # sig1 diff cross section at 1 keV for grain size 0.1 microns
    # ts scattering angle radians
    # Loop to integrate over grain size distribution
    na=20
    dela=(amax-amin)/na
    # Normalisation for the integral of N(a)=A.a**qa
    q1=1.0-qa
    aint=dela*q1/(amax**q1-amin**q1)
    # Peak value of (j1(x)/x)**2 distribution is 1/9
    # Factor 1.0e6 because sig1 for a=0.1 microns
    anorm=9.0e6*sig1*aint
    # Note convert keV to wavelength in microns lam=12.4e-4/ekv
    aol=2.0*np.pi*ts*ekv/12.4e-4
    sig=0.0
    for k in range(na):
        agr=amin+dela*(k+0.5)
        sig=sig+scatfun(aol*agr)*agr**(6-qa)
    return sig*anorm
def gaussian(x,sig):
    return np.exp(-x**2/2.0/sig**2)
def rings(data,derr,dsou,ekv,srate,tdlo,tdhi,zd,amin,amax,qa,sig1):
    # Model flux cts/s in expanding rings
    # data        array of rings flux data cts/s [nrings,ndata]
    # derr        array of errors on data
    # dsou        distance to source PC
    # ekv         array of energies keV (equally spaced across observed band)
    # srate       source spectrum cts/s/keV
    # tdlo,tdhi   delay times of observation secs
    # zd          fraction of source distance to rings
    # amin,amax   grain size radius range microns
    # qa          grain size distribution index N(a)=A.a^-qa
    # sig1        differential cross-section of 1 grain cm2, 1 keV, 0.1 microns
    # return (angs,ndust,model,chisq,ndof)
    # angs        observed angles of rings radians
    # ndust       column density of dust for each ring grains cm-2
    # edust       error on column density of dust for each ring grains cm-2
    # model       array of predicted cts/s for rings
    # chisq       Chi-squared between data and model array
    # ndof        number of degrees of freedom
    # arrays for results
    nrings=len(zd)
    ndata=len(tdlo)
    ndust=np.empty(nrings)
    edust=np.empty(nrings)
    model=np.empty([nrings,ndata])
    angs=np.empty([nrings,ndata])
    # 
    aobs=np.empty(ndata)
    drate=np.empty(ndata)
    esam=ekv[1]-ekv[0]
    ne=len(ekv)
    # Loop for rings
    for i in range(nrings):
        # Loop for data at delay times
        for k in range(ndata):
            # fitted observed angle range for delay time
            asclo=thetascat(tdlo[k],dsou,zd[i])
            aschi=thetascat(tdhi[k],dsou,zd[i])
            # observed ring angle
            aobs[k]=(asclo+aschi)*0.5*(1-zd[i])
            drate[k]=0.0
            # Loop over energy
            for j in range(ne):
                # differential cross-section for delay time
                xlo=diffsec(ekv[j],amin,amax,qa,sig1,asclo)
                xhi=diffsec(ekv[j],amin,amax,qa,sig1,aschi)
                # integrate over solid angle of ring
                solidangle=(aschi**2-asclo**2)*np.pi
                sigma=(xlo+xhi)*0.5*solidangle
                drate[k]=drate[k]+sigma*srate[j]
            drate[k]=drate[k]*esam/(tdhi[k]-tdlo[k])
        sdat=sum(data[i])
        ndust[i]=sdat/sum(drate)
        edust[i]=ndust[i]*np.sqrt(sum(derr[i]**2))/sdat
        model[i]=drate*ndust[i]
	angs[i]=aobs
    # Calculate Chi-squared
    chisq=sum((data.flatten()-model.flatten())**2/derr.flatten()**2)
    ndof=ndata*nrings-(nrings+3)
    return (angs,ndust,edust,model,chisq,ndof)
class stat:
    def __init__(self,data,derr,dsou,ekv,srate,tdlo,tdhi,zd,sig1):
        self.data=data
        self.derr=derr
        self.dsou=dsou
        self.ekv=ekv
        self.srate=srate
        self.tdlo=tdlo
        self.tdhi=tdhi
        self.zd=zd
        self.sig1=sig1
    def __call__(self,pars):
        print("parameters",pars)
        angs,ndust,edust,model,chisq,ndof=rings(self.data,self.derr,self.dsou,
        self.ekv,self.srate,self.tdlo,self.tdhi,self.zd,
        pars[0],pars[1],pars[2],self.sig1)
	print("chisq",chisq)
	print("ndust",ndust)
        return chisq
