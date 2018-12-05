## @package astro
# Functions for astronomical and astrophysical calculations and analysis
import numpy as np
from astrofor import *
# Initialisation and reseting
def init():
    """Initialise common blocks on import"""
    astrofor.qrx_init()
def reset():
    """Reset common blocks to initial condition"""
    astrofor.qrx_init()
init()
# Bin edges and centres
def ctob(x):
    """Convert n bin centres to n+1 bin boundaries
        x      bin centres length n
    return bin boundaries length n+1
    """
    n=len(x)
    xx=(x[1:n]+x[0:n-1])/2.0
    x0=np.array([2*x[0]-xx[0]])
    xn=np.array([2*x[n-1]-xx[n-2]])
    xout=np.concatenate((x0,xx))
    return np.concatenate((xout,xn))
def btoc(x):
    """Convert n+1 bin boundaries to n bin centres
        x      bin boundaries length n+1
    return bin centres length n
    """
    n=len(x)
    return (x[1:n]+x[0:n-1])/2.0
# Continua spectra
def brems(ekev,t):
    """Bremsstrahlung spectrum
        ekev    array of photon energies keV
        t       temperature keV
    return  Bremsstrahlung continuum photons/keV at energies ekev
    """
    return astrofor.qrx_brems(ekev,t)
# Interstellar absorption
def habs(cd,ekev):
    """X-ray absorption by a Hydrogen column density
        cd      hydrogen column 10**21 cm-2
        ekev    array of photon energies keV
    return  array of absorption factors
    """
    return astrofor.qra_habs(cd,ekev)
def ismtau(nh,z,ekev):
    """X-ray optical depth of cold ISM gas
        nh      column density 10**21 cm-3 at redshift z
        z       redshift of source
        ekev    array of photon energies keV
    return  optical depth for each energy
    """
    a=astrofor.qrx_ismtau(nh,z,ctob(ekev))
    return a[0]
def iismtau(nh,z,tk,pl,ist,ekev):
    """X-ray optical depth of ionized ISM gas
        nh      column density 10**21 cm-3 at redshift z
        z       redshift of source
        tk      temperature Kelvin
        pl      powerlaw index of continuum spectrum
        ist     ionization state = L/nR^2
        ekev    array of photon energies keV
    return  optical depth for each energy
    """
    a=astrofor.qrx_iismtau(nh,z,tk,pl,ist,ctob(ekev))
    return a[0]
def igmtau(n0,z,h0,omegam,omegal,ekev):
    """X-ray optical depth of IGM gas
        n0      number density cm-3 at z=0
        z       redshift of source
        h0      cosmological H0 km s-1 Mpc-1
        omegam  cosmological omegam
        omegal  cosmological omegal
        ekev    array of photon energies keV
    return  optical depth for each energy
    """
    a=astrofor.qrx_igmtau(n0,z,h0,omegam,omegal,ctob(ekev))
    return a[1]
def iigmtau(n0,pl,tk,ist,z,h0,omegam,omegal,ekev):
    """ X-ray optical depth of ionized IGM gas
        n0     number density cm-3 at z=0
        pl     powerlaw index of continuum spectrum
        tk     temperature Kelvin
        ist    ionization state L/nR^2
        z      redshift of source
        h0     cosmological H0 km s-1 Mpc-1
        omegam cosmological omegam
        omegal cosmological omegal
        ekev   array of photon energies keV
    return optical depth for each energy
    """
    a=astrofor.qrx_iigmtau(n0,pl,tk,ist,z,h0,omegam,omegal,ctob(ekev))
    return a[1]
def lyftau(n0,dind,m0,mind,z,h0,omegam,omegal,ekev):
    """X-ray optical depth of Lyman Forest
        n0     number density cm-3 at z=0
        dind   density index wrt z+1
        m0     metalicity log[X/H] at z=0
        mind   metalicity log[X/H] index wrt z
        z      redshift of source
        h0     cosmological H0 km s-1 Mpc-1
        omegam cosmological omegam
        omegal cosmological omegal
        ekev   array of photon energies keV
    return optical depth for each energy
    """
    a=astrofor.qrx_lyftau(n0,dind,m0,mind,z,h0,omegam,omegal,ctob(ekev))
    return a[1]
def lyftauvz(n0,dind,m0,mind,h0,omegam,omegal,ekev,z):
    """X-ray optical depth of Lyman Forest vs. z
        n0     number density cm-3 at z=0
        dind   density index wrt z+1
        m0     metalicity log[X/H] at z=0
        mind   metalicity log[X/H] index wrt z
        h0     cosmological H0 km s-1 Mpc-1
        omegam cosmological omegam
        omegal cosmological omegal
        ekev   photon energy keV
        z      array of redshifts
    return optical depth at ekev for each redshift
    """
    return astrofor.qrx_lyftauvz(n0,dind,m0,mind,h0,omegam,omegal,ekev,z)
def iigmtauvz(n0,dind,m0,mind,pl,tk,ist,h0,omegam,omegal,ekev,z):
    """X-ray optical depth of ionized IGM gas vs. redshift
        n0     number density cm-3 at z=0
        dind   density index wrt z+1
        m0     metalicity log[X/H] at z=0
        mind   metalicity log[X/H] index wrt z
        pl     powerlaw index of continuum spectrum
        tk     temperature Kelvin
        ist    ionization state = L/nR^2
        h0     cosmological H0 km s-1 Mpc-1
        omegam cosmological omegam
        omegal cosmological omegal
        ekev   energy keV
        z      array of redshift values
    return optical depth at ekev for each redshift
    """
    return astrofor.qrx_iigmtauvz(n0,dind,m0,mind,pl,tk,ist,h0,omegam,omegal,ekev,z)
def setabund(abun,amet):
    """Set abundances for XSPEC routines
        abun   XSPEC abundance table
        amet   metalicity
    """
    astrofor.qrx_setabund(abun,amet)
# Cosmology
class cosmic: pass
def cosmo(h0,omegam,omegal,zmax):
    """Calculation of cosmological quantities using equations from David Hogg astro-ph/9905116
        h0     cosmological H0 km s-1 Mpc-1
        omegam cosmological omegam
        omegal cosmological omegal
        zmax   maximum redshift
    return list with following at samples dz=0.01
        z      redshift values
        ez     scaling function
        dc     comoving line-of-sight distance
        dm     transverse comoving distance
        da     angular diameter distance
        dl     luminosity distance
        distm  distance modulus
        dvc    comoving volume element
        tlbak  look back time
        vc     integrated comoving volume over whole sky
        thsec  Hubble time in seconds
        thgyr  Hubble time in Giga years
        dhmpc  Hubble length in Mpc
    """
    dz=0.01
    maxnum=int(zmax/dz)
    a=astrofor.qra_cosmo(h0,omegam,omegal,dz,maxnum)
    b=cosmic()
    b.h0=h0
    b.omegam=omegam
    b.omegal=omegal
    b.z=a[0]
    b.ez=a[1]
    b.dc=a[2]
    b.dm=a[3]
    b.da=a[4]
    b.dl=a[5]
    b.distm=a[6]
    b.dvc=a[7]
    b.tlbak=a[8]
    b.vc=a[9]
    b.thsec=a[10]
    b.thgyr=a[11]
    b.dhmpc=a[12]
    return b
class kcorr: pass
def kcorrb(e1src,e2src,e1obs,e2obs,z,gamma1,gamma2,ecobs):
    """K-correction using numerical integration of Band function
        e1src  lower source frame energy
        e2src  upper source frame energy
        e1obs  array of lower observed energies
        e2obs  array of upper observed energies
        z      array of redshifts
        gamma1 array of observed photon indices below Ec
        gamma2 array of observed photon indices above Ec
        ecobs  array of observed Ec energies
    return list with the following for each object
        kcorr  K-correction factor 
        obint  integral over observed band
        boint  integral over Eiso band in observed frame
    """
    a=astrofor.qra_kcorrb(e1src,e2src,e1obs,e2obs,z,gamma1,gamma2,ecobs)
    b=kcorr()
    b.e1src=e1src
    b.e2src=e2src
    b.e1obs=e1obs
    b.e2obs=e2obs
    b.z=z
    b.gamma1=gamma1
    b.gamma2=gamma2
    b.ecobs=ecobs
    b.kcorr=a[0]
    b.obint=a[1]
    b.boint=a[2]
    return b
