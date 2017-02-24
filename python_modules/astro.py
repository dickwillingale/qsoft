# Astronomy routies for qpy
import numpy as np
import astrofor
# Initialisation and reseting
def init():
    astrofor.qrx_init()
def reset():
    astrofor.qrx_init()
init()
# Bin edges and centres
def ctob(x):
    n=len(x)
    xx=(x[1:n]+x[0:n-1])/2.0
    x0=np.array([2*x[0]-xx[0]])
    xn=np.array([2*x[n-1]-xx[n-2]])
    xout=np.concatenate((x0,xx))
    return np.concatenate((xout,xn))
def btoc(x):
    n=len(x)
    return (x[1:n]+x[0:n-1])/2.0
# Interstellar absorption
def habs(cd,ekev):
    return astrofor.qra_habs(cd,ekev)
def ismtau(nh,z,ekev):
    a=astrofor.qrx_ismtau(nh,z,ctob(ekev))
    return a[0]
def iismtau(nh,z,tk,pl,ist,ekev):
    a=astrofor.qrx_iismtau(nh,z,tk,pl,ist,ctob(ekev))
    return a[0]
def igmtau(n0,z,h0,omegam,omegal,ekev):
    a=astrofor.qrx_igmtau(n0,z,h0,omegam,omegal,ctob(ekev))
    return a[1]
def iigmtau(n0,pl,tk,ist,z,h0,omegam,omegal,ekev):
    a=astrofor.qrx_iigmtau(n0,pl,tk,ist,z,h0,omegam,omegal,ctob(ekev))
    return a[1]
def lyftau(n0,dind,m0,mind,z,h0,omegam,omegal,ekev):
    a=astrofor.qrx_lyftau(n0,dind,m0,mind,z,h0,omegam,omegal,ctob(ekev))
    return a[1]
def lyftauvz(n0,dind,m0,mind,h0,omegam,omegal,ekev,z):
    return astrofor.qrx_lyftauvz(n0,dind,m0,mind,h0,omegam,omegal,ekev,z)
def iigmtauvz(n0,dind,m0,mind,pl,tk,ist,h0,omegam,omegal,ekev,z):
    return astrofor.qrx_iigmtauvz(n0,dind,m0,mind,pl,tk,ist,h0,omegam,omegal,ekev,z)
def setabund(abun,amet):
    astrofor.qrx_setabund(abun,amet)
# Cosmology
class cosmic: pass
def cosmo(h0,omegam,omegal,zmax):
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
