# X-ray Scattering Code for Python
from __future__ import print_function
import xscatfor
# Initialisation
xscatfor.qrs_init()
# Shareble object version of Wiscombe Mie scattering subroutine
# and Rayleigh-Gans approximation for scattering by dust
class mievs:
    def list(self):
        print("qext",self.qext)
        print("qsca",self.qsca)
        print("gqsc",self.gqsc)
        print("pmom",self.pmom)
        print("sforw",self.sforw)
        print("sback",self.sback)
        print("s1",self.s1)
        print("s2",self.s2)
        print("tforw",self.tforw)
        print("tback",self.tback)
        print("spike",self.spike)
def miev0(xx,crefin,perfct,mimcut,anyang,xmu,nmom,ipolzn,momdim,prnt):
    a=xscatfor.miev0(xx,crefin,perfct,mimcut,anyang,xmu,nmom,ipolzn,momdim,prnt)
    b=mievs()
    b.qext=a[0]
    b.qsca=a[1]
    b.gqsc=a[2]
    b.pmom=a[3]
    b.sforw=a[4]
    b.sback=a[5]
    b.s1=a[6]
    b.s2=a[7]
    b.tforw=a[8]
    b.tback=a[9]
    b.spike=a[10]
    return b
#print(xrayfor.miev0.__doc__)
# X-ray scattering from interstellar dust
def dustrings(data,derr,dsou,ekv,srate,dts,td,zd,amin,amax,qa,sig1):
    return xscatfor.qra_dustrings(data,derr,dsou,ekv,srate,dts,td,zd,
    amin,amax,qa,sig1)
class duststat:
    def __init__(self,data,derr,dsou,ekv,srate,dts,td,zd,sig1):
        self.data=data
        self.derr=derr
        self.dsou=dsou
        self.ekv=ekv
        self.srate=srate
        self.dts=dts
        self.td=td
        self.zd=zd
        self.sig1=sig1
        self.ncall=0
    def __call__(self,pars):
        angs,ndust,edust,model,chisq,ndof=dustrings(self.data,self.derr,
        self.dsou, self.ekv,self.srate,self.dts,self.td,self.zd,
        pars[0]*pars[1],pars[1],pars[2],self.sig1)
        self.ncall=self.ncall+1
        print(self.ncall,pars,chisq)
        return chisq
def dustthetascat(td,ds,zd):
    return xscatfor.qra_dustthetascat(td,ds,zd)
# X-ray optical properties
class xoptdat: pass
def xopt(mspec,rho,ekev,itype):
    a=xscatfor.qrt_xopt(len(mspec),mspec,rho,ekev,itype)
    b=xoptdat() 
    b.alpha=a[0]
    b.gamma=a[1]
    b.absl=a[2]
    b.f1=a[3]
    b.f2=a[4]
    return b
def xfresnel(alpha,gamma,angs):
    a=xscatfor.qrt_xfresnel(alpha,gamma,angs)
    b=xoptdat()
    b.rs=a[0]
    b.rp=a[1]
    b.runp=a[2]
    return b
def mlayer(angs,ekev,nr,ni,d,nper):
    a=xscatfor.qrt_mlayer(angs,ekev,nr,ni,d,nper)
    b=xoptdat()
    b.rsig=a[0]
    b.rpi=a[1]
    b.tsig=a[2]
    b.tpi=a[3]
    return b
