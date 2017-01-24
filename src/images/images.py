# Image processing routines
from __future__ import print_function
import imagesfor 
import numpy as np
import matplotlib.pylab as plt
# Initialisation and reseting
def init():
    imagesfor.qri_init()
def reset():
    imagesfor.qri_init()
init()
def rectangles(x,y,w,h,t):
    cth=np.cos(t)
    sth=np.sin(t)
    print("length x",len(x))
    print("length t",len(t))
    xa=w*cth
    xb=w*sth
    ya=h*cth
    yb=h*sth
    x1=x+xa+yb
    x2=x+xa-yb
    x3=x-xa-yb
    x4=x-xa+yb
    y1=y+xb-ya
    y2=y+xb+ya
    y3=y-xb+ya
    y4=y-xb-ya
    for i in range(len(x1)):
        plt.plot([x1[i],x2[i],x3[i],x4[i],x1[i]],
        [y1[i],y2[i],y3[i],y4[i],y1[i]],
        color='k',linestyle='-',linewidth=1)
def binxy(x,y,iq,w,xleft,xright,ybot,ytop,nx,ny):
    a=imagesfor.qri_binxy(x,y,iq,w,xleft,xright,ybot,ytop,nx,ny)
    return a
def quaderr(x0,y0,x1,y1,y):
    # Quadratic estimator for confidence limit
    # x0,y0         parameter and statistic at minimum
    # x1,y1         parameter and statistic near minimum
    # y             required statistic value
    # returns       estimate of parameter corresponding to y
    a=(y1-y0)/(x1-x0)**2
    b=-2.0*a*x0
    c=y0-a*x0**2-b*x0
    ba=b**2-4.0*a*(c-y)
    if ba>=0 and a!=0:
        return np.sqrt(ba)/(2.0*a)
    else:
        return 0.0
def srchmin(pars,pl,ph,stat,delstat,derr):
    # Search for minimum statistic and return best fit parameters and
    # confidence  limits of the parameters
    # pars        initial parameter values
    # pl          hard lower limit of parameter values
    # ph          hard upper limit of parameter values
    # stat        the statistic function to be minimised
    # delstat     the change in statistic for confidence limits
    # derr        initial error estimates for parameters
    #             0 fixed, >0 to estimate confidence range
    # returns     list from optim() plus confidence limits parlo and parhi
    npa=len(pars)
    # make local copies of hard limits
    prl=np.copy(pl)
    prh=np.copy(ph)
    # find statistic minimum and estimate hessian matrix
    bnds=np.swapaxes([prl,prh],0,1)
    ft=optimize.minimize(stat,pars,method='L-BFGS-B',bounds=bnds)
    # estimate errors on parameters
    ft.delstat=delstat
    # If hessian matrix available use to estimate steps
    #dig= np.absolute(np.diagonal(ft.hess))
    # If not estimate 2nd derivitive from initial error estimate
    dig=np.zeros(npa)
    for k in range(npa):
        if derr[k]>0:
            dig[k]= 2.0*delstat/derr[k]**2
    delpar= np.empty(npa)
    # fix hard limits of parameters for which we don't want error estimate
    for k in range(npa):
        # Trap zero on hessian diagonal
        if dig[k]==0:
            delpar[k]=0
            if derr[k]==0:
                prl[k]=ft.x[k]-ft.x[k]/200
                prh[k]=ft.x[k]+ft.x[k]/200
        else:
            delpar[k]=np.sqrt(2.0*delstat/dig[k])
            if derr[k]==0:
                prl[k]=ft.x[k]-delpar[k]/200
                prh[k]=ft.x[k]+delpar[k]/200
    print("Min statistic",ft.fun)
    print("best fit parameters",ft.x)
    print("diag",dig)
    print("delpar",delpar)
    nfr= np.sum(derr>0)
    print("Find limits for ",nfr," parameters")
    parhi=np.empty(npa)
    parlo=np.empty(npa)
    for k in range(npa):
        if derr[k]>0 and delpar[k]>0:
            print("searching for error range",k,ft.x[k],delpar[k],prl[k],prh[k])
            # upper limit
            epar= 0
            while epar==0:
                pval=min(ft.x[k]+delpar[k],prh[k])
                part=np.copy(ft.x)
                partl=np.copy(prl)
                parth=np.copy(prh)
                part[k]=pval
                partl[k]=pval-delpar[k]/200
                parth[k]=pval+delpar[k]/200
                bnds=np.swapaxes([partl,parth],0,1)
                ftp=optimize.minimize(stat,part,method='L-BFGS-B',bounds=bnds)
                if ftp.fun<ft.fun:
                    ft.x=ftp.x
                    ft.fun=ftp.fun
                    print("new min found",ft.fun)
                    print("new fit parameters",ftp.x)
                    if ftp.x[k]>prh[k]:
                        epar=999
                else:
                    if ftp.x[k]>prh[k]:
                        epar= 999
                    else:
                        epar=quaderr(ft.x[k],ft.fun,pval,ftp.fun,ft.fun+delstat)
            print("upper",ft.x[k],ft.fun,pval,ftp.fun,ft.fun+delstat,epar)
            parhi[k]=min(ft.x[k]+epar,prh[k])
            # lower limit
            epar= 0
            while epar==0:
                pval=max(ft.x[k]-delpar[k],prl[k])
                part=np.copy(ft.x)
                partl=np.copy(prl)
                parth=np.copy(prh)
                part[k]=pval
                partl[k]=pval-delpar[k]/200
                parth[k]=pval+delpar[k]/200
                bnds=np.swapaxes([partl,parth],0,1)
                ftp=optimize.minimize(stat,part,method='L-BFGS-B',bounds=bnds)
                if ftp.fun<ft.fun:
                    ft.x=ftp.x
                    ft.fun=ftp.fun
                    print("new min found",ft.fun)
                    print("new fit parameters",ftp.x)
                    if ftp.x[k]>prh[k]:
                        epar=999
                else:
                    if ftp.x[k]>prh[k]:
                        epar= 999
                    else:
                        epar=quaderr(ft.x[k],ft.fun,pval,ftp.fun,ft.fun+delstat)
            print("lower",ft.x[k],ft.fun,pval,ftp.fun,ft.fun+delstat,epar)
            parlo[k]=max(ft.x[k]-epar,prl[k])
        else:
            parhi[k]=0.0
            parlo[k]=0.0
    ft.parlo=parlo
    ft.parhi=parhi
    return ft
class bdata: pass
def beam(arr,rbeam,blev,bvar):
    a=imagesfor.qri_beam(arr,rbeam,blev,bvar)
    b=bdata()
    b.nsam=a[0]
    b.bflux=a[1]
    b.bsigma=a[2]
    b.flux=a[3]
    b.fsigma=a[4]
    b.peak=a[5]
    b.cen=a[6]
    b.tha=a[7]
    b.rmsa=a[8]
    b.rmsb=a[9]
    b.fwhm=a[10]
    b.hew=a[11]
    b.w90=a[12]
    b.fwhmp=a[13]
    b.hewp=a[14]
    b.w90p=a[15]
    b.fwhmc=a[16]
    b.hewc=a[17]
    b.w90c=a[18]
    return b
def setfield(nx,xleft,xright,ny,ybot,ytop):
    imagesfor.qri_setfield(nx,xleft,xright,ny,ybot,ytop)
def setsky(xtodeg,ytodeg,ipr,mjd,ra,dec,roll):
    imagesfor.qri_setsky(xtodeg,ytodeg,ipr,mjd,ra,dec,roll)
def setpos(ipos,p):
    imagesfor.qri_setpos(ipos,p)
class iposition: pass
def getpos():
    a=imagesfor.qri_getpos()
    b=iposition()
    b.pix=a[0]
    b.xyl=a[1]
    b.aes=a[2]
    b.equ=a[3]
    b.ecl=a[4]
    b.gal=a[5]
    return b
