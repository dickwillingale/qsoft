## @package images
# Functions for image display and manipulation
from __future__ import print_function
import imagesfor 
import locator
import numpy as np
import matplotlib.pylab as plt
from scipy.stats import chi2
from scipy import optimize
# Initialisation and reseting
def init():
    """Initialise common blocks on import"""
    imagesfor.qri_init()
def reset():
    """Reset common blocks to initial condition"""
    imagesfor.qri_init()
init()
def rectangles(x,y,w,h,t):
    """Draw rectangles

    Args:
        x: centres in local x
        y: centres in local y
        w: widths 
        h: heights
        t: rotation angles radians
    """
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
    """x-y event binning to form an image

    Args:
        x:      array of x positions
        Y:      array of y positions
        iq:     array of quality values (0 for OK)
        w:      array of weights
        xleft:  minimum x value (left edge) for image array
        xright: maximum x value (rigth edge) for image array
        ybot:   minimum y value (bottom edge) for image array
        ytop:   maximum y value (top edge) for image array
        nx:     size of 1st dimension of image array (number of rows)
        ny:     size of 2nd dimension of image array (number of columns)

    Returns:
        | image array
        | R version returns an image object - image array object$data_array
    """
    a=imagesfor.qri_binxy(x,y,iq,w,xleft,xright,ybot,ytop,nx,ny)
    return a
def quaderr(x0,y0,x1,y1,y):
    """Quadratic estimator for confidence limit

    Args:
        x0: parameter at minimum
        y0: statistic at minimum
        x1: parameter near minimum
        y1: statistic near minimum
        y:  required statistic value

    Returns:
        estimate of parameter corresponding to y
    """
    a=(y1-y0)/(x1-x0)**2
    b=-2.0*a*x0
    c=y0-a*x0**2-b*x0
    ba=b**2-4.0*a*(c-y)
    if ba>=0 and a!=0:
        return np.sqrt(ba)/(2.0*a)
    else:
        return 0.0
def srchmin(pars,pl,ph,stat,delstat,derr):
    """Search for minimum statistic and return best fit parameters and confidence limits of the parameters

    Args:
        pars:        initial parameter values
        pl:          hard lower limit of parameter values
        ph:          hard upper limit of parameter values
        stat:        the statistic function to be minimised
        delstat:     the change in statistic for confidence limits
        derr:        initial error estimates for parameters, 0 fixed, >0 to estimate confidence range

    Returns:
        list from optim() plus confidence limits parlo and parhi

    The statistic function call must return the value of a statistic with call of form stat(pars) (e.g. images.peakchisq(pars))
    """
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
    ft.par=ft.x
    ft.parlo=parlo
    ft.parhi=parhi
    return ft
class bdata: pass
def peakchisq(fpars):
    """Chi-squared for image peak fitting

    Args:
        fpars:       array of fitting parameters

    | par 0: normalisation (value at peak)
    | par 1: x-centre (pixel position)
    | par 2: y-centre (pixel position)
    | par 3: Lorentian width (pixels) 24th Sept. 2017 RW
    
    Returns:
        Chi-squared value

    Used by function beam() which loads Fortran common block BEAMFIT with image data.
    """
    chis=imagesfor.qri_peakchisq(fpars)
    return chis
def lorentzian(x,p):
    """Lorentzian profile

    Args:
        x:          array of x values
        p:          array of fitting parameters

    | par 0: normalisation (value at peak)
    | par 1: x-centre (pixel position)
    | par 2: Lorentian width (pixels)
    
    Returns:
        array of function values evaluated at x
    """
    return p[0]/(1+np.square((x-p[1])*2.0/p[2]))
def king_profile(x,p):
    """King function (modified Lorentzian) profile

    Args:
        x:          array of x values
        p:          array of fitting parameters

    | par 0: normalisation (value at peak)
    | par 1: x-centre (pixel position)
    | par 2: full width half maximum (pixels)
    | par 3: power index (1 for Lorentzian)

    Returns:
        array of function values evaluated at x
    """
    G=p[2]/np.sqrt(np.power(2.0,1.0/p[3])-1.0)
    return p[0]/np.power(1+np.square((x-p[1])*2.0/G),p[3])
def gaussian(x,p):
    """Gaussian profile

    Args:
        x:          array of x values
        p:          array of fitting parameters
    
    | par 0: normalisation (value at peak)
    | par 1: x-centre (pixel position)
    | par 2: Gaussian full width half maximum (pixels)

    Returns:
        array of function values evaluated at x
    """
    return p[0]*np.exp(-np.square(x-p[1])/(2*np.square(p[2]/2.37)))
def sqbeam(arr,hbeam,blev,bvar):
    """Analysis of source above background within a square beam

    Args:
        arr:     image array
        hbeam:   half width of square beam in pixels
        blev:    average background level per pixel (to be subtracted)
        bvar:    variance on blev (-ve for counting statistics)

    Returns:
        | list with the following
        | **nx,ny**:    dimension of beam pixels (truncated if falls off edge)
        | **xpi,xpr**:  x output arrays, pixel position and flux
        | **ypi,ypr**:  y output arrays, pixel position and flux
        | **bflux**:    background in beam (e.g. counts)
        | **bsigma**:   standard deviation of background
        | **flux**:     source flux above background in beam (e.g. counts)
        | **fsigma**:   standard deviation of source flux
        | **peak**:     source x,y peak position
        | **cen**:      source x,y centroid position
        | **rmsx**:     rms width in x (pixels) about centroid
        | **rmsy**:     rms width in y (pixels) about centroid
        | **pi5**:      5% x,y position
        | **pi25**:     25% x,y position
        | **med**:      median (50%) x,y position
        | **pi75**:     75% x,y position
        | **pi95**:     95% x,y position
        | **hewx**:     HEW (half energy width) x (pixels)
        | **hewy**:     HEW (half energy width) y (pixels)
        | **w90x**:     W90 (90% width) x (pixels)
        | **w90y**:     W90 (90% width) y (pixels)
        | **fitx**:     parameters from x profile fit using king_profile() 
        | **fity**:     parameters from y profile fit using king_profile() 
    |
    | Fits performed if bvar!=0.
    | Parameters are saved in the lists fitx and fity
    |   0: peak value (no error range calculated)
    |   1: peak X or Y pixel position (no error range calculated)
    |   2: Lorentzian width including 90% upper and lower bounds
    |   3: power index (1 for Lorentzian)
    | The position of the sqbeam is the current position within the image.
    | Use function setpos() to set the current position.
    """
    npp=1000
    a=imagesfor.qri_sqbeam(arr,hbeam,blev,bvar,npp)
    b=bdata()
    b.xpi=a[0]
    b.xpr=a[1]
    b.ypi=a[2]
    b.ypr=a[3]
    b.buf=a[4]
    b.nx=a[5]
    b.ny=a[6]
    b.bflux=a[7]
    b.bsigma=a[8]
    b.flux=a[9]
    b.fsigma=a[10]
    b.peak=a[11]
    b.cen=a[12]
    b.rmsx=a[13]
    b.rmsy=a[14]
    b.pi5=a[15]
    b.pi25=a[16]
    b.med=a[17]
    b.pi75=a[18]
    b.pi95=a[19]
    b.hewx=a[20]
    b.hewy=a[21]
    b.w90x=a[22]
    b.w90y=a[23]
    b.xpi=b.xpi[0:b.nx]
    b.xpr=b.xpr[0:b.nx]
    b.ypi=b.ypi[0:b.ny]
    b.ypr=b.ypr[0:b.ny]
    if bvar!=0:
        delstat= chi2.isf(0.1,4)
        derr=np.array([False,False,False,False])
        def xchisq(fpars):
            #xm=lorentzian(b.xpi,fpars)
            #xm=gaussian(b.xpi,fpars)
            xm=king_profile(b.xpi,fpars)
            return np.sum(np.square(b.xpr-xm)/bvar/b.ny)
        pval=np.amax(b.xpr[0:b.nx])
        spars=np.array([pval,b.med[0],b.hewx/2.,1.0])
        lpars=np.array([pval/2,b.cen[0]-b.hewx/2,b.hewx/2.,0.3])
        upars=np.array([pval*2,b.cen[0]+b.hewx/2,b.hewx*2.,3.0])
        b.fitx=srchmin(spars,lpars,upars,xchisq,delstat,derr)
        def ychisq(fpars):
            #ym=lorentzian(b.ypi,fpars)
            #ym=gaussian(b.ypi,fpars)
            ym=king_profile(b.ypi,fpars)
            return np.sum(np.square(b.ypr-ym)/bvar/b.nx)
        pval=np.amax(b.ypr[0:b.ny])
        spars=np.array([pval,b.med[1],b.hewy/2.,1.0])
        lpars=np.array([pval/2,b.cen[1]-b.hewy/2,b.hewy/2.,0.3])
        upars=np.array([pval*2,b.cen[1]+b.hewy/2,b.hewy*2.,3.0])
        b.fity=srchmin(spars,lpars,upars,ychisq,delstat,derr)
    else:
        b.fitx=False
        b.fity=False
    return b
def beam(arr,rbeam,blev,bvar):
    """Analysis of source above background within a circular beam

    Args:
        arr:     image array
        rbeam:   radius of beam in pixels
        blev:    average background level per pixel (to be subtracted)
        bvar:    variance on blev (-ve for counting statistics)

    Returns:
        | list with the following
        | **nsam**:    number of pixels in beam
        | **bflux**:   background in beam (e.g. counts)
        | **bsigma**:  standard deviation of background
        | **flux**:    source flux above background in beam (e.g. counts)
        | **fsigma**:  standard deviation of source flux
        | **peak**:    source x,y peak position
        | **cen**:     source x,y centroid position
        | **tha**:     angle (degrees) of major axis wrt x (x to y +ve)
        | **rmsa**:    source max rms width (major axis) (pixels)
        | **rmsb**:    source min rms width (minor axis) (pixels)
        | **fwhm**:    full width half maximum (pixels) about beam centre
        | **hew**:     half energy width (pixels) about beam centre
        | **w90**:     W90 (90% width) (pixels) about beam centre
        | **fwhmp**:   full width half maximum (pixels) about peak
        | **hewp**:    half energy width (pixels) about peak
        | **w90p**:    W90 (90% width) (pixels) about peak
        | **fwhmc**:   full width half maximum (pixels) about centroid
        | **hewc**:    half energy width (pixels) about centroid
        | **w90c**:    W90 (90% width) (pixels) about centroid
        | **fit**:     parameters from peak fit using peakchisq() 
    |
    | Fit performed if bvar!=0.
    | Parameters are saved in the list fit (see function srchmin())
    |   0:      peak value (no error range calculated)
    |   1:      peak X pixel position with 90% error range
    |   2:      peak Y pixel position with 90% error range
    |   3:      Lorentzian width including 90% upper and lower bounds
    |
    | The position of the beam is the current position within the image.
    | Use function setpos() to set the current position.
    """
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
    if bvar!=0:
        delstat= chi2.isf(0.1,4)
        iix=int(np.rint(b.cen[0]))
        iiy=int(np.rint(b.cen[1]))
        pval=arr[iix,iiy]
        spars=np.array([pval,b.cen[0],b.cen[1],b.fwhmc/2.])
        lpars=np.array([pval/2,b.cen[0]-b.fwhmc/2,
        b.cen[1]-b.fwhmc/2,b.fwhmc/2./2])
        upars=np.array([pval*2,b.cen[0]+b.fwhmc/2,
        b.cen[1]+b.fwhmc/2,b.fwhmc/2.*2])
        derr=np.array([False,True,True,True])
        b.fit=srchmin(spars,lpars,upars,peakchisq,delstat,derr)
    else:
        b.fit=False
    return b
def lecbeam(arr,s,h,g,blev,bvar,nt):
    """Analysis of source above background in a lobster eye cross-beam
   
    Args:
        arr:     image array
        s:       size of cross-beam square area in pixels
        h:       height of cross-arm quadrant in pixels (=2d/L)
        g:       fractional length of cross-arm (length=g*2d/L)
        blev:    average background level per pixel (to be subtracted)
        bvar:    variance on blev (-ve for counting statistics)
        nt:      dimension of output quadrant flux distribution in pixels
    
    Returns:
        list of the following
    | **qua**:     quadrant surface brightness distribution array nt by nt
    | **quan**:    quadrant  pixel occupancy array nt by nt
    | **nsam**:    number of pixels in beam
    | **bflux**:   background in beam (e.g. counts)
    | **bsigma**:  standard deviation of background
    | **flux**:    source flux above background in beam (e.g. counts)
    | **fsigma**:  standard deviation of source flux
    | **peak**:    source x,y peak position
    | **cen**:     source x,y centroid position
    | **hew**:     half energy width (pixels)
    | **w33**:     W33 33% width (pixels)
    | **ahe**:     half energy area (sq pixels)
    | **aw33**:    AW33 33% width area (sq pixels)
    | **fpeak**:   flux in peak pixel
    | **fit**:     results from fitting the quadrant distribution
    | **norm**:    fitted peak value
    | **G**:       fitted peak value Lorentzian width G pixels
    | **eta**      fitted ratio of cross-arms to peak
    | **alpha**    fitted index (1 for Lorentzian)
    |
    | The fitted parameters can also be found in fit.x
    | x[0] normalisation
    | x[1] width G pixels
    | x[2] eta brighness of cross-arms wrt spot
    | x[3] King index alpha 
    |
    | The position of the beam is the current position within the image.
    | Use function setpos() to set the current position.
    """
    a=imagesfor.qri_lecbeam(arr,s,h,blev,bvar,nt)
    b=bdata()
    b.qua=a[0]
    b.quan=a[1]
    b.nsam=a[2]
    b.bflux=a[3]
    b.bsigma=a[4]
    b.flux=a[5]
    b.fsigma=a[6]
    b.peak=a[7]
    b.cen=a[8]
    b.hew=a[9]
    b.w90=a[10]
    b.ahew=a[11]
    b.aw90=a[12]
    b.fpeak=a[13]
    def quafun(x,y,p):
        """Lobster eye quadrant using King function (modified Lorentzian)

        Args:
            x:          array of x values
            y:          array of x values
            p:          array of fitting parameters

        | par 0: normalisation (value at peak)
        | par 1: width of King profile G (not quite FWHM)
        | par 2: eta strength of cross-arm wrt central spot
        | par 3: power index alpha (1 for Lorentzian)

        Returns:
            function evaluated at x,y
        """
        eg=p[2]*p[1]/h
        x1=1/np.power(1+np.square(x*2.0/p[1]),p[3])
        x2=eg*(1-np.square(x/h))
        y1=1/np.power(1+np.square(y*2.0/p[1]),p[3])
        y2=eg*(1-np.square(y/h))
        return (x1*y1+x1*y2+x2*y1+x2*y2)*p[0]/np.square(1+eg)
    zqua=b.qua
    zqua.shape=(nt,nt)
    nqua=b.quan
    nqua.shape=(nt,nt)
    def quastat(p):
        stat=0.0
        for i in range(nt):
            xqua=i+0.5
            for j in range(nt):
                if nqua[i,j]>0:
                    yqua=j+0.5
                    mm=quafun(xqua,yqua,p)
                    stat=stat+np.square(zqua[i,j]-mm)/mm
        return stat
    delstat= chi2.isf(0.1,2)
    derr=np.array([False,False,False,False])
    spars=np.array([zqua[0,0],b.hew,1.0,1.0])
    lpars=np.array([zqua[0,0]/10,b.hew/10,0.1,0.1])
    upars=np.array([zqua[0,0]*10,b.hew*10,10.0,5.0])
    b.fit=srchmin(spars,lpars,upars,quastat,delstat,derr)
    b.norm=b.fit.x[0]
    b.G=b.fit.x[1]
    b.eta=b.fit.x[2]
    b.alpha=b.fit.x[3]
    return b
def lecimage(s,h,b,xcen,ycen,nx,ny):
    """Create an image array of the lobster eye cross-beam
    Args:
        s:      size of cross-beam square area in pixels
        h:      height of cross-arm triangle in pixels (=2d/L)
        b:      jwidth of cross-arm triangle in pixels
        xcen:   centre pixel (see coords below)
        ycen:   centre pixel (see coords below)
        nx:     first dimension of array
        ny:     second dimension of array
    
    Returns:
        image array

    Coordinate system for xcen,ycen is:
        | x runs from 0.0 on left to nx on right
        | Y runs from 0.0 on bottom to ny on top
        | Centre of bottom left pixel is therefore 0.5,0.5
        | Centre of top right pixel is nx-0.5,ny-0.5
    """
    a=imagesfor.qri_lecimage(s,h,b,xcen,ycen,nx,ny)
    return a
def lepsf(s,h,g,eta,alp,xcen,ycen,nx,ny):
    """Create an image of the lobster eye PSF

    Args:
        s:      size of cross-beam square area in pixels
        h:      height of cross-arm triangle in pixels (=2d/L)
        g:      width of Lorentzian central spot in pixels
        eta:    cross-arm to peak ratio at centre
        alp:    index (1 for Lorentzian)
        xcen:   centre of PSF (see below for coords. system)
        ycen:   centre of PSF (see below for coords. system)
        nx:     first dimension of array
        ny:     second dimension of array
    
    Returns:
        image array
    
    Coordinate system for xcen,ycen is:
        | x runs from 0.0 on left to nx on right
        | Y runs from 0.0 on bottom to ny on top
        | Centre of bottom left pixel is therefore 0.5,0.5
        | Centre of top right pixel is nx-0.5,ny-0.5
    """
    a=imagesfor.qri_lepsf(s,h,g,eta,alp,xcen,ycen,nx,ny)
    return a
def lebin(xe,ye,s,h,g,eta,nx,ny):
    """Create an image from an event list binning using lobster eye psf
    
    The image is effectively a cross-correlation of the event list with the lobster eye cross-beam.

    Args:
        xe:      x event positions, pixels
        ye:      y event positions, pixels
        s:      size of cross-beam square area in pixels
        h:      height of cross-arm triangle in pixels (=2d/L)
        g:      width of Lorentzian central spot in pixels
        eta:    cross-arm to peak ratio at centre
        nx:     first dimension of array
        ny:     second dimension of array
    
    Returns:
        image array
    
    |   Event positions assumed to run:
    |       x: 0 left edge to nx right edge
    |       y: 0 bottom edge to ny top edge
    |   Therefore centre of left bottom pixel is 0.5,0.5
    """
    a=imagesfor.qri_lebin(xe,ye,s,h,g,eta,nx,ny)
    return a
def setfield(nx,xleft,xright,ny,ybot,ytop):
    """Set local coordinates for image field

    Args:
        nx:     number of columns
        xleft:  local x left edge
        xright: local x right edge
        ny:    number of rows
        ybot:  local y bottom edge
        ytop:  local y top edge
    """
    imagesfor.qri_setfield(nx,xleft,xright,ny,ybot,ytop)

class iposition: pass
def getpos():
    """Get current position in image field

    Returns:
        list with following

    |   **pix**:     pixel position
    |   **xyl**:     local position
    |   **aes**:     local azimuth,elevation degrees
    |   **equ**:     Celestial RA,DEC degrees J2000
    |   **ecl**:     Ecliptic EA,EL degrees
    |   **gal**:     Galactic LII,BII degrees
    """
    a=imagesfor.qri_getpos()
    b=iposition()
    b.pix=a[0]
    b.xyl=a[1]
    b.aes=a[2]
    b.equ=a[3]
    b.ecl=a[4]
    b.gal=a[5]
    return b
def setpos(ipos,p):
    """Set current position in image field

    Args:
        ipos:   coordinate index
        p:      position coordinate pair

    ipos values:
        | 1: pixel 0-NCOLS, 0-NROWS
        | 2:  local X,Y
        | 3:  local azimuth,elevation degrees
        | 4:  Celestial RA,DEC degrees J2000
        | 5:  Ecliptic EA,EL degrees
        | 6:  Galactic LII,BII degrees
    
    Returns:
        current position using getpos()
    """
    imagesfor.qri_setpos(ipos,p)
    return getpos()
def setsky(xtodeg,ytodeg,ipr,mjd,ra,dec,roll):
    """Set up sky coordinates for image field

    Args:
        xtodeg: scale from local x to degrees (usually -ve)
        ytodeg: scale from local y to degrees
        ipr:    projection between local XY and local spherical
        mjd:    Modified Julian date
        ra:     Right Ascension (degrees J2000 at local origin)
        dec:    Declination (degrees J2000 at local origin)
        roll:   Roll angle (degrees from North to +ve elev. +ve clockwise)

    ipr values:
        |   0: Plate Carre
        |   1: Aitoff (Hammer)
        |   2: Lambert (equatorial aspect of Azimuthal equal-area projection)
    """
    imagesfor.qri_setsky(xtodeg,ytodeg,ipr,mjd,ra,dec,roll)
# Python QSOFT locator
def plt_show_locator(fg,npos):
    """Get local coordinate positions using the cursor

    Args:
        fg:     figure id returned by plt.figure()
        npos:   number of positions to be retured (clicked)
    
    Returns:
        xx,yy,  2 arrays containing npos local coordinates
    
    This routine provides a basic level of interaction with a displayed figure.
    It is used in place of a simple plt.show() call so that npos local coordinate
    positions can be selected interactively using the cursor from a displayed image
    and returned in arrays to the user.
    """
    return locator.plt_show_locator(fg,npos)
# Convert position to local XY coords
def toxy(ipos,p):
    """Convert position to local xy coordinates

    Args:
        ipos:   coordinate index

    | 1: pixel 0-NCOLS, 0-NROWS
    | 2: local X,Y
    | 3: local azimuth,elevation degrees
    | 4: Celestial RA,DEC degrees J2000
    | 5: Ecliptic EA,EL degrees
    | 6: Galactic LII,BII degrees
    
    Returns:
        position as local xy
    """
    n=len(p[:,0])
    pout=np.empty([n,2])
    for i in range(n):
        setpos(ipos,p[i])
        a=getpos()
        pout[i]=a.xyl
    return pout
# Plot a Hammer projection grid on figure pic
def hamgrid(pic):
    """Plot a Hammer projection grid on figure

    Args:
        pic:    figure object
    """
    n=45
    gap=10
    fixe=np.ones(n)*gap
    fixa=np.ones(n)*gap
    wa=360
    sama=np.arange(n)*wa/(n-1)-wa/2
    we=180
    same=np.arange(n)*we/(n-1)-we/2
    pp=np.empty([n,2])
    na=int(wa/2/gap)
    for i in range(-na,na):
        pp[:,0]=fixa*i
        pp[:,1]=same
        xy=toxy(3,pp)
        pic.plot(xy[:,0],xy[:,1],"w--")
    ne=int(we/2/gap)
    for i in range(-ne,ne):
        pp[:,0]=sama
        pp[:,1]=fixe*i
        xy=toxy(3,pp)
        pic.plot(xy[:,0],xy[:,1],"w--")
# Plot a Lambert projection grid on figure pic
def lamgrid(pic):
    """Plot a Lambert projection grid on figure
       
       Args:
           pic:    figure object
    """
    n=45
    gap=10
    fixe=np.ones(n)*gap
    fixa=np.ones(n)*gap
    wa=180
    sama=np.arange(n)*wa/(n-1)-wa/2
    we=180
    same=np.arange(n)*we/(n-1)-we/2
    pp=np.empty([n,2])
    na=int(wa/2/gap)
    for i in range(-na,na):
        pp[:,0]=fixa*i
        pp[:,1]=same
        xy=toxy(3,pp)
        pic.plot(xy[:,0],xy[:,1],"w--")
    ne=int(we/2/gap)
    for i in range(-ne,ne):
        pp[:,0]=sama
        pp[:,1]=fixe*i
        xy=toxy(3,pp)
        pic.plot(xy[:,0],xy[:,1],"w--")
