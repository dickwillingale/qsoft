# X-ray sequential ray tracing
import xsrtfor 
import numpy as np
# Initialisation and reseting
def init():
    xsrtfor.qrt_init()
def reset():
    xsrtfor.qrt_init()
init()
def rseed(iseed):
    xsrtfor.qr_rseed(iseed)
# Proton/electron path tracing
class field: pass
def bfield(dm,pdx,pdy,pdz,ddx,ddy,ddz,px,py,pz):
    a=xsrtfor.qrt_bfield(dm,pdx,pdy,pdz,ddx,ddy,ddz,px,py,pz)
    b=field()
    dx=a[0]
    dy=a[1]
    dz=a[2]
    b.bf=np.sqrt(dx**2+dy**2+dz**2)
    b.dx=dx/b.bf
    b.dy=dy/b.bf
    b.dz=dz/b.bf
    b.rmin=a[3]
    return b
class eptrace: pass
def prtathena(dm,pdx,pdy,pdz,ddx,ddy,ddz,ekv,tsig,xaper,xdiv,
    rrings,trings,drings,rdet,maxst):
    a=xsrtfor.qrt_prtathena(dm,pdx,pdy,pdz,ddx,ddy,ddz,ekv,tsig,xaper,xdiv,
    rrings,trings,drings,rdet,maxst)
    b=eptrace()
    b.npath=a[3]
    b.iqual=a[4]
    b.xp=a[0][:b.npath]
    b.yp=a[1][:b.npath]
    b.zp=a[2][:b.npath]
    return b
def eltmxt(dm,pdx,pdy,pdz,ddx,ddy,ddz,ekv,tsig,xaper,xdiv,
    apy,apz,apsiz,wdiv,ddiv,sdet,maxst):
    a=xsrtfor.qrt_eltmxt(dm,pdx,pdy,pdz,ddx,ddy,ddz,ekv,tsig,xaper,xdiv,
    apy,apz,apsiz,wdiv,ddiv,sdet,maxst)
    b=eptrace()
    b.npath=a[3]
    b.iqual=a[4]
    b.xp=a[0][:b.npath]
    b.yp=a[1][:b.npath]
    b.zp=a[2][:b.npath]
    return b
# Sequential Ray Tracing
class srtdat: pass
def srtreset():
    xsrtfor.qrt_init()
def srtlist():
    xsrtfor.qrt_list()
def shift(iss,pl):
    xsrtfor.qrt_shift(iss,pl)
def rotate(iss,pl,ax,angle):
    xsrtfor.qrt_rotate(iss,pl,ax,angle)
def kbs(pcen,pnor,raxi,ipack,rmin,rmax,flen,csize,pitch,wall,
    plmin,plmax,idf,iq):
    nmax=2000
    a=xsrtfor.qrt_kbs(pcen,pnor,raxi,ipack,rmin,rmax,flen,csize,
    pitch,wall,plmin,plmax,idf,iq,nmax)
    return a
def surface(iss,it,ekev,srgh,fmin,pind,alpha,gamma,
    angs,refs,gpitch,dhub,order):
    nref=len(angs)
    xsrtfor.qrt_surface(iss,it,ekev,srgh,fmin,pind,alpha,gamma,
    nref,angs,refs,gpitch,dhub,order)
def mirror(idd,idf,iq,anml,arfx,apos,alim,nsurf):
    xsrtfor.qrt_mirror(idd,idf,iq,anml,arfx,apos,alim,nsurf)
def moa(pcen,pno,rax,rcur,xyap,pitch,wall,plen,idf,iq):
    xsrtfor.qrt_moa(pcen,pno,rax,rcur,xyap,pitch,wall,plen,idf,iq)
def baffle(xmin,xmax,rad,ax,ar,rp,iq):
    xsrtfor.qrt_baffle(xmin,xmax,rad,ax,ar,rp,iq)
def lens(idd,idf,iq,anml,arfx,apos,rap,r1,r2,refind,thickq):
    xsrtfor.qrt_lens(idd,idf,iq,anml,arfx,apos,rap,r1,r2,refind,thickq)
def prism(idd,idf,iq,anml,arfx,apos,rap,d1,d2,refind,thick):
    xsrtfor.qrt_prism(idd,idf,iq,anml,arfx,apos,rap,d1,d2,refind,thick)
def opgrat(idd,defi,iq,al,fpos,zpos,gpos,alim):
    nlim=len(alim)
    a=qrt_opgrat(idd,defi,iq,al,fpos,zpos,gpos,nlim,alim)
    return a
def c1nest(ax,ar,ff,iq,ib,idf,pl,ph,hl,hh,rpl,rph,rhl,rhh,tin,tj,tout):
    wrk=np.array([(len(rpl)-1)*2+12])
    xsrtfor.qrt_c1nest(ax,ar,ff,iq,ib,idf,pl,ph,hl,hh,rpl,rph,rhl,rhh,tin,tj,
    tout,wrk)
def w1nest(xj,rj,ra,pl,ph,hl,hh,tin,tj,tout,ax,ar,ff,defi,iq,ib):
    xsrtfor.qrt_w1nest(xj,rj,ra,pl,ph,hl,hh,tin,tj,tout,ax,ar,ff,defi,iq,ib)
def wolter2(rp,gp,rh,gh,rm,fovr,ax,ar,ff,idf,iq):
    xsrtfor.qrt_wolter2(rp,gp,rh,gh,rm,fovr,ax,ar,ff,idf,iq)
def spider(cone,apos,anml,arfx,nsec,cwid,awid):
    xsrtfor.qrt_spider(cone,apos,anml,arfx,nsec,cwid,awid)
def trace(ideb,riris,iopt):
    a=xsrtfor.qrt_trace(ideb,riris,iopt)
    b=srtdat()
    b.area=a[0]
    b.dshft=a[1]
    b.ybar=a[2]
    b.zbar=a[3]
    b.rms=a[4]
    return b
def source(it,sd,sp,ap,an,ar,al,apry,nr,idef):
    xsrtfor.qrt_source(it,sd,sp,ap,an,ar,al,apry,nr,idef)
def detector(idd,dpos,dnml,drfx,dlim,radet):
    xsrtfor.qrt_detector(idd,dpos,dnml,drfx,dlim,radet)
def fresnel(nreal,kimag,angs):
    a=xsrtfor.qrt_fresnel(nreal,kimag,angs)
    b=srtdat()
    b.rs=a[0]
    b.rp=a[1]
    b.runp=a[2]
    return b
def deform(idd,it,nm,nx,ny):
    xsrtfor.qrt_defs(idd,it,nm,nx,ny)
def defmat(idd,im,xsam,ysam,zdef):
    nw=len(xsam)*len(ysam)
    w1=np.empty(nw,dtype=float)
    w2=np.empty(nw,dtype=float)
    xsrtfor.qrt_mats(idd,im,xsam,ysam,zdef,w1,w2)
def sqpore(pcen,pnorm,raxis,rcur,ipack,rap,pitch,wall,plen,idf,
    iq,plmin,plmax,fibre,ar):
    xsrtfor.qrt_sqpore(pcen,pnorm,raxis,rcur,ipack,rap,pitch,wall,plen,
    idf,iq,plmin,plmax,fibre,ar)
def sqmpoarr(pcen,pnorm,raxis,rcur,hwid,idf,ar):
    xsrtfor.qrt_sqmpoarr(pcen,pnorm,raxis,rcur,hwid,idf,ar)
def sipore(pcen,pnorm,raxis,flen_bn,rpitch,apitch,wall,
    rm,pm,tm,wm,hm,am,cm,gm,wfr,a2j,idf,iq):
    xsrtfor.qrt_sipore(pcen,pnorm,raxis,flen_bn,rpitch,apitch,wall,
    rm,pm,tm,wm,hm,am,cm,gm,wfr,a2j,idf,iq)
def aperture(idd,idf,ap,an,ar,alim,nsurf):
    xsrtfor.qrt_aperture(idd,idf,ap,an,ar,alim,nsurf)
def elips(org,axs,cen,xmin,xmax,amin,amax,smb,rab,ide,iq):
    xsrtfor.qrt_elips(org,axs,cen,xmin,xmax,amin,amax,smb,rab,ide,iq)
