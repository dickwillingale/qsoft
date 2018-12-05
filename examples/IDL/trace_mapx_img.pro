PRO TRACE_MAPX_IMG
  
  ; useful vectors
  sn=[1,0,0]
  nn=[0,0,0]
  rx=[0,1,0]

  ; set focal length, radius of MCP aperture and slumping curvature
  f=50.
  rcur1=1e9
  ekev=1.
  psca= 1
  foc = 0 ;if 1 then focusing experiment, else, collimator

  ; define pores
  pitch=0.0258
  wall=0.0058
  pl = 0.5

  ;radii of modules
  rmm= REPLICATE(rcur1,1)
  ; length of pores in MPOs
  lm= REPLICATE(pl,1)
  ; Number of MPOs
  nd= N_ELEMENTS(lm)
  ; Set frame radius of curvature to mean
  rocf= 1E9
  
  ; MPO aperture dimension mm
  wplate= 39.
  wm= REPLICATE(wplate,nd)
  hm= wm
  ; Gap size between apertures mm
  gap= 2.
  ; Pore size mm
  ps= REPLICATE(pitch-wall,nd)
  ; loverd of pores
  loverd= lm/ps
  ; Open fraction
  openf= REPLICATE(0.6,nd)
  ; Wall thickness to give open fraction
  tm= REPLICATE(wall,nd)
  ; Pitch of pores
  pmm= REPLICATE(pitch,nd)
  ; Multifibre size
  fm= (ps+tm)*35.
  ; Rotation of plates
  am= REPLICATE(0,nd)
  ; Positions of plates
  xm= 0;[1.]*wm/2.
  ym= 0;[1.]*wm/2.
  ; half width of full aperture
  hwid= 40./2.

  ; set up random arrays of bias angles rms 0.25 arcmins
  ;bias= 1.*!dpi/180.0/60.0
  seed = 1.234
  ;bu=RANDOMN(seed,nd)*bias
  ;bv=RANDOMN(seed,nd)*bias
  bu=[5.]*!dpi/180.0/60.0
  bv=[4.]*!dpi/180.0/60.0
  
  ; set up deformations etc. for the plates
  xd= LINDGEN(nd)
  yd= 1
  ; ds is scaling factor for intrinsic slumping errors
  ds=1
  zx=REPLICATE(ds,nd)
  ; de is the maximum thermoelelastic pore axial pointing error radians
  ; If -ve then fixed pattern tilt
  de=  -1.0*!dpi/180.0/60.0
  zy=REPLICATE(de,nd)
  ; trms is the rms pore rotation error about pore axis in radians
  trms= 0.0*!dpi/180
  za=REPLICATE(trms,nd)
  ; dshr is the shear distance for each pore mm used in multifibre model
  ; of the shear error
  dshr= ps(0)*0.01
  zd=RANDOMN(seed*2.0,nd)*0.0+dshr
  ; qrms +ve rms Gaussian tilt/figure errors
  ; qrms -ve Cauchy (Lorentzian) tilt/figure errors 2*qrms=FWHM
  qrms= 2.*!dpi/180.0/60.0
  fd=REPLICATE(qrms,nd)
  ;
  ; reflecting surface definition
  rough= 11.;13.0
  rind=1.6
  brk= 5;350.0
  gam= rind-1.0
  nrough= rough^2/((1+1/gam)*brk^(-gam))
  rho=22.65
  mspec="Ir"
  sq= REPLICATE(1,nd)
  ; spare parameter
  spar= REPLICATE(0,nd)
  ;
  plt= [xm,ym,am,wm,hm,lm,rmm,fm,pmm,tm,sq,bu,bv,spar,spar]
  
  idef=1
  hdet= 19.2/2
  np= 512.
  dpix= 2.*hdet/np
  
  
  dlim=[-hdet,-hdet,hdet,hdet]
  hmxt= hdet/2
  ; surface optical properties
  xop=QRT_XOPT(mspec,rho,ekev,1)
  is=1
  ; set up source
  ccy=SIN(0.0*!dpi/180)
  ccz=SIN(0.0*!dpi/180)
  di=[-SQRT(1.0-ccy^2-ccz^2),ccy,ccz]
  slim=[-hwid,-hwid,hwid,hwid]
  nray=500000
  
    psrc=[f*2.,0,0]
    ; set position of centre of square pore MCP array
    pmcp=[f,0,0]
    pap=[f*1.03,0,0]

 
  alim=[-100.,-1.,100.,1.]
  ; rotate detector axes wrt the plate reference axes
  thetadet=-0.*!dpi/180.
  rxdet=[0,cos(thetadet),sin(thetadet)]
  ;
  qrt_q=QRT_RESET()
  qrt_q=QRT_SURFACE(is,1,ekev,nrough,brk,rind,xop.alpha,xop.gamma,0,0,0,0,0)
 ; qrt_q=QRT_SOURCE(4,di,psrc,pap,sn,rx,slim,0,nray,0)

  im = FLTARR(100,100)
  im(45:55,40:60) = 25.
  xsam = FINDGEN(100)
  ysam = FINDGEN(100)
  zsam = im
  deff = 2
  qrt_q=QRT_DEFORM(deff,1,1,100,100)
  qrt_q=QRT_DEFMAT(deff,1,xsam,ysam,zsam)
  
 qrt_q=QRT_SOURCE(4,di,psrc,pap,sn,rx,slim,0,nray,deff)
  qrt_q=QRT_DEFORM(idef,1,5,nd,1)
  qrt_q=QRT_DEFMAT(idef,1,xd,yd,zx)
  qrt_q=QRT_DEFMAT(idef,2,xd,yd,zy)
  qrt_q=QRT_DEFMAT(idef,3,xd,yd,za)
  qrt_q=QRT_DEFMAT(idef,4,xd,yd,zd)
  qrt_q=QRT_DEFMAT(idef,5,xd,yd,fd)
  qrt_q=QRT_SQMPOARR(pmcp,sn,rx,rocf,hwid,idef,plt)
  qrt_q=QRT_DETECTOR(2,nn,sn,rxdet,dlim,0)
  ;qrt_list()
  results=QRT_TRACE(0,2,-1)
  ; get detected positions
  detpos=READ_QRT('detected.dat',/det)
  ; create image
  a=QRI_BINXY(detpos.ydet,detpos.zdet,0,detpos.area,-hdet,hdet,np,-hdet,hdet,np)

  ;pdf(paste(model,"_psf_",ekev,".pdf",sep=""))
  a.xlab="mm"
  a.ylab="mm"
  cgLOADCT, 1,/reverse
  WINDOW,2,xsize=768,ysize=768
  ;convert to cm^2
  a.data_array = a.data_array/100.
  a.zlim = a.zlim/100.
  PRINT, TOTAL(a.data_array), ' = rough effarea cm^2" 
  ;dpix/10 is pixel size in cm
  ;a.data_array=a.data_array/(dpix/10.)^2; normalise pixel size i.e. cm^2(collectingarea)/cm^2(pixel area in focal plane)
  ;a.zlim(1) = MAX(a.data_array)
  QRI_ZOOMIMAGE,a,[-hdet,hdet],[-hdet,hdet],[0,1.5E-3];a.zlim[1]/psca]

  ;lines(c(-hmxt,hmxt,hmxt,-hmxt,-hmxt),c(-hmxt,-hmxt,hmxt,hmxt,-hmxt),col="red")
  ;title(main=paste(model,"Ray Tracing Simulation",ekev,"keV"))
  ;dev.off()
  QRI_s=QRI_SETFIELD(np,-hdet,hdet,np,-hdet,hdet)
  QRI_sp=QRI_SETPOS(2,[0,0])
  ; beam size of plate
  rpix=hdet/dpix
;  b=QRI_BEAM(a.data_array,rpix,0,-1)


RETURN
END