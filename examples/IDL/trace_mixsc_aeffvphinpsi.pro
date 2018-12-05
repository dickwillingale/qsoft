PRO trace_mixsc_aEffvPHInPSI
  
  ; useful vectors
  sn=[1,0,0]
  nn=[0,0,0]
  rx=[0,1,0]

  ; set focal length, radius of MCP aperture and slumping curvature
  f=550./2.
  rcur1=f*2.
  ekev = 4.09
  psca= 1
  foc = 0 ;if 1 then focusing experiment, else, collimator

  ; define pores
  pitch=0.0258
  wall=0.0058
  pl = 1.2

  ;radii of modules
  rmm= REPLICATE(rcur1,4)
  ; length of pores in MPOs
  lm= REPLICATE(pl,4)
  ; Number of MPOs
  nd= N_ELEMENTS(lm)
  ; Set frame radius of curvature to mean
  rocf= 550.
  
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
  openf= REPLICATE(0.6,4)
  ; Wall thickness to give open fraction
  tm= REPLICATE(wall,4)
  ; Pitch of pores
  pmm= REPLICATE(pitch,4)
  ; Multifibre size
  fm= (ps+tm)*35.
  ; Rotation of plates
  am= REPLICATE(0,nd)
  ; Positions of plates
  xm= [1.,1.,-1.,-1.]*wm/2. + [1.,1.,-1.,-1.]
  ym= [1.,-1.,-1.,1.]*wm/2. + [1.,-1.,-1.,1.]
  ; half width of full aperture
  hwid= 40.

  ; set up random arrays of bias angles rms 0.25 arcmins
  ;bias= 1.*!dpi/180.0/60.0
  seed = 1.234
  ;bu=RANDOMN(seed,nd)*bias
  ;bv=RANDOMN(seed,nd)*bias
  bu=[5.,5.0,-4,-4.]*!dpi/180.0/60.0
  bv=[4.,-4.,-3.,4.]*!dpi/180.0/60.0
  
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
  rough= 11.0
  rind=1.6
  brk= 5.0
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
  hdet= 19.2/2.
  np= 512.
  dpix= 2.*hdet/np
  
  
  dlim=[-hdet,-hdet,hdet,hdet]
  hmxt= hdet/2
  is=1

  slim=[-hwid,-hwid,hwid,hwid]
  nray=100000
  
  IF foc EQ 1 THEN BEGIN
    psrc=[f*1.05,0,0]
    ; set position of centre of square pore MCP array
    pmcp=[f,0,0]
    pap=[f*1.03,0,0]
  ENDIF ELSE BEGIN
     psrc=[rcur1*1.05,0,0]
     ; set position of centre of square pore MCP array
     pmcp=[rcur1,0,0]
     pap=[rcur1*1.03,0,0]
  ENDELSE
 
  alim=[-100.,-1.,100.,1.]
  ; rotate detector axes wrt the plate reference axes
  thetadet=-0.*!dpi/180.
  rxdet=[0,cos(thetadet),sin(thetadet)]
  ;
  Arange = 12.
  PhiAngs = (FINDGEN(50)+1)*Arange/50. - Arange/2.
  PsiAngs = (FINDGEN(50)+1)*Arange/50. - Arange/2.
  Aeff = FLTARR(N_ELEMENTS(phiAngs), N_ELEMENTS(phiAngs))
  
  FOR phiC = 0,N_ELEMENTS(phiAngs) - 1 DO BEGIN
    phi = PhiAngs(phic)
    IF Phi MOD 10 EQ 0 THEN PRINT, Phi
    FOR psiC = 0,N_ELEMENTS(psiAngs) - 1 DO BEGIN
        psi = PsiAngs(psic)
        ; set up source
        ccy=SIN(phi*!dpi/180)
        ccz=SIN(psi*!dpi/180)
        di=[-SQRT(1.0-ccy^2-ccz^2),ccy,ccz]

        qrt_q=QRT_RESET()
        ; surface optical properties
        xop=QRT_XOPT(mspec,rho,ekev,1)
        qrt_q=QRT_SURFACE(is,1,ekev,nrough,brk,rind,xop.alpha,xop.gamma,0,0,0,0,0)
        qrt_q=QRT_SOURCE(2,di,nn,psrc,sn,rx,slim,0,nray,0)
        qrt_q=QRT_DEFORM(idef,1,5,nd,1)
        qrt_q=QRT_DEFMAT(idef,1,xd,yd,zx)
        qrt_q=QRT_DEFMAT(idef,2,xd,yd,zy)
        qrt_q=QRT_DEFMAT(idef,3,xd,yd,za)
        qrt_q=QRT_DEFMAT(idef,4,xd,yd,zd)
        qrt_q=QRT_DEFMAT(idef,5,xd,yd,fd)
        qrt_q=QRT_APERTURE(4,0,pap,sn,[0,SIN(!DPI/4.),COS(!DPI/4.)],alim,0)
        qrt_q=QRT_APERTURE(4,0,pap,sn,[0,-SIN(!DPI/4.),COS(!DPI/4.)],alim,0)
        qrt_q=QRT_APERTURE(6,0,pap,sn,rx,[4.1,100.,0,2.*!DPI],0)
        qrt_q=QRT_SQMPOARR(pmcp,sn,rx,rocf,hwid,idef,plt)
        qrt_q=QRT_DETECTOR(2,nn,sn,rxdet,dlim,0)
        ;qrt_list()
        results=QRT_TRACE(0,2,-1)
        ; get detected positions
        detpos=READ_QRT('detected.dat',/det)
        ; create image
        a=QRI_BINXY(detpos.ydet,detpos.zdet,0,detpos.area/100.,-hdet,hdet,np,-hdet,hdet,np)
        Aeff(phic,psic) = TOTAL(a.data_array)
    ENDFOR
  ENDFOR

SAVE, PhiAngs, PsiAngs, Aeff, FILENAME='/Users/Adrian/Desktop/mixs/MIXSCAeffVsPhiAndPsi.sav'
cgLOADCT, 13
surface_area = SURFACE(Aeff,PhiAngs, PsiAngs, TEXTURE_IMAGE=Aeff)

RETURN
END