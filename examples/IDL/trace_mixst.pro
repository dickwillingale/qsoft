PRO TRACE_MIXST

; useful vectors
  sn=[1,0,0]
  nn=[0,0,0]
  rx=[0,1,0]

; set focal length, radius of MCP aperture and slumping curvature
  f=1000.
  rmcp=110.
  rcur1=f*4.
  rcur2=f*4./3.

; define pores
  pitch=0.0258
  wall=0.0058
  ipack=6
  fb=1.902

; the following is ignored for ipack=6 but need to be set
  pl=1

; parameters for MIXS aperture
  nsec1=6
  nsec2=12
  nsec3=18
  gap=4.0

; Common parameters for front and back plates
  rhi1=41.0
  rhi2=71.0
  rhi3=105.0
  rlo1=20.0
  rlo2=49.0
  rlo3=79.0
  pl1=2.3
  pl2=1.3
  pl3=0.9

; define the large arrayd for qrt_radmpoarr
; widht in radial dimension
  ww1=rhi1-rlo1
  ww2=rhi2-rlo2
  ww3=rhi3-rlo3
;centre in radial dimansion
  rr1=(rhi1+rlo1)/2.
  rr2=(rhi2+rlo2)/2.
  rr3=(rhi3+rlo3)/2.
;azimuthal extent of each sector
  dd1=2.*!dpi/nsec1
  dd2=2.*!dpi/nsec2
  dd3=2.*!dpi/nsec3

;input surface quality indices
  iq1 = 1
  iq2 = 2
  
  ;concatenate parameters for inner middle and outer rings into single array
  
  XP =  [REPLICATE(rr1,nsec1), REPLICATE(rr2,nsec2), REPLICATE(rr3,nsec3)] ;input x position of each MPO (R for mixs geometry)
  YP =  [(FINDGEN(nsec1))*dd1-!dpi+0.5*dd1, (FINDGEN(nsec2))*dd2-!dpi+0.5*dd2, (FINDGEN(nsec3))*dd3-!dpi+0.5*dd3] ;input y position of each MPO
  TC =  REPLICATE(0.0,nsec1+nsec2+nsec3) ;input rotation angle of each MPO
  WC =  [REPLICATE(ww1,nsec1), REPLICATE(ww2,nsec2), REPLICATE(ww3,nsec3)] ;input width of each MPO delx  OR radial extent of sector for radial MPO array (in mm)
  HC =  [REPLICATE(dd1,nsec1), REPLICATE(dd2,nsec2), REPLICATE(dd3,nsec3)] ;input height of each MPO dely OR azimuthal extent of sector for radial MPO array (in radians)
  AC =  [REPLICATE(pl1,nsec1), REPLICATE(pl2,nsec2), REPLICATE(pl3,nsec3)] ;input axial length of each MPO (thickness)
  CU1 =  REPLICATE(rcur1,nsec1+nsec2+nsec3) ;input radius of curvature of each MPO
  CU2 =  REPLICATE(rcur2,nsec1+nsec2+nsec3) ;input radius of curvature of each MPO
  MF =  REPLICATE(fb,nsec1+nsec2+nsec3) ;input size of multifibres in MPO
  PP =  REPLICATE(pitch,nsec1+nsec2+nsec3) ;input pitch of pores (pore size + wall thickness)
  WA =  REPLICATE(wall,nsec1+nsec2+nsec3) ;input wall thickness between pores
  SQ =  [REPLICATE(iq1,nsec1), REPLICATE(iq2,nsec2), REPLICATE(iq2,nsec3)] ;input reflecting surface quality index for pores
  BU =  REPLICATE(0.0,nsec1+nsec2+nsec3);input bias angle x radians
  BV =  REPLICATE(0.0,nsec1+nsec2+nsec3);input bias angle y radians
  SPar =  REPLICATE(0.0,nsec1+nsec2+nsec3);input spare parameter!
 ; print, hc
  nmod = nsec1+nsec2+nsec3


  a1=[XP,YP,TC,WC,HC,AC,CU1,MF,PP,WA,SQ,BU,BV,SPar,SPar]
  a2=[XP,YP,TC,WC,HC,AC,CU2,MF,PP,WA,SQ,BU,BV,SPar,SPar]

xm=xp*COS(YP)
ym=xp*SIN(YP)
wm=wc
loverd=ac/(pp-wa)

IF 0 THEN BEGIN
  hwid = (rr3(0)+wm(-1)/2.)*2.
  graphic=PLOT((FINDGEN(100)-50)*hwid/100, (FINDGEN(100)-50)*hwid/100,linestyle=6,aspect_ratio=1,name='aa')
    
  FOR i = 0,N_ELEMENTS(xp)-1 DO BEGIN
    polyg = QR_SECTORS(xp(i),yp(i),wm(i),hc(i),spar(i),color='green',wind='aa')
    asd=TEXT(xm(i),ym(i)+8,STRTRIM(i,2),/data,alignment=0.5,font_size=6)         
    blurb= STRTRIM(sq(i),2)
    asdf=TEXT(xm(i),ym(i),blurb,/data,alignment=0.5,font_size=6)
    blurb= STRMID(ac(i),STRPOS(ac(i),'.')-1,4) + ' ' + STRMID(loverd(i),STRPOS(loverd(i),'.')-3,5)
    asdfg=TEXT(xm(i),ym(i)-8,blurb,/data,alignment=0.5,font_size=6)
  ENDFOR
ENDIF

;Position of front plate is defined by slump radius and inner ring thickness.
;adjust radius of curvature of front plate by thickness of inner annulus
  rcur1=rcur1;+pl1

;Offset of entire front optic towards source to give a larger interplate gap
  ofmcp1=2.5

; positions of MCP apertures
  pmcp1=[pl1+ofmcp1,0,0]
  pmcp2=[0,0,0]

; define aperture for source illumination
  sapin=rlo1-2.0
  sapout=rmcp+2.0
  spos=[27000.0,0,0]
  apos=[pl1+ofmcp1+2.0,0,0]
  slim=[sapin,sapout]

;either specify number of rays over aperture or area of aperture per ray
  nray=5000000
  apry=0.0
 ;nray=0
 ;apry=0.02
  stype=3

  ; set detector
  dhw=rmcp+10.
  pix=0.2
  par=pix/f
  pam=par*180.0/!pi*60.
  dlim=[-dhw,-dhw,dhw,dhw]
  rbeam=5.
  npix=FIX(dhw*2./pix)
  dtype=2
  pdet=[-1000.,0,0]
  
; aperture definition
;  th1=270.*!dpi/180.
;  th2=360.*!dpi/180.
;  ra1=rlo2
;  ra2=rhi2
;  radlim=[ra1,ra2,th1,th2]
;  apos=[pl1+ofmcp1+1.0,0,0]
  
; RAY TRACING BIT
arange = 0.5 ;degrees
offang = DINDGEN(6)*arange/5.

; open window and set !p.multi for multiple plots
cgDISPLAY,768,768
;!P.multi=[0,3,2]

FOR i = 0,0 DO BEGIN

  a=QRT_RESET()
  b=QR_RSEED(5)
    
; define source & geometry
  ccy=-SIN(offang(i)*!dpi/180.)
  ccz=0.0
  ccx=SQRT(1.0-ccy^2-ccz^2)
  di=[-ccx,-ccy,ccz]

    c=QRT_SOURCE(1,di,spos,apos,sn,rx,slim,apry,nray,0) 

  ekev=1.5
  rough= 13.0
  rind=1.4
  brk= 350.0
  gam= rind-1.0
  nrough= rough^2/((1+1/gam)*brk^(-gam))
  nrough = 0
    Rf=QRT_XOPT('Ir',22.42,ekev,1) 
    Rf2=QRT_XOPT('Si40 O98 K8 Na7 Pb6 Bi2',3.3,ekev,1); input bare glass

    d=QRT_SURFACE(iq1,1,ekev,nrough,brk,rind,Rf.alpha,Rf.gamma,0,0,0,0,0)
    dd=QRT_SURFACE(iq2,1,ekev,nrough,brk,rind,Rf2.alpha,Rf2.gamma,0,0,0,0,0)
    
    idef=0
    ; define radial aperture blocks between sectors.
    pap = [pl1+ofmcp1+1.,0,0]
    ;sextants
    alim = [-rmcp*1.2,-gap/2.,rmcp*1.2,gap/2.]
    FOR gg= 0,5 DO qrt_q=QRT_APERTURE(4,0,pap,sn,[0,SIN(gg*2.*!DPI/6.),COS(gg*2.*!DPI/6.)],alim,0)
    ;Middle sectors
    alim = [rhi1+2.,-gap/2.,rhi2+2.,gap/2.]
    FOR gg= 0,11 DO qrt_q=QRT_APERTURE(4,0,pap,sn,[0,SIN(gg*2.*!DPI/12.),COS(gg*2.*!DPI/12.)],alim,0)
    ;Outer Sectors
    alim = [rhi2+2.,-gap/2.,rhi3+2.,gap/2.]
    FOR gg= 0,17 DO qrt_q=QRT_APERTURE(4,0,pap,sn,[0,SIN(gg*2.*!DPI/18.),COS(gg*2.*!DPI/18.)],alim,0)    

     ff=QRT_APERTURE(6,0,pap,sn,[0,0,1],[20,110,dd1,2.*dd1],1)   
    e=QRT_RADMPOARR(pmcp1,sn,[0,0,1],rcur1,rmcp,idef,1,a1)
    f=QRT_RADMPOARR(pmcp2,sn,[0,0,1],rcur2,rmcp,idef,2,a2)
    


    g=QRT_DETECTOR(dtype,pdet,sn,rx,dlim,0)
;ll = QRT_LIST()
; trace rays
    results=QRT_TRACE(0,dhw*2.,-2)
    
; get detected positions and rays
  detpos = READ_QRT('detected.dat',/det,/traced)

; create image
  ihw = 55
  npix = FIX(ihw/0.05)   
  h = QRI_BINXY(detpos.YDET,detpos.ZDET,0,detpos.AREA,-ihw,ihw,npix,-ihw,ihw,npix)
  ;convert the image from mm^2 to cm^2 & set plot axes titles
  h.data_array = h.data_array/100.
  h.zlim = h.zlim/100.
  h.xlab ='y - mm'
  h.ylab = 'z - mm'
  h.title = 'Angle = ' + STRMID(offang(i),STRPOS(OFFANG(I),'.')-1,3) + ' Degrees'

; set color table limits /100 for mm^2 to cm^2 and 10000 
; to increase contrast by saturating the color table  
IF i EQ 0 THEN ctmax = h.zlim/(100.)
  QRI_ZOOMIMAGE,h,h.xlim,h.ylim,ctmax
ENDFOR
dp = detpos
dispRays = QRI_DISPLAYRAYS(dp,100,0,0,0, lims=[-1100,200,-120,120])
dp=detpos
refMap = QRI_REFMAP(dp,500.,500., lims=[-rmcp,rmcp,-rmcp,rmcp])

;reset !P.multi 
!P.multi=0
stop
RETURN
END