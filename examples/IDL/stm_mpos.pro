PRO STM_MPOS


  ; Set up array of Ir coated 40 micron STM plates
  ppp= READ_ASCII("/Users/Adrian/IDLWorkspace/am136/qIDL/examples/stm_mpos.dat",header=head, $
          data_start=1, template=ASCII_TEMPLATE("/Users/Adrian/IDLWorkspace/am136/qIDL/examples/stm_mpos.dat"))
  pp = {plate:ppp.field1,ix:ppp.field2, iy:ppp.field3, rocp:ppp.field4, rocx:ppp.field5,thick:ppp.field6,openf:ppp.field7}
  plab= pp.plate
  ; radius of curvature of MPOs
  rm= pp.rocx
  PRINT, "standard deviation of RoC",stddev(rm)
    ; length of pores in MPOs
  lm= pp.thick
  ; Number of MPOs
  nd= N_ELEMENTS(lm)
  ; Set frame radius of curvature to mean
  rocf= MEAN(rm)
  PRINT,"mean RoC",rocf
  
  ; MPO aperture dimension mm
  wplate= 38.
  wm= REPLICATE(wplate,nd)
  hm= wm
  ; Gap size between apertures mm
  gap= 4
  ; Pore size mm
  ps= REPLICATE(0.04,nd)
  ; loverd of pores
  loverd= lm/ps
  ; Open fraction
  openf= pp.openf
  ; Wall thickness to give open fraction
  tm= ps/SQRT(openf)-ps
  ; Pitch of pores
  pm= ps+tm
  ; Multifibre size
  fm= (ps+tm)*25.
  ; Rotation of plates
  am= REPLICATE(0,nd)
  ; Positions of plates
  ix= pp.ix
  iy= pp.iy
  ; half width of full aperture
  hwid= 110.
  ;
  xm=DBLARR(nd)
  ym=DBLARR(nd)
  FOR i = 0,nd-1 DO BEGIN
      xm(i)= (ix(i)-3)*(wplate+gap)
      ym(i)= (iy(i)-3)*(wplate+gap)
  ENDFOR

; set up random arrays of bias angles rms 0.25 arcmins
bias= 0.25*!dpi/180.0/60.0
seed = 1.234
bu=RANDOMN(seed,nd)*bias
bv=RANDOMN(seed,nd)*bias

; set up deformations etc. for the plates
xd= LINDGEN(nd)
yd= 1
; ds is scaling factor for intrinsic slumping errors
ds=1
zx=REPLICATE(ds,nd)
; de is the maximum thermoelelastic pore axial pointing error radians
; If -ve then fixed pattern tilt
de=  -5.0*!dpi/180.0/60.0
zy=REPLICATE(de,nd)
; trms is the rms pore rotation error about pore axis in radians
trms= 0.0*!dpi/180
za=REPLICATE(trms,nd)
; dshr is the shear distance for each pore mm used in multifibre model
; of the shear error
dshr= ps(0)*0.02
zd=RANDOMN(seed*2.0,nd)*0.0+dshr
; qrms +ve rms Gaussian tilt/figure errors
; qrms -ve Cauchy (Lorentzian) tilt/figure errors 2*qrms=FWHM
qrms= -1.0*!dpi/180.0/60.0
fd=REPLICATE(qrms,nd)
;
; reflecting surface definition
rough= 11.0
rind=1.4
brk= 10.0
gam= rind-1.0
nrough= rough^2/((1+1/gam)*brk^(-gam))
rho=22.65
mspec="Ir"
sq= REPLICATE(1,nd)
; spare parameter
sp= REPLICATE(0,nd)
;
plt= [xm,ym,am,wm,hm,lm,rm,fm,pm,tm,sq,bu,bv,sp,sp]

; IDL save file
SAVE,rocf,hwid,nd,xd,yd,zx,zy,za,zd,fd,plt,nrough,rind,brk,mspec,rho,plab,filename="/Users/Adrian/IDLWorkspace/am136/qIDL/examples/stm_mpos.sav"

graphic=PLOT((FINDGEN(100)-50)*hwid*2/100, (FINDGEN(100)-50)*hwid*2/100,linestyle=6,aspect_ratio=1)
;qr_rectangles(xm[1:nd],ym[1:nd],wm[1:nd]/2,hm[1:nd]/2,am[1:nd])
;qr_rectangles(0,0,hwid,hwid,0)

polyg1 = POLYGON([-HWID,-HWID,HWID,HWID],[-HWID,HWID,HWID,-HWID],/data,color='r')
FOR i = 0,nd-1 DO BEGIN
  
  asd=TEXT(xm(i),ym(i)+10,STRING(plab(i)),/data,alignment=0.5,font_size=10)
  polyg = POLYGON([xm(i)-wm(i)/2.,xm(i)-wm(i)/2.,xm(i)+wm(i)/2.,xm(i)+wm(i)/2.],$
    [ym(i)-hm(i)/2.,ym(i)+hm(i)/2.,ym(i)+hm(i)/2.,ym(i)-hm(i)/2.],/data,color='g')
  blurb= STRTRIM(rm(i),2)
  asdf=TEXT(xm(i),ym(i),blurb,/data,alignment=0.5,font_size=10)
  blurb= STRMID(lm(i),STRPOS(lm(i),'.')-1,4) + ' ' + STRMID(loverd(i),STRPOS(loverd(i),'.')-3,5)
  asdfg=TEXT(xm(i),ym(i)-10,blurb,/data,alignment=0.5,font_size=8)
ENDFOR



RETURN
END