  ;+
  ;
  ; NAME:
  ;   QRI_BEAM
  ;
  ; PURPOSE:
  ;   This Library implememnts th QSOFT analysis package in IDL using calls to sharable objects
  ;   external to IDL. Original routines were in general written in Fortan 77 and are compiled 
  ;   into sharable objects.
  ;    
  ;   The purpose of this routine is to analyse a beam or point source in an image. 
  ;
  ;
  ; :Categories:
  ;    Image analysis
  ;
  ; :Keywords:
  ;    raw_arr: in, optional, type=dblarr(2)
  ;       Set this keyword to tell QSOFT that you are analysing a standard array 
  ;       (not one generated in QSOFT), such that other variables in the fortran 
  ;       common blocks are initialised correctly.
  ;
  ; :Examples:
  ;   A=QRI_BEAM(arr, rbeam, blev, bvar, RAW_ARR = [X0,Y0])
  ;     INPUTS: 
  ;       arr = image to be analysed should be double precision 
  ;       rbeam = radius of beam to analyse
  ;       blev = ??
  ;       bvar = ??
  ;
  ; :Author:
  ;       Adrian Martindale
  ;
  ; :History:
  ;     Change History::
  ;       Written by: Adrian Martindale, 15/10/2018 to implement QSOFT in IDL.
  ;     
  ;-
  
FUNCTION QRI_BEAM, arr,rbeam,blev,bvar, raw_arr = raw

IF KEYWORD_SET(raw) THEN BEGIN
  init = QRI_INIT()
  set_pos = QRI_SETPOS(2,[raw(0)*!DTOR,raw(1)*!DTOR])
ENDIF

  nsam=0l
  bflux=0.d
  
  
  bsigma=0.d
  flux=0.d
  fsigma=0.d
  peak=[0.d,0.d]
  cen=[0.d,0.d]
  tha=0.d
  rmsa=0.d
  rmsb=0.d
  fwhmm=0.d
  hew=0.d
  w90=0.d
  fwhmp=0.d
  hewp=0.d
  w90p=0.d
  fwhmc=0.d
  hewc=0.d
  w90c=0.d
qsoft = GETENV('QSOFT')
A=CALL_EXTERNAL(qsoft + '/R_libraries/imagesR.so', "qri_beam_", $
  LONG(N_ELEMENTS(arr(0,*))),$
  LONG(N_ELEMENTS(arr(*,0))),$
  DOUBLE(arr),$
  DOUBLE(rbeam),$
  DOUBLE(blev),$
  DOUBLE(bvar),$
  nsam,$
  bflux,$
  bsigma,$
  flux,$
  fsigma,$
  peak,$
  cen,$
  tha,$
  rmsa,$
  rmsb,$
  fwhmm,$
  hew,$
  w90,$
  fwhmp,$
  hewp,$
  w90p,$
  fwhmc,$
  hewc,$
  w90c,$
/AUTO_GLUE)

; Do peak fit
IF bvar NE 0 THEN BEGIN
  delstat = qchisq(0.9,4)
  
  pval = arr(a.cen(0),a.cen(1))
  spars = [pval,a.cen(0),a.cen(1),a.fwhmc/2.]
  lpars = [pval/2,a.cen(0)-a.fwhmc/2,a.cen(1)-a.fwhmc/2,a.fwhmc/2./2]
  upars = [pval*2,a.cen(0)+a.fwhmc/2,a.cen(1)+a.fwhmc/2,a.fwhmc/2.*2]
  derr = ['F','T','T','T']
  f = qr_srchmin(spars,lpars,upars,'qri_peakchisq',delstat,derr)
ENDIF ELSE f = 'F'

;# The fit parameters are saved in the list fit returned
;# 1 peak value (no error range calculated)
;# 2 peak X pixel position including 90% upper and lower bounds
;# 3 peak Y pixel position including 90% upper and lower bounds
;# 4 Gaussian sigma including 90% upper and lower bounds
;# add parameters to structure of results

beam = {nsam: nsam,$
bflux:bflux,$
bsigma:bsigma,$
flux:flux,$
fsigma:fsigma,$
peak:peak,$
cen:cen,$
tha:tha,$
rmsa:rmsa,$
rmsb:rmsb,$
fwhm:fwhmm,$
hew:hew,$
w90:w90,$
fwhmp:fwhmp,$
hewp:hewp,$
w90p:w90p,$
fwhmc:fwhmc,$
hewc:hewc,$
w90c:w90c, $
fit:f $
}

RETURN, beam
END