;+
;
; NAME:
;   QRT_OPGRAT
;
; PURPOSE:
;   This Library implememnts th QSOFT analysis package in IDL using calls to sharable objects
;   external to IDL. Original routines were in general written in Fortan 77 and are compiled
;   into sharable objects.
;
;   The purpose of this routine is to fit a lorentzian to the peak of the PSF
;
;
; :Categories:
;    ray tracing
;
; :Keywords:
;
; :Examples:
;     A = QRT_OPGRAT
;
;
; :Author:
;       Adrian Martindale
;
; :History:
;     Change History::
;       Written by: Adrian Martindale, 07/02/2018 to implement QSOFT in IDL.
;
;-
;

FUNCTION QRT_OPGRAT,id,idf,iq,alpha,fp,zp,gp,al

  adir=DBLARR(3)
  rdir=DBLARR(3)
  hpos=DBLARR(3)
  graz=DBLARR(1)
  gammaa=DBLARR(1)
  dhub=DBLARR(1)
  
  qsoft = GETENV('QSOFT')
  A=CALL_EXTERNAL(qsoft + '/R_libraries/xsrtR.so', "qrt_opgrat_", $
      LONG(id),$
      LONG(idf),$
      LONG(iq),$
      DOUBLE(alpha),$
      DOUBLE(fp),$
      DOUBLE(zp),$
      DOUBLE(gp),$
      LONG(N_ELEMENTS(al)),$
      DOUBLE(al),$
      adir,$
      rdir,$
      hpos,$
      graz,$
      gammaa,$
      dhub, $
    /AUTO_GLUE)


  RETURN,A
END
