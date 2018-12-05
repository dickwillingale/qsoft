;+
;
; NAME:
;   QRT_PRISM
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
;     A = QRT_PRISM
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

FUNCTION QRT_PRISM,id,idf,iq,anml,arfx,apos,rap,d1,d2,refind,thick


  qsoft = GETENV('QSOFT')
  A=CALL_EXTERNAL(qsoft + '/R_libraries/xsrtR.so', "qrt_prism_", $
      LONG(id),$
      LONG(idf),$
      LONG(iq),$
      DOUBLE(anml),$
      DOUBLE(arfx),$
      DOUBLE(apos),$
      DOUBLE(rap),$
      DOUBLE(d1),$
      DOUBLE(d2),$
      DOUBLE(refind),$
      DOUBLE(thick),$
    /AUTO_GLUE)


  RETURN,A
END
