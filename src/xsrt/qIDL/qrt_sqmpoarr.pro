;+
;
; NAME:
;   QRT_SQMPOARR
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
;     A = QRT_SQMPOARR(fpars)
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

FUNCTION QRT_SQMPOARR,pcen,pnorm,raxis,rcur,hwid,idf,ar


  qsoft = GETENV('QSOFT')
  A=CALL_EXTERNAL(qsoft + '/R_libraries/xsrtR.so', "qrt_sqmpoarr_", $
      DOUBLE(pcen), $
      DOUBLE(pnorm), $
      DOUBLE(raxis), $
      DOUBLE(rcur), $
      DOUBLE(hwid), $
      LONG(idf), $
      LONG(N_ELEMENTS(ar)), $
      DOUBLE(ar), $
    /AUTO_GLUE)





RETURN,A
END
