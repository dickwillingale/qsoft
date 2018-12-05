;+
;
; NAME:
;   QRT_WOLTER2
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
;     A = QRT_WOLTER2(rp,gp,rh,gh,rm,fovr,ax,ar,ff,idf,iq)
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

FUNCTION QRT_WOLTER2,rp,gp,rh,gh,rrm,fovr,ax,ar,ff,idf,iq


  qsoft = GETENV('QSOFT')
  A=CALL_EXTERNAL(qsoft + '/R_libraries/xsrtR.so', "qrt_wolter2_", $
        DOUBLE(rp),$
        DOUBLE(gp),$
        DOUBLE(rh),$
        DOUBLE(gh),$
        DOUBLE(rrm),$
        DOUBLE(fovr),$
        DOUBLE(ax),$
        DOUBLE(ar),$
        DOUBLE(ff),$
        LONG(idf),$
        LONG(iq),$
    /AUTO_GLUE)


  RETURN,A
END
