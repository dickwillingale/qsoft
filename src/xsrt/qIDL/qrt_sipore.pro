;+
;
; NAME:
;   QRT_SIPORE
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
;     A = QRT_SIPORE(pcen,pnorm,raxis,flen,rpitch,apitch,wall,rrm,ppm,tm,wm,hm,am,cm,gm,wfr,a2j,idf,iq)
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

FUNCTION QRT_SIPORE,pcen,pnorm,raxis,flen,rpitch,apitch,wall,rrm,ppm,tm,wm,hm,am,cm,gm,wfr,a2j,idf,iq


  qsoft = GETENV('QSOFT')
  A=CALL_EXTERNAL(qsoft + '/R_libraries/xsrtR.so', "qrt_sipore_", $
        DOUBLE(pcen),$
        DOUBLE(pnorm),$
        DOUBLE(raxis),$
        DOUBLE(flen),$
        DOUBLE(rpitch),$
        DOUBLE(apitch),$
        DOUBLE(wall),$
        LONG(N_ELEMENTS(rm)),$
        DOUBLE(rm),$
        DOUBLE(pm),$
        DOUBLE(tm),$
        DOUBLE(wm),$
        DOUBLE(hm),$
        DOUBLE(am),$
        DOUBLE(cm),$
        DOUBLE(gm),$
        DOUBLE(wfr),$
        DOUBLE(a2j),$
        LONG(idf),$
        LONG(iq),$
    /AUTO_GLUE)


  RETURN,A
END
