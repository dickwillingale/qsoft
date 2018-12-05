;+
;
; NAME:
;   QRT_W1NEST
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
;     A = QRT_W1NEST
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

FUNCTION QRT_W1NEST,xj,rj,ra,pl,ph,hl,hh,tin,tj,tout,ax,ar,ff,idf,iq,ib


  qsoft = GETENV('QSOFT')
  A=CALL_EXTERNAL(qsoft + '/R_libraries/xsrtR.so', "qrt_w1nest_", $
      LONG(N_ELEMENTS(rj)),$
      DOUBLE(xj),$
      DOUBLE(rj),$
      DOUBLE(ra),$
      DOUBLE(pl),$
      DOUBLE(ph),$
      DOUBLE(hl),$
      DOUBLE(hh),$
      DOUBLE(tin),$
      DOUBLE(tj),$
      DOUBLE(tout),$
      DOUBLE(ax),$
      DOUBLE(ar),$
      DOUBLE(ff),$
      LONG(idf),$
      LONG(iq),$
      LONG(ib),$
    /AUTO_GLUE)


  RETURN,A
END
