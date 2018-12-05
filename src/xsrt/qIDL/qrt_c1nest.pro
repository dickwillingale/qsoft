;+
;
; NAME:
;   QRT_C1NEST
;
; PURPOSE:
;   This Library implememnts th QSOFT analysis package in IDL using calls to sharable objects
;   external to IDL. Original routines were in general written in Fortan 77 and are compiled
;   into sharable objects.
;
;
;
; :Categories:
;    ray tracing
;
; :Keywords:
;
; :Examples:
;     A = QRT_C1NEST(xj,pl,ph,hl,hh,rpl,rph,rhl,rhh,tin,tj,tout,ax,ar,ff,idf,iq,ib)
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

FUNCTION QRT_C1NEST,xj,pl,ph,hl,hh,rpl,rph,rhl,rhh,tin,tj,tout,ax,ar,ff,idf,iq,ib

nw=(N_ELEMENTS(rpl)-1)*2+12

  qsoft = GETENV('QSOFT')
  A=CALL_EXTERNAL(qsoft + '/R_libraries/xsrtR.so', "qrt_c1nest_", $
      DOUBLE(ax),$
      DOUBLE(ar),$
      DOUBLE(ff),$
      LONG(iq),$
      LONG(ib),$
      LONG(idf),$
      LONG(N_ELEMENTS(rpl)),$
      DOUBLE(xj),$
      DOUBLE(pl),$
      DOUBLE(ph),$
      DOUBLE(hl),$
      DOUBLE(hh),$
      DOUBLE(rpl),$
      DOUBLE(rph),$
      DOUBLE(rhl),$
      DOUBLE(rhh),$
      DOUBLE(tin),$
      DOUBLE(tj),$
      DOUBLE(tout),$
      LONG(nw),$
      DBLARR(nw),$
    /AUTO_GLUE)


  RETURN,A
END
