;+
;
; NAME:
;   QRT_LENS
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
;     A = QRT_LENS(id,idf,iq,anml,arfx,apos,rap,r1,r2,refind,thick)
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

FUNCTION QRT_LENS,id,idf,iq,anml,arfx,apos,rap,r1,r2,refind,thick


  qsoft = GETENV('QSOFT')
  A=CALL_EXTERNAL(qsoft + '/R_libraries/xsrtR.so', "qrt_lens_", $
        LONG(id),$
        LONG(idf),$
        LONG(iq),$
        DOUBLE(anml),$
        DOUBLE(arfx),$
        DOUBLE(apos),$
        DOUBLE(rap),$
        DOUBLE(r1),$
        DOUBLE(r2),$
        DOUBLE(refind),$
        DOUBLE(thick),$
    /AUTO_GLUE)


  RETURN,A
END
