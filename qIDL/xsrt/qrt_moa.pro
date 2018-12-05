;+
;
; NAME:
;   QRT_MOA
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
;     A = QRT_MOA(pcen,pno,rax,rcur,xyap,pitch,wall,plen,idf,iq)
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

FUNCTION QRT_MOA,pcen,pno,rax,rcur,xyap,pitch,wall,plen,idf,iq


  qsoft = GETENV('QSOFT')
  A=CALL_EXTERNAL(qsoft + '/R_libraries/xsrtR.so', "qrt_moa_", $
        DOUBLE(pcen), $
        DOUBLE(pno), $
        DOUBLE(rax), $
        DOUBLE(rcur), $
        DOUBLE(xyap), $
        DOUBLE(pitch), $
        DOUBLE(wall), $
        DOUBLE(plen), $
        LONG(idf), $
        LONG(iq), $
    /AUTO_GLUE)


  RETURN,A
END
