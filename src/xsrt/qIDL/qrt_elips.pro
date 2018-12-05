;+
;
; NAME:
;   QRT_ELIPS
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
;     A = QRT_ELIPS(org,axs,cen,xmin,xmax,amin,amax,smb,rab,ide,iq)
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

FUNCTION QRT_ELIPS,org,axs,cen,xmin,xmax,amin,amax,smb,rab,ide,iq


  qsoft = GETENV('QSOFT')
  A=CALL_EXTERNAL(qsoft + '/R_libraries/xsrtR.so', "qrt_elips_", $
        DOUBLE(org),$
        DOUBLE(axs),$
        DOUBLE(cen),$
        DOUBLE(xmin),$
        DOUBLE(xmax),$
        DOUBLE(amin),$
        DOUBLE(amax),$
        DOUBLE(smb),$
        DOUBLE(rab),$
        LONG(ide),$
        LONG(iq),$
    /AUTO_GLUE)


  RETURN,A
END
