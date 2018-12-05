;+
;
; NAME:
;   QRI_GETPOS
;
; PURPOSE:
;   This Library implememnts th QSOFT analysis package in IDL using calls to sharable objects
;   external to IDL. Original routines were in general written in Fortan 77 and are compiled
;   into sharable objects.
;
;   The purpose of this routine is to find QSOFT's current position
;
;
; :Categories:
;    Image analysis
;
; :Keywords:
;
; :Examples:
;       A = QRI_GETPOS()
;     OUTPUTS: structure containing position in following coordianates:
;                  1 pixel 0-NCOLS, 0-NROWS
;                  2 local X,Y
;                  3 local azimuth,elevation degrees
;                  4 Celestial RA,DEC degrees J2000
;                  5 Ecliptic EA,EL degrees
;                  6 Galactic LII,BII degrees
;
; :Author:
;       Adrian Martindale
;
; :History:
;     Change History::
;       Written by: Adrian Martindale, 15/10/2018 to implement QSOFT in IDL.
;
;-
;

FUNCTION QRI_GETPOS

PIX=[0.d,0.d]
XYL=[0.d,0.d]
AES=[0.d,0.d]
EQU=[0.d,0.d]
ECL=[0.d,0.d]
GAL=[0.d,0.d]

qsoft = GETENV('QSOFT')
A=CALL_EXTERNAL(qsoft + '/R_libraries/imagesR.so', "qri_getpos_", $
  DOUBLE(PIX),$
  DOUBLE(XYL),$
  DOUBLE(AES),$
  DOUBLE(EQU),$
  DOUBLE(ECL),$
  DOUBLE(GAL),$
/AUTO_GLUE)

out = {PIX:pix,XYL:xyl,AES:aes,EQU:equ,ECL:ecl,GAL:ecl}

RETURN, out
END