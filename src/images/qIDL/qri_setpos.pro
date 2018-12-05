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
;       A = QRI_GETPOS(ipos,p)
;
;  IPOS = integer that determines the behaviour:
;           input  1 pixel 0-NCOLS, 0-NROWS
;                  2 local X,Y
;                  3 local azimuth,elevation degrees
;                  4 Celestial RA,DEC degrees J2000
;                  5 Ecliptic EA,EL degrees
;                  6 Galactic LII,BII degrees
;           output
;                   P = Doublearr(2) containing the position to set:
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

FUNCTION QRI_SETPOS,ipos,p



qsoft = GETENV('QSOFT')
A=CALL_EXTERNAL(qsoft + '/R_libraries/imagesR.so', "qri_setpos_", $
  LONG(ipos),$
  DOUBLE(p),$
/AUTO_GLUE)

RETURN, 'y'
END