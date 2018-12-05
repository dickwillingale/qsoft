;+
;
; NAME:
;   QRI_SETSKY
;
; PURPOSE:
;   This Library implememnts th QSOFT analysis package in IDL using calls to sharable objects
;   external to IDL. Original routines were in general written in Fortan 77 and are compiled
;   into sharable objects.
;
;   The purpose of this routine is to insert the sky coordinates for an image field into the fortran common blocks
;
;
; :Categories:
;    Image analysis
;
; :Keywords:
;
; :Examples:
;     A = QRI_SETSKY(xtodeg,ytodeg,ipr,mjd,ra,dec,roll)
;       INPUTS:  
;         *XTODEG     input    scale from local X to degrees (usually -ve)
;         *YTODEG     input    scale from local Y to degrees
;         *IPR        input    projection between local XY and local spherical
;         *                0 Plate Carre
;         *                1 Aitoff (Hammer)
;         *                2 Lambert (equatorial aspect of Azimuthal equal-area
;         *                  projection)
;         *MJD        input    Modified Julian date
;         *R          input    Right Ascension (degrees J2000 at local origin)
;         *DEC        input    Declination (degrees J2000 at local origin)
;         *ROLL       input    Roll angle (degrees from North to +ve elev. +ve clockwise)

; :Author:
;       Adrian Martindale
;
; :History:
;     Change History::
;       Written by: Adrian Martindale, 15/10/2018 to implement QSOFT in IDL.
;
;-
;

FUNCTION QRI_SETSKY, xtodeg,ytodeg,ipr,mjd,ra,dec,roll

qsoft = GETENV('QSOFT')
A=CALL_EXTERNAL(qsoft + '/R_libraries/imagesR.so', "qri_setsky_", $
  DOUBLE(xtodeg),$
  DOUBLE(ytodeg),$
  LONG(ipr),$
  DOUBLE(mjd),$
  DOUBLE(ra),$
  DOUBLE(dec),$
  DOUBLE(roll),$
/AUTO_GLUE)

RETURN,A
END
