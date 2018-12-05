;+
;
; NAME:
;   QR_XYZTRANS
;
; PURPOSE:
;   This Library implememnts th QSOFT analysis package in IDL using calls to sharable objects
;   external to IDL. Original routines were in general written in Fortan 77 and are compiled
;   into sharable objects.
;
;   The purpose of this routine is to perform a 3-axis transformation given the x,y,z coordinates 
;   of the current position, a new origin and the azimuth, elevation and roll of a new coordinate 
;   system
;
;
; :Categories:
;    coordinate transforms
;
; :Keywords:
;
; :Examples:
;     A = QR_XYZTRANS(x, y, z, cen, aer)
;       inputs:
;           x,y,z = double array of lenght n listing x, y and z positions in current coordinate system
;           cen = centre of new coordiante system
;           aer = azimuth, elevation and roll of new x-axis in radians roll angle from Z-axis to +ve elev. +ve clockwise
;       outputs:
;          structure containing tags xp,yp and zp which are the projected positions in the new coordinate system.
; :Author:
;       Adrian Martindale
;
; :History:
;     Change History::
;       Written by: Adrian Martindale, 15/10/2018 to implement QSOFT in IDL.
;
;-
;

FUNCTION QR_XYZTRANS, x, y, z, cen, aer

len = N_ELEMENTS(x)
  xp=DBLARR(len)
  yp=DBLARR(len)
  zp=DBLARR(len)

qsoft = GETENV('QSOFT')
A=CALL_EXTERNAL(qsoft + '/R_libraries/imagesR.so', "qr_xyztrans_", $
  LONG(N_ELEMENTS(x)),$
  DOUBLE(x),$
  DOUBLE(y),$
  DOUBLE(z),$
  DOUBLE(cen),$
  DOUBLE(aer),$
  DOUBLE(xp),$
  DOUBLE(yp),$
  DOUBLE(zp),$
/AUTO_GLUE)

out = {xp:xp,yp:yp,zp:zp}

RETURN,OUT
END
