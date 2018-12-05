;+
;
; NAME:
;   QR_QUADERR
;
; PURPOSE:
;   This Library implememnts th QSOFT analysis package in IDL using calls to sharable objects
;   external to IDL. Original routines were in general written in Fortan 77 and are compiled
;   into sharable objects.
;
;   The purpose of this routine is to bin data into 2d images and optionally apply a weighting 
;   to the binning. 
;
;
; :Categories:
;    Image analysis
;
; :Keywords:
;
; :Examples:
;       A = QR_QUADERR(x0,y0,x1,y1,y)
;   INPUTS:
;  # Quadratic estimator for confidence limit
;  # x0,y0   parameter and statistic at minimum
;  # x1,y1   parameter and statistic near minimum
;  # y   required statistic value
;  # returns estimate of parameter corresponding to y
;       
; :Author:
;       Adrian Martindale
;
; :History:
;     Change History::
;       Written by: Adrian Martindale, 15/01/2018 to implement QSOFT in IDL.
;
;-
;


FUNCTION QR_QUADERR,x0,y0,x1,y1,y

  a=(y1-y0)/(x1-x0)^2
  b=-2*a*x0
  c=y0-a*x0^2-b*x0
  ba=b^2-4*a*(c-y)

  IF (ba GE 0) AND (a NE 0) THEN RETURN, sqrt(ba)/(2*a) ELSE RETURN, 0

END