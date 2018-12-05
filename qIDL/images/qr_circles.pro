;+
;
; NAME:
;   QR_CIRCLES
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
;       color=color of rectangle? string e.g. 'black', 'green'
;       fill=fill the rectangle? 1=yes,0=no
;       wind=name of graphic to add rectanges too. string eg 'plotname'
; :Examples:
;       A = QR_CIRCLES(x,y,r,color=color,fill=fill,wind=wind)
;   INPUTS:
;       x = array of centre of circles in x
;       y = array of centre of circles in y
;       r = array of radii of circles
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


FUNCTION QR_CIRCLES,x,y,r,color=color,fill=fill,wind=wind

IF ~KEYWORD_SET(color) THEN poly_col='black' ELSE poly_col=color
IF ~KEYWORD_SET(fill) THEN poly_fill=0 ELSE poly_fill=1  
  
  IF KEYWORD_SET(wind) THEN poly_wind=wind ELSE BEGIN
    newplot=plot([(x+r),(x-r)]*1.1,[(y+r),(y-r)]*1.1,linestyle=6,aspect_ratio=1)
    poly_wind=newplot
  ENDELSE
  
  FOR i = 0,N_ELEMENTS(x)-1 DO BEGIN
    cir=ELLIPSE(x(i),y(i),major=r(i),/data,color=poly_col,FILL_BACKGROUND=poly_fill,TARGET=poly_wind)
  ENDFOR

RETURN, cir
END