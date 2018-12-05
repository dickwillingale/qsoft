;+
;
; NAME:
;   QR_RECTANGLES
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
;       A = QR_RECTANGLES(x,y,w,h,t,color=color,fill=fill,wind=wind)
;   INPUTS:
;       x = array of centre of rectangles in x
;       y = array of centre of rectangles in y
;       w = array of widths of rectangles
;       h = array of heights of rectangle
;       t = theta rotations of rectangles
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


FUNCTION QR_RECTANGLES,x,y,w,h,t,color=color,fill=fill,wind=wind

IF ~KEYWORD_SET(color) THEN poly_col='black' ELSE poly_col=color
IF ~KEYWORD_SET(fill) THEN poly_fill=0 ELSE poly_fill=1


  cth=COS(t)
  sth=SIN(t)
  xa=w*cth
  xb=w*sth
  ya=h*cth
  yb=h*sth
  x1=x+xa+yb
  x2=x+xa-yb
  x3=x-xa-yb
  x4=x-xa+yb
  y1=y+xb-ya
  y2=y+xb+ya
  y3=y-xb+ya
  y4=y-xb-ya
  
  IF KEYWORD_SET(wind) THEN poly_wind=wind ELSE BEGIN
    newplot=plot([x1,x2,x3,x4]*1.1,[y1,y2,y3,y4]*1.1,linestyle=6,aspect_ratio=1)
    poly_wind=newplot
  ENDELSE
  
  FOR i = 0,N_ELEMENTS(x)-1 DO BEGIN
    pol=POLYGON([x1(i),x2(i),x3(i),x4(i)],[y1(i),y2(i),y3(i),y4(i)],color=poly_col,FILL_BACKGROUND=poly_fill,TARGET=poly_wind,/data)
  ENDFOR

RETURN, pol
END