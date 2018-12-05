;+
;
; NAME:
;   QR_SECTORS
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
;       A = QR_SECTORS(r,y,w,h,t,color=color,fill=fill,wind=wind)
;   INPUTS:
;       r = array of centre of sectors in r
;       phi = array of centre of sectors in phi
;       w = array of widths of sectors in r
;       dphi = array of width of sectors in phi
;       t = theta rotations of sectors
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


FUNCTION QR_SECTORS,r,phi,w,dphi,t,color=color,fill=fill,wind=wind

IF ~KEYWORD_SET(color) THEN poly_col='black' ELSE poly_col=color
IF ~KEYWORD_SET(fill) THEN poly_fill=0 ELSE poly_fill=1

Rin = R-w/2.
Rout = R+w/2.
PhiLow = phi-dphi/2.
PhiHi = phi+dphi/2


  
  IF KEYWORD_SET(wind) THEN poly_wind=wind ELSE BEGIN
    newplot=plot([0,Rin*cos(philow),rin*cos(philow),Rin*cos(phihi),rin*cos(phihi),rout*cos(phihi),rout*cos(phihi),rout*cos(philow),rout*cos(philow)]*1.1,$
      [0,Rin*sin(philow),rin*sin(philow),Rin*sin(phihi),rin*sin(phihi),rout*sin(phihi),rout*sin(phihi),rout*sin(philow),rout*sin(philow)]*1.1,linestyle=6,aspect_ratio=1)
    poly_wind=newplot
  ENDELSE
  
  FOR i = 0,N_ELEMENTS(r)-1 DO BEGIN
    OutArcX = Rout(i)*COS(PhiLow(i)+FINDGEN(1000)*dphi(i)/1000.)
    InArcX = Rin(i)*COS(PhiLow(i)+FINDGEN(1000)*dphi(i)/1000.)
    OutArcY = Rout(i)*SIN(PhiLow(i)+FINDGEN(1000)*dphi(i)/1000.)
    InArcY = Rin(i)*SIN(PhiLow(i)+FINDGEN(1000)*dphi(i)/1000.)    
    polyg = POLYGON([Rin(i)*COS(PhiLow(i)), Rout(i)*COS(PhiLow(i)),OutArcX, Rout(i)*COS(PhiHi(i)), Rin(i)*COS(PhiHi(i)), REVERSE(InArcX)],$
      [Rin(i)*SIN(PhiLow(i)), Rout(i)*SIN(PhiLow(i)), OutArcY, Rout(i)*SIN(PhiHi(i)), Rin(i)*SIN(PhiHi(i)), REVERSE(InArcY)],$
      color=poly_col,FILL_BACKGROUND=poly_fill,TARGET=poly_wind,/data)
  ENDFOR

RETURN, polyg
END