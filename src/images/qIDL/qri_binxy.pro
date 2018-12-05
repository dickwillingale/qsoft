;+
;
; NAME:
;   QRI_BINXY
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
;       A = QRI_BINXY(x,y,iq,w,xleft,xright,nx,ybot,ytop,ny)
;   INPUTS:
;       x = array of x-values (double)
;       y = array of y-values (double)
;       iq = surface qualities to accept
;       w = weights (either scalar or same size as x and y) 
;       xleft = x-axis low linit
;       xright = x-axis high limit
;       nx = number of bins in x
;       ybot = yaxis low limit
;       ytop = y-axis high limit
;       ny = number of bins in y
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


FUNCTION QRI_BINXY,x,y,iq,w,xleft,xright,nx,ybot,ytop,ny

data_array = DBLARR(nx,ny)

qsoft = GETENV('QSOFT')
A= CALL_EXTERNAL(qsoft + '/R_libraries/imagesR.so', 'qri_binxy_', $
  LONG(N_ELEMENTS(x))       ,$ ;input        number of events
  DOUBLE(X)                   ,$ ;input        array of x positions
  DOUBLE(Y)                   ,$ ;input        array of y positions
  LONG(N_ELEMENTS(IQ))      ,$ ;input        number of quality values (1 or N)
  LONG(IQ)                  ,$ ;input        array of quality values (0 for OK)
  LONG(N_ELEMENTS(W))       ,$ ;input        number of weights (1 or N)
  DOUBLE(W)                   ,$ ;input        array of weights
  DOUBLE(XLEFT)               ,$ ;input        minimum x value
  DOUBLE(XRIGHT)              ,$ ;input        maximum x value
  DOUBLE(YBOT)                ,$ ;input        minimum y value
  DOUBLE(YTOP)                ,$ ;input        maximum y value
  LONG(NX)                  ,$ ;input        dimensions of output array
  LONG(NY)                  ,$ ;input        dimensions of output array
  DOUBLE(DATA_ARRAY)          ,$
/AUTO_GLUE)

B = QRI_MAKEIMAGE(data_array,DOUBLE(xleft),DOUBLE(xright),DOUBLE(ybot),DOUBLE(ytop))

RETURN, B
END