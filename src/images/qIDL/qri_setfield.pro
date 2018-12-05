;+
;
; NAME:
;   QRI_SETFIELD
;
; PURPOSE:
;   This Library implememnts th QSOFT analysis package in IDL using calls to sharable objects
;   external to IDL. Original routines were in general written in Fortan 77 and are compiled
;   into sharable objects.
;
;   The purpose of this routine is to insert the field parameters for an image into the fortran common block
;
;
; :Categories:
;    Image analysis
;
; :Keywords:
;
; :Examples:
;     A = QRI_SETFIELD(nx, xleft, xright, ny, ybot, ytop)
;       inputs:
;         *NX        input        number of columns
;         *XLEFT     input        local x left edge
;         *XRIGHT    input        local y right edge
;         *NY        input        number of rows
;         *YBOT      input        local y bottom edge
;         *YTOP      input        local y top edge

; :Author:
;       Adrian Martindale
;
; :History:
;     Change History::
;       Written by: Adrian Martindale, 15/10/2018 to implement QSOFT in IDL.
;
;-
;

FUNCTION QRI_SETFIELD, nx, xleft, xright, ny, ybot, ytop

   nx = LONG(nx)
   xleft = DOUBLE(xleft)
   xright = DOUBLE(xright)
   ny = LONG(ny)
   ybot = DOUBLE(ybot)
   ytop = DOUBLE(ytop)

qsoft = GETENV('QSOFT')
A=CALL_EXTERNAL(qsoft + '/R_libraries/imagesR.so', "qri_setfield_", $
  nx,$
  xleft,$
  xright,$
  ny,$
  ybot,$
  ytop,$
/AUTO_GLUE)

xsam = (xright-xleft)/nx
ysam = (ytop-ybot)/ny

B = {nx:nx,xleft:xleft,xright:xright,xsam:xsam,ny:ny,ybot:ybot,ytop:ytop,ysam:ysam, $
        xp: DINDGEN(nx)*xsam/nx + xleft + xsam/2, $
        yp: DINDGEN(ny)*ysam/ny + ybot + ysam/2}

RETURN,B
END
