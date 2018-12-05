
;+
;
; NAME:
;   QR_RAWREAD
;
; PURPOSE:
;   This Library implememnts th QSOFT analysis package in IDL using calls to sharable objects
;   external to IDL. Original routines were in general written in Fortan 77 and are compiled
;   into sharable objects.
;
;   The purpose of this routine is to analyse a beam or point source in an image.
;
;
; :Categories:
;    Image analysis
;
; :Keywords:
;    
;
; :Examples:
;   A=QR_RAWREAD
;
; :Author:
;       Adrian Martindale
;
; :History:
;     Change History::
;       Written by: Adrian Martindale, 15/10/2018 to implement QSOFT in IDL.
;
;-



FUNCTION QR_RAWREAD,filename,scale,linear = linear,linfile


; Read .raw MCP event files
; Adrian Martindale's version 21st July 2015
bfile = filename

hsize=512
head = BYTARR(hsize)

file_inf = FILE_INFO(filename)
idat = INTARR(4,(file_inf.size - hsize)/(4.*2.))

OPENR, unit, filename, /GET_LUN
READU, unit, head
READU, unit, idat

idat=FLOAT(idat)
n2 = N_ELEMENTS(idat)/4

;Dick's algorithm
; x= idat[3,0:n2-1]/(idat[1,0:n2-1]+idat[3,0:n2-1])-0.5
; y= idat[2,0:n2-1]/(idat[2,0:n2-1]+idat[4,0:n2-1])-0.5
;Image_Display algorithm
x = REFORM(idat(1,*)/(idat(1,*)+idat(3,*))-0.5)
y = REFORM(idat(2,*)/(idat(0,*)+idat(2,*))-0.5)

; linearise the data?
IF KEYWORD_SET(linear) THEN BEGIN
      ;format x and y as done by IDL (no shift and times number of pixels)
      x = (x+0.5)*512.0
      y = (y+0.5)*512.0
      ;read linearisation data file
      RESTORE, linfile
      
      ;reset scale from unlinearised estimate to real number
      mmperpix = pitch/scale
      scal = mmperpix*512.0
      ;NB 512/2 is half the width of the plot in image_display
      PRINT, "Scale set to = ",scal, " = ", mmperpix, " mm per pixel for an image_display pixel"
    
    nx = FLTARR(n2)
    ny = FLTARR(n2)
    FOR k = 0, nterms DO BEGIN
      FOR j = 0, nterms DO BEGIN
        tmp = x^k * y^j
        nx += xfit[j,k] * tmp
        ny += yfit[j,k] * tmp
      ENDFOR
    ENDFOR
    x = (nx/512.0)-0.5
    y = (ny/512.0)-0.5

ENDIF ELSE BEGIN
    scal=scale
ENDELSE

cs= cos(!dpi/4.)
ss= sin(!dpi/4.)

; need x and y to be the same for R and IDL. NB becaus of row vs column major the x-axis needs to be flipped
; to maintain consistency with IDL definitions which are implicit in the linearisation file.
xp= (y*ss-x*cs)*scal
yp= (x*ss+y*cs)*scal

ph= TOTAL(idat,1)
out= {xp:xp,yp:yp,ph:ph,header:head}
RETURN,OUT
END
