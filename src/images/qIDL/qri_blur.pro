;+
;
; NAME:
;   QRI_BLUR
;
; PURPOSE:
;   This Library implememnts th QSOFT analysis package in IDL using calls to sharable objects
;   external to IDL. Original routines were in general written in Fortan 77 and are compiled
;   into sharable objects.
;
;   The purpose of this routine is to print the positions held in the common blocks to the screen
;
;
; :Categories:
;    Image analysis
;
; :Keywords:
;
; :Examples:
;     A = QRI_BLUR()
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

FUNCTION QRI_BLUR,a,m


  nx=N_ELEMENTS(a(0,*))
  ny=N_ELEMENTS(a(*,0))
  mx=N_ELEMENTS(m(0,*))
  my=N_ELEMENTS(m(*,0))
  arr=DBLARR(nx*ny)
  
  qsoft = GETENV('QSOFT')
  A=CALL_EXTERNAL(qsoft + '/R_libraries/imagesR.so', "qri_blur_", $
        LONG(nx),$
        LONG(ny),$
        DOUBLE(a),$
        LONG(mx),$
        LONG(my),$
        DOUBLE(m),$
        DOUBLE(arr), $
    /AUTO_GLUE)

RETURN,arr
END
