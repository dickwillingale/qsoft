;+
;
; NAME:
;   QRI_DISPAYIMAGE
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
;     A = QRI_DISPLAYIMAGE()
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

FUNCTION QRI_DISPLAYIMAGE,ima
qri_zoomimage,ima.data_array,ima.xlim,ima.ylim,ima.zlim
RETURN,'y'
END
