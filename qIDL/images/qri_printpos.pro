;+
;
; NAME:
;   QRI_PRINTPOS
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
;     A = QRI_PRINTPOS()
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

FUNCTION QRI_PRINTPOS

a=QRI_GETPOS()

PRINT, 'PIXEL       ', a.pix , FORMAT ='(A12,2F10.5)'
PRINT, 'XYLOC       ', a.xyl , FORMAT ='(A12,2F10.5)'
PRINT, 'SPLOC       ', a.aes , FORMAT ='(A12,2F10.5)'
PRINT, 'EQUATORIAL  ', a.equ , FORMAT ='(A12,2F10.5)'
PRINT, 'ECLIPTIC    ', a.ecl , FORMAT ='(A12,2F10.5)'
PRINT, 'GALACTIC    ', a.gal , FORMAT ='(A12,2F10.5)'
RETURN,1
END
