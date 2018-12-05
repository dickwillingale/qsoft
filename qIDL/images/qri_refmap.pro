;+
;
; NAME:
;   QRI_REFMAP
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
;       A = QRI_REFMAP()
;   INPUTS:
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

FUNCTION QRI_REFMAP,aa,nbinsy,nbinsz,lims=lims

  IF ~KEYWORD_SET(lims) THEN BEGIN
    ymin = MIN(aa.ryp)
    ymax = MAX(aa.ryp)
    zmin = MIN(aa.rzp)
    zmax = MAX(aa.rzp)
  ENDIF ELSE BEGIN
    ymin = lims(0)
    ymax = lims(1)
    zmin = lims(2)
    zmax = lims(3)
  ENDELSE

  

  ; filter by interception type in order to find detected rays etc
  optic_intercept_indices = WHERE(aa.IQU EQ -2)
  Ray_stops = SHIFT(optic_intercept_indices - 1,-1);WHERE((aa.IQU EQ -1) OR (aa.IQU EQ 0))
  ;set last elements to a real number rather than -1 to make sure nhit1 is correct
  Ray_stops(-1) = N_ELEMENTS(aa.IQU)
  detected_indices = WHERE(aa.IQU EQ -1)
  last_ref_indices = WHERE(aa.IQU EQ -1) - 1; NOTE -1 after the "where" to find the last reflection rather than detection

  ; calculate number of hits per ray and store in array of the same length a iqu
  nhit1 = FLTARR(N_ELEMENTS(aa.iqu))
  nhit1(ray_stops) = DOUBLE(Ray_stops) - DOUBLE(optic_intercept_indices) - 1.d
  

  ; create the images
  ; @last reflection point
  Last_Ref_area = QRI_BINXY(aa.ryp(last_ref_indices),aa.rzp(last_ref_indices),0,aa.rarea(last_ref_indices)/100.,ymin,ymax,nbinsy,zmin,zmax,nbinsz)
  optrays = QRI_BINXY(aa.ryp(last_ref_indices),aa.rzp(last_ref_indices),0,1.,ymin,ymax,nbinsy,zmin,zmax,nbinsz)
  optrays_nhit = QRI_BINXY(aa.ryp(last_ref_indices),aa.rzp(last_ref_indices),0,nhit1(detected_indices),ymin,ymax,nbinsy,zmin,zmax,nbinsz) ; NB USING DETECTED_INDICES FOR NHIT! IS DELIBERATE!!!
  optrays_avg = DOUBLE(optrays_nhit.data_array)/DOUBLE(optrays.data_array)

  ;@optic aperture
  Optic_aperture = QRI_BINXY(aa.ryp(optic_intercept_indices),aa.rzp(optic_intercept_indices),0,aa.rarea(optic_intercept_indices)/100., $
    ymin,ymax,nbinsy,zmin,zmax,nbinsz)


  ;@detector
  ; effective area image detector plane
  detimg = QRI_BINXY(aa.ryp(detected_indices),aa.rzp(detected_indices),0,aa.rarea(detected_indices)/100.,ymin,ymax,nbinsy,zmin,zmax,nbinsz)
  detrays = QRI_BINXY(aa.ryp(detected_indices),aa.rzp(detected_indices),0,1.,ymin,ymax,nbinsy,zmin,zmax,nbinsz)
  detrays_nhit = QRI_BINXY(aa.ryp(detected_indices),aa.rzp(detected_indices),0,nhit1(detected_indices),ymin,ymax,nbinsy,zmin,zmax,nbinsz)
  detrays_avg = DOUBLE(detrays_nhit.data_array)/DOUBLE(detrays.data_array)
  
  
  ;load ray hits color table
  colourT = RAY_HITS_CT()
  window,11,xsize=2^9.5,ysize=2^9.5
  mult=!P.multi
  !P.multi=[0,2,2]
  QRI_ZOOMIMAGE,detrays_avg,[ymin,ymax],[zmin,zmax],[0,7]
    clut=qr_collut(/REVERSE)
  QRI_ZOOMIMAGE,detimg.data_array,[ymin,ymax],[zmin,zmax],[0,MAX(detimg.data_array)/(10000.)]
    colourT = RAY_HITS_CT()
  QRI_ZOOMIMAGE,optrays_avg,[ymin,ymax],[zmin,zmax],[0,7]
     clut=qr_collut(/REVERSE)
  QRI_ZOOMIMAGE,optic_aperture.data_array,[ymin,ymax],[zmin,zmax],[MIN(optic_aperture.data_array),MAX(optic_aperture.data_array)] 

  !P.multi=mult

  RETURN, 'y'
END