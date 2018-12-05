;+
;
; NAME:
;   QRI_DISPLAYRAYS
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
;       A = QRI_DISPLAYRAYS()
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


FUNCTION QRI_DISPLAYRAYS,a,nmax,az,el,rl,lims=lims

  aer = [az,el,rl]*!dpi/180.
  b = QR_XYZTRANS(a.RXP,a.RYP,a.RZP,[0,0,0],aer)
  
  IF ~KEYWORD_SET(lims) THEN BEGIN
    xmin = MIN(b.xp)
    xmax = MAX(b.xp)
    ymin = MIN(b.yp)
    ymax = MAX(b.yp)
  ENDIF ELSE BEGIN
    xmin = lims(0)
    xmax = lims(1)
    ymin = lims(2)
    ymax = lims(3)
  ENDELSE


ray = plot([xmin,xmax],[ymin,ymax], yrange= [ymin,ymax]*1.2, xrange = [xmin,xmax]*1.1, xtitle="x mm",ytitle="y mm", linestyle=6, aspect_ratio=1, name = 'rays')

np = N_ELEMENTS(a.IQU)
nplot = 0
is = 0
i = 0

colors = ['black','green','red','blue','yellow','cyan','magenta','purple']

WHILE (nplot LT nmax) AND (i LE np) DO BEGIN
  ii = i+1
  IF (a.IQU(ii) EQ -2) OR (i EQ np) THEN BEGIN
      IF i-is-1 LE N_ELEMENTS(colors)-1 THEN raycol = colors(i-is-1) ELSE raycol = colors(N_ELEMENTS(colors)-1)
      lns = POLYLINE(b.xp(is:i),b.yp(is:i), color=raycol, target = 'rays', /DATA)
      nplot = nplot+1
      is = ii
  ENDIF
  i = ii
ENDWHILE

FOR j=0,N_ELEMENTS(colors)-1 DO BEGIN
  IF j EQ 0 THEN a=TEXT(XMIN+(XMAX-XMIN)*(j+2)/N_ELEMENTS(colors)/2,YMAX,'N HITS = ', color='black',/data,font_size=20,target='rays')
  b=TEXT(XMIN+(XMAX-XMIN)*(j+5)/N_ELEMENTS(colors)/2,YMAX,STRTRIM(j,2), color=colors(j),/data,font_size=20,target='rays')
ENDFOR


RETURN, 'y'
END