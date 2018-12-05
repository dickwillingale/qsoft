FUNCTION QR_COLLUT, reverse=reverse

;read QR colour table
qsoft = GETENV('QSOFT')
infile = qsoft + '/data/lut6.dat'


  OPENR, lun, infile, /Get_Lun
  header = ''
  nl = FILE_LINES(infile)-1
  data = {i:FIX(0),r:0.0,g:0.0,b:0.0}
  t = REPLICATE(data,nl)
  READF, lun, header, t
  FREE_LUN, lun


IF KEYWORD_SET(reverse) THEN BEGIN
  t.r=TEMPORARY(REVERSE(BYTE(BYTSCL(t.r))))
  t.g=TEMPORARY(REVERSE(BYTE(BYTSCL(t.g))))
  t.b=TEMPORARY(REVERSE(BYTE(BYTSCL(t.b))))
ENDIF ELSE BEGIN
  t.r=TEMPORARY(BYTE(BYTSCL(t.r)))
  t.g=TEMPORARY(BYTE(BYTSCL(t.g)))
  t.b=TEMPORARY(BYTE(BYTSCL(t.b)))
ENDELSE


  TVLCT, t.r, t.g, t.b
  tab = FLTARR(3,256)
  tab(0,*)=t.r
  tab(1,*)=t.g
  tab(2,*)=t.b

RETURN,tab
END