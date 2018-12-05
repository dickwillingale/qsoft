PRO plot_mixst005_mcp

  infile = '/Users/Adrian/Documents/MIXS/Modelling and Data Analysis/Calibration documents/momt005 - FM/Al/imgt04279.raw'
;infile = '/Users/Adrian/Documents/MIXS/Modelling and Data Analysis/Calibration documents/sext020 - FM/images/Al/imgt03549.raw'
linfile= '/Users/Adrian/Documents/MIXS/Modelling and Data Analysis/Calibration documents/imgt03275_new.lin'
;infile = '/Users/Adrian/Desktop/images/se001-a14/images/lob/imgt01848.raw'
  aa= QR_RAWREAD(infile,4.,/linear,linfile)
  
  bb= STRSPLIT(STRING(aa.header),/extract)
  cc= STRSPLIT(bb(4),':',/extract)
 ; filestart = FLOAT(cc(0))*60.*60.+FLOAT(cc(1))*60.+FLOAT(cc(2))
 ; stopt = ''
  ;print,bb
  

  hwid = 35.
  np = 512.
  pix = (hwid*2./np)/10. ; pixel size in cm
  
  id=QRI_BINXY(aa.xp,aa.yp,0,1,-hwid,hwid,np,-hwid,hwid,np)
  ;id=QRI_BINXY(aa.xp,aa.yp,0,1,-20,20,np,0,40,np)
  WINDOW,1,xsize=2^9.5,ysize=2^9.5
  cgLOADCT,1,/reverse
  
 id.xlab = 'x (mm)'
 id.ylab = 'y (mm)'
  QRI_ZOOMIMAGE,id,[-hwid,hwid],[-hwid,hwid],[0,MAX(id.data_array)]
 ; QRI_ZOOMIMAGE,id,[-20,20],[-0,40],[0,MAX(id.data_array)]

  IF 0 THEN BEGIN
  
  PRINT, 'SELECT CENTRE OF PATCH'
  CURSOR, xcen, ycen, /down,/data
  print, xcen,ycen
  
  xmin = (xcen - 19.2/2.+hwid)/(pix*10.)
  xmax = (xcen + 19.2/2.+hwid)/(pix*10.)
  ymin = (ycen - 19.2/2.+hwid)/(pix*10.)
  ymax = (ycen + 19.2/2.+hwid)/(pix*10.)
  
  plots,[xcen - 19.2/2.,xcen - 19.2/2.,xcen + 19.2/2.,xcen + 19.2/2.,xcen - 19.2/2.], $
    [ycen - 19.2/2.,ycen + 19.2/2.,ycen + 19.2/2.,ycen - 19.2/2.,ycen - 19.2/2.],color=cgcolor('red')
  
  AefT = TOTAL(id.data_array(xmin:xmax,ymin:ymax))
  PRINT, 'MIXS Effective area = ', AefT
  ENDIF
  
  cgzplot, findgen(np)*2.*hwid/np-hwid, TOTAL(id.data_array,2),psym=10,xtitle='X (mm)',ytitle = 'Intensity (arb)'
  cgzplot, findgen(np)*2.*hwid/np-hwid, TOTAL(id.data_array,1),psym=10,xtitle='X (mm)',ytitle = 'Intensity (arb)'


RETURN
END