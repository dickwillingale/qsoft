PRO plot_MCP_data

infile = '/Users/Adrian/Documents/MIXS/Modelling and Data Analysis/Calibration documents/momc006 - FM/images/Al/collimating mode/imgt03436.raw'
linfile= '/Users/Adrian/Documents/MIXS/Modelling and Data Analysis/Calibration documents/imgt03275_new.lin'

infie = '/Users/Adrian/Desktop/images/imgt00012.raw'

  aa= QR_RAWREAD(infile,4.,/linear,linfile)
  
  bb= STRSPLIT(STRING(aa.header),/extract)
  cc= STRSPLIT(bb(4),':',/extract)
  filestart = FLOAT(cc(0))*60.*60.+FLOAT(cc(1))*60.+FLOAT(cc(2))
  stopt = ''
  print,bb
  READ,stopt,prompt = 'Supply stoptime'
  stoptime = STRSPLIT(stopt,':',/extract)
  filestop = FLOAT(stoptime(0))*60.*60. + FLOAT(stoptime(1))*60. + FLOAT(stoptime(2))
  time =filestop-filestart

  hwid = 15.
  np = 512.
  pix = (hwid*2./np)/10. ; pixel size in cm
  
  id=QRI_BINXY(aa.xp,aa.yp,0,1,-hwid,hwid,np,-hwid,hwid,np)
  WINDOW,1,xsize=768,ysize=768
  cgLOADCT,1,/reverse
  
  inf=''
  READ,inf, prompt = 'supply PN diode information for calibration (epsilon, tPN, nPN)'
  split_inf = STRSPLIT(inf, ',', /extract)
  epsilon = FLOAT(split_inf(0))
  tPN = FLOAT(split_inf(1))
  nPN = FLOAT(split_inf(2))
  
 ; convert to effective area per pixel using PN diode data and known noise in MCP detecotr ~30cps over 9x9cm
 noise = 0;0.37;counts/cm^2/s
 noisePerPixel = noise*time*pix^2
  id.data_array = (id.data_array - noiseperpixel)/(time) * epsilon * (tPN * 0.080425)/nPN
  
  id.zlim(1) = MAX(id.data_array)
  QRI_ZOOMIMAGE,id,[-hwid,hwid],[-hwid,hwid],[0,1.5E-4]
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
  

RETURN
END