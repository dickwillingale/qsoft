FUNCTION READ_QRT,infile, det=det, traced=traced

  IF KEYWORD_SET(det) THEN BEGIN
    ;READ the detected file
    OPENR, lun, infile, /Get_Lun
    nl = FILE_LINES(infile)-1
    header = ''
    dataStr = {XD: 0.0,YD:0.0, ZD:0.0,XC:0.0,YC:0.0,ZC:0.0,XR:0.0,YR:0.0,ZR:0.0,YDET:0.0,ZDET:0.0,AREA:0.0,NHIT1:0}
    srt_in = REPLICATE(dataStr,nl)
    READF, lun, header, srt_in
    FREE_LUN, lun
    
    srtdata = CREATE_STRUCT('xd',srt_in.xd, 'yd',srt_in.yd, 'zd',srt_in.zd, $
      'xc',srt_in.xc, 'yc',srt_in.yc, 'zc',srt_in.zc, $
      'xr',srt_in.xr, 'yr',srt_in.yr, 'zr',srt_in.zr, $
      'ydet',srt_in.ydet,'zdet',srt_in.zdet, $
      'area',srt_in.area, 'nhit1',srt_in.nhit1)
    
  ENDIF

  IF KEYWORD_SET(traced) THEN BEGIN

  traced_file = STRJOIN( StrSplit(infile, 'detected', /Regex, /Extract, /Preserve_Null), 'traced')
  ;READ the traced file
  OPENR, lun, traced_file, /Get_Lun
  nl = FILE_LINES(traced_file)-1
  header = ''
  dataStr = {rxp:0.0,ryp:0.0,rzp:0.0,area:0.0,iqu:0}
  srt_traced_in = REPLICATE(dataStr,nl)
  READF, lun, header, srt_traced_in
  FREE_LUN, lun
  
      IF ISA(srtdata) THEN BEGIN
        srtdata = CREATE_STRUCT(srtdata, 'rxp',srt_traced_in.rxp, 'ryp',srt_traced_in.ryp,'rzp',srt_traced_in.rzp, $
            'rarea',srt_traced_in.area, 'iqu',srt_traced_in.iqu)
      ENDIF ELSE BEGIN
        srtdata = CREATE_STRUCT('rxp',srt_traced_in.rxp, 'ryp',srt_traced_in.ryp,'rzp',srt_traced_in.rzp, $
            'rarea',srt_traced_in.area, 'iqu',srt_traced_in.iqu)
      ENDELSE
  ENDIF

;  IF ISA(srtdata) THEN srt_file = CREATE_STRUCT('srtdata',srtdata) ELSE RETURN, !values.F_nan
  IF ISA(srtdata) THEN srt_file = srtdata ELSE RETURN, !values.F_nan
RETURN, srt_file
END