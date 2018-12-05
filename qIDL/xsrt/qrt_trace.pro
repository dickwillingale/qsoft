FUNCTION QRT_TRACE,ideb,riris,iadju

;  *IDEB        input        debugging level (0 none)
;  *RIRIS       input        radius about centre of detector for analysis
;  *                if 0.0 then no analysis of detected distribution
;  *IOPT        input  <0 save traced.dat and detected.dat files
;  *                   >0 adjust focus and save traced.dat and detected.dat
;  *                   only rays with IOPT reflections are used in adjustment
;  *                   =0 don't save or adjust focus
;  *AREA        output        detected area within RIRIS
;  *DSHFT       output        axial shift to optimum focus (0.0 if IOPT<=0)
;  *YBAR        output        y centroid of detected distribution
;  *ZBAR        output        z centroid of detected distribution
;  *RMS         output        rms radius of detected distribution


  area=0d
  dshft=0d
  ybar=0d
  zbar=0d
  rms=0d

qsoft = GETENV('QSOFT')
  A=CALL_EXTERNAL(qsoft + '/R_libraries/xsrtR.so', "qrt_trace_", $
        LONG(ideb), $
        DOUBLE(riris), $
        LONG(iadju), $
        DOUBLE(area), $
        DOUBLE(dshft), $
        DOUBLE(ybar), $
        DOUBLE(zbar), $
        DOUBLE(rms), $
    /AUTO_GLUE)
    
    rtn = {area:area,dshft:dshft,ybar:ybar,zbar:zbar,rms:rms}

RETURN, rtn
END