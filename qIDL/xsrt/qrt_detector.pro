FUNCTION QRT_DETECTOR,id,dpos,dnml,drfx,dlim,radet

;  *ID        input        detector type
;  *DPOS      input        detector position
;  *DNML      input        detector normal
;  *DRFX      input        detector reference axis
;  *DLIM      input        detector limits
;  *RADET     input        radius of curvature of spherical detector


qsoft = GETENV('QSOFT')
  A=CALL_EXTERNAL(qsoft + '/R_libraries/xsrtR.so', "qrt_detector_", $
        LONG(id),$
        DOUBLE(dpos),$
        DOUBLE(dnml),$
        DOUBLE(drfx),$
        DOUBLE(dlim),$
        DOUBLE(radet),$
    /AUTO_GLUE)
    


RETURN, a
END