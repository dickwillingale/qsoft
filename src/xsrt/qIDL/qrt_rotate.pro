FUNCTION QRT_ROTATE,is,pl,ax,angle

;   *IS        input        surface element number integer
;   *PL        input        position of rotation centre dblarr(3)
;   *AX        input        rotation axis dblarr(3)
;   *ANGLE     input        rotation angle degrees
;   
; note call_external writes to stdout not the IDL command line unit. 
; the output may be in th eterminal from which you started idlde if you are working in the IDE

qsoft = GETENV('QSOFT')
  A=CALL_EXTERNAL(qsoft + '/R_libraries/xsrtR.so', "qrt_rotate_", $
    LONG(is), $
    DOUBLE(pl), $
    DOUBLE(ax), $
    DOUBLE(angle), $
    /AUTO_GLUE)

RETURN, 'y'
END