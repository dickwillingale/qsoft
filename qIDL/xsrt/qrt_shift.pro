FUNCTION QRT_SHIFT,is,pl

;  *IS        input        surface element number integer
;  *PL        input        vector shift dblarr(3)

; note call_external writes to stdout not the IDL command line unit. 
; the output may be in th eterminal from which you started idlde if you are working in the IDE

qsoft = GETENV('QSOFT')
  A=CALL_EXTERNAL(qsoft + '/R_libraries/xsrtR.so', "qrt_shift_", $
    LONG(is), $
    DOUBLE(pl), $
    /AUTO_GLUE)

RETURN, 'y'
END