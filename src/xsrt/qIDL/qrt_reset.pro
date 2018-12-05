FUNCTION QRT_RESET

; note call_external writes to stdout not the IDL command line unit. 
; the output may be in th eterminal from which you started idlde if you are working in the IDE

qsoft = GETENV('QSOFT')
  A=CALL_EXTERNAL(qsoft + '/R_libraries/xsrtR.so', "qrt_init_", $
    /AUTO_GLUE)

RETURN, 'QRT Reset'
END