FUNCTION QRT_INIT

qsoft = GETENV('QSOFT')
  A=CALL_EXTERNAL(qsoft + '/R_libraries/xsrtR.so', "qrt_init_", $
    /AUTO_GLUE)

RETURN, 'Initialised QRT'
END