FUNCTION QRS_INIT

qsoft = GETENV('QSOFT')
  A=CALL_EXTERNAL(qsoft + '/R_libraries/xscatR.so', "qrs_init_", $
    /AUTO_GLUE)

RETURN, 'Initialised QRS'
END