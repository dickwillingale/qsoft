FUNCTION QRI_INIT

qsoft = GETENV('QSOFT')
  A=CALL_EXTERNAL(qsoft + '/R_libraries/imagesR.so', "qri_init_", $
    /AUTO_GLUE)

RETURN, 'Initialised QRI'
END