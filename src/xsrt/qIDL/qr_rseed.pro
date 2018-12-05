FUNCTION QR_RSEED,iseed

qsoft = GETENV('QSOFT')
  A=CALL_EXTERNAL(qsoft + '/R_libraries/xsrtR.so', "qr_rseed_", $
    LONG(iseed),$
    /AUTO_GLUE)

RETURN, 'y'
END