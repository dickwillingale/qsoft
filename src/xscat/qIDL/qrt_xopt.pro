FUNCTION QRT_XOPT,mspec,rho,ekev,itype

  alpha=DBLARR(N_ELEMENTS(ekev))
  gammaa=DBLARR(N_ELEMENTS(ekev))
  absl=DBLARR(N_ELEMENTS(ekev))
  f1=DBLARR(N_ELEMENTS(ekev))
  f2=DBLARR(N_ELEMENTS(ekev))


; format string mspec into a byte array as IDL strings are difficult to pass.
chars = BYTE(mspec)

qsoft = GETENV('QSOFT')
  A=CALL_EXTERNAL(qsoft + '/R_libraries/xscatR.so', "qrt_xopt_", $
          LONG(STRLEN(mspec)),$
          chars,$
          DOUBLE(rho),$
          LONG(N_ELEMENTS(ekev)),$
          DOUBLE(ekev),$
          LONG(itype),$
          alpha,$
          gammaa,$
          absl,$
          f1,$
          f2,$
    /AUTO_GLUE)

out = {alpha:alpha,gamma:gammaa,absl:absl,f1:f1,f2:f2}

RETURN, out
END