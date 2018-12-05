FUNCTION QRT_XFRESNEL, alpha,gammaa,angs

  rs=DBLARR(N_ELEMENTS(angs))
  rp=DBLARR(N_ELEMENTS(angs))
  runp=DBLARR(N_ELEMENTS(angs))

qsoft = GETENV('QSOFT')
  A=CALL_EXTERNAL(qsoft + '/R_libraries/xscatR.so', "qrt_xfresnel_", $
          DOUBLE(alpha),$
          DOUBLE(gammaa),$
          LONG(N_ELEMENTS(angs)),$
          DOUBLE(angs),$
          DOUBLE(rs),$
          DOUBLE(rp),$
          DOUBLE(runp),$
    /AUTO_GLUE)

out = {rs:rs,rp:rp,runp:runp}

RETURN, out
END