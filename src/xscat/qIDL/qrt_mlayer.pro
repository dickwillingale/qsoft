FUNCTION QRT_MLAYER, angs,ekev,nr,ni,d,nper

          rsig=DBLARR(N_ELEMENTS(angs))
          rpi=DBLARR(N_ELEMENTS(angs))
          tsig=DBLARR(N_ELEMENTS(angs))
          tpi=DBLARR(N_ELEMENTS(angs))

qsoft = GETENV('QSOFT')
  A=CALL_EXTERNAL(qsoft + '/R_libraries/xscatR.so', "qrt_mlayer_", $
          LONG(N_ELEMENTS(angs)),$
          DOUBLE(angs),$
          DOUBLE(ekev),$
          LONG(N_ELEMENTS(nr)),$
          DOUBLE(nr),$
          DOUBLE(ni),$
          DOUBLE(d),$
          LONG(nper),$
          DOUBLE(rsig),$
          DOUBLE(rpi),$
          DOUBLE(tsig),$
          DOUBLE(tpi),$
    /AUTO_GLUE)

out = {rs:rs,rp:rp,runp:runp}

RETURN, out
END