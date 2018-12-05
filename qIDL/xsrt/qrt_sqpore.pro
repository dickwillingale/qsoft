FUNCTION QRT_SQPORE,pcen,pnorm,raxis,rcur,ipack,rap,pitch,wall,plen,idf,iq,plmin,plmax,fibre,ar

qsoft = GETENV('QSOFT')
  A=CALL_EXTERNAL(qsoft + '/R_libraries/xsrtR.so', "qrt_sqpore_", $
          DOUBLE(pcen),$
          DOUBLE(pnorm),$
          DOUBLE(raxis),$
          DOUBLE(rcur),$
          LONG(ipack),$
          DOUBLE(rap),$
          DOUBLE(pitch),$
          DOUBLE(wall),$
          DOUBLE(plen),$
          LONG(idf),$
          LONG(iq),$
          DOUBLE(plmin),$
          DOUBLE(plmax),$
          DOUBLE(fibre),$
          LONG(N_ELEMENTS(ar)),$
          DOUBLE(ar),$
    /AUTO_GLUE)

RETURN, A
END