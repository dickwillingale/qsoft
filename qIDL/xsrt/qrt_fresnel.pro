FUNCTION QRT_FRESNEL,nreal,kimag,angs

;  *NREAL     input        real part of refractive index
;  *KIMAG     input        imaginary part of refractive index
;  *NANGS     input        number of incidence angles
;  *ANGS      input        incidence angles (degrees range 0-90)
;  *RS        output        sigma reflectivity
;  *RP        output        pi reflectivity
;  *RUNP      output        unpolarized reflectivity

  rs=DBLARR(N_ELEMENTS(angs))
  rp=DBLARR(N_ELEMENTS(angs))
  runp=DBLARR(N_ELEMENTS(angs))

qsoft = GETENV('QSOFT')
  A=CALL_EXTERNAL(qsoft + '/R_libraries/xsrtR.so', "qrt_fresnel_", $
        DOUBLE(nreal),$
        DOUBLE(kimag),$
        LONG(N_ELEMENTS(angs)),$
        DOUBLE(angs),$
        DOUBLE(rs),$
        DOUBLE(rp),$
        DOUBLE(runp),$
    /AUTO_GLUE)
    
out = {rs:rs,rp:rp,runp:runp}

RETURN, out
END