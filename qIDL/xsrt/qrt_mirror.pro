FUNCTION QRT_MIRROR,id,idf,iq,anml,arfx,apos,alim

;  *ID        input        aperture type
;  *IDF       input        deformation index
;  *IQ        input        surface quality index
;  *ANML      input        surface normal
;  *ARFX      input        surface reference axis
;  *APOS      input        surface reference position
;  *N         input        number of limits
;  *ALIM      input        limits array
;  *NSURF     input        number of subsequent surfaces ID=2

qsoft = GETENV('QSOFT')
  A=CALL_EXTERNAL(qsoft + '/R_libraries/xsrtR.so', "qrt_mirror_", $
        LONG(id),$
        LONG(idf),$
        LONG(iq),$
        DOUBLE(anml),$
        DOUBLE(arfx),$
        DOUBLE(apos),$
        LONG(N_ELEMENTS(alim)),$
        DOUBLE(alim),$
    /AUTO_GLUE)
    

RETURN, a
END