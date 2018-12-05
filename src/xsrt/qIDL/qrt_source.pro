FUNCTION QRT_SOURCE,it,sd,sp,ap,an,ar,al,apry,nr,idef

;  *IT        input        source type
;  *                1 point source at infinity               radial limits
;  *                2 point source at infinity               cartesian limits
;  *                3 point source at finite distance        radial limits
;  *                4 point source at finite distance        cartesian limits
;  *SD        input        source direction cosines
;  *SP        input        source position
;  *AP        input        aperture position
;  *AN        input        aperture normal
;  *AR        input        aperture reference axis
;  *AL        input        aperture limits
;  *APRY      input        area per ray - if 0 then use NRAY
;  *NR        input        number of rays
;  *IDEF      input        deformation index


qsoft = GETENV('QSOFT')
  A=CALL_EXTERNAL(qsoft + '/R_libraries/xsrtR.so', "qrt_source_", $
        LONG(it),$
        DOUBLE(sd),$
        DOUBLE(sp),$
        DOUBLE(ap),$
        DOUBLE(an),$
        DOUBLE(ar),$
        DOUBLE(al),$
        DOUBLE(apry),$
        LONG(nr),$
        LONG(idef),$
    /AUTO_GLUE)
    


RETURN, a
END