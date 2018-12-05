FUNCTION QRT_DEFORM,id,it,nm,nx,ny

;  *ID        input        deformation index
;  *IT        input        deformation type (1 matrix)
;  *NM        input        number of sub-matrices
;  *NX        input        number of x samples
;  *NY        input        number of y samples

qsoft = GETENV('QSOFT')
  A=CALL_EXTERNAL(qsoft + '/R_libraries/xsrtR.so', "qrt_defs_", $
        LONG(id),$
        LONG(it),$
        LONG(nm),$
        LONG(nx),$
        LONG(ny),$   
    /AUTO_GLUE)

RETURN, A
END