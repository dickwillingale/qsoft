FUNCTION QRT_DEFMAT,id,im,x,y,z

;  *ID        input        deformation index
;  *IM        input        sub-matrix index
;  *KNX       input        number of x samples
;  *KNY       input        number of y samples
;  *XSAM      input        x values
;  *YSAM      input        y values
;  *ZDEF      input        deformation matrix
;  *NW        input        size of work arrays (at least MAX(KNX,KNY))
;  *W1        input        work array
 ; *W2        input        work array

qsoft = GETENV('QSOFT')
  A=CALL_EXTERNAL(qsoft + '/R_libraries/xsrtR.so', "qrt_mats_", $
          LONG(id),$
          LONG(im),$
          LONG(N_ELEMENTS(x)),$
          LONG(N_ELEMENTS(y)),$
          DOUBLE(x),$
          DOUBLE(y),$
          DOUBLE(z),$
          LONG(N_ELEMENTS(x)*N_ELEMENTS(y)),$
          DBLARR(N_ELEMENTS(x)*N_ELEMENTS(y)),$
          DBLARR(N_ELEMENTS(x)*N_ELEMENTS(y)),$
    /AUTO_GLUE)

RETURN, A
END