FUNCTION QRT_SURFACE,is,it,ekev,srgh,fmin,pind,alpha,gammaa,angs,refs,gpitch,dhub,order

;  *IS      input    surface quality index
;  *IT      input    surface type
;  *             1 refl. (Fresnel), 2 refl. (look-up), 3 refract, 4 diffract
;  *EKEV    input    photon energy keV
;  *SRGH    input    Specific roughness (A**2 mm)
;  *             if -ve then rms figure gradient error radians
;  *FMIN    input    Minimum surface spatial frequency
;  *PIND    input    Surface roughness power spectrum index
;  *ALPHA   input    real part of diel. constant or refractive index ratio N1/N2
;  *GAMMA   input    imaginary part of dielectric constant
;  *NREF    input    number of ANGS and REFS pairs
;  *ANGS    input    incidence angles (degrees increasing) (QTYPE 2 and 4)
;  *REFS    input    reflectivity values (QTYPE 2 and 4)
;  *GPITCH  input    d-spacing for grating mm (QTYPE 4)
;  *DHUB    input    distance from ruling hub to surface reference point (QTYPE 4)
;  *             if <1.0 mm then d-spacing gradient across ruling (QTYPE 4)
;  *ORDER   input    diffraction order (QTYPE 4)

qsoft = GETENV('QSOFT')
  A=CALL_EXTERNAL(qsoft + '/R_libraries/xsrtR.so', "qrt_surface_", $
        LONG(is),$
        LONG(it),$
        DOUBLE(ekev),$
        DOUBLE(srgh),$
        DOUBLE(fmin),$
        DOUBLE(pind),$
        DOUBLE(alpha),$
        DOUBLE(gammaa),$
        LONG(N_ELEMENTS(angs)),$
        DOUBLE(angs),$
        DOUBLE(refs),$
        DOUBLE(gpitch),$
        DOUBLE(dhub),$
        DOUBLE(order),$
    /AUTO_GLUE)
    

RETURN, a
END