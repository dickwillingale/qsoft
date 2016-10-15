*+QRT_SURFACE Set surface quality parameters
        SUBROUTINE QRT_SURFACE(IS,IT,EKEV,SRGH,FMIN,PIND,ALPHA,GAMMA,
     +  NREF,ANGS,REFS,GPITCH,DHUB,ORDER)
*IS      input    surface quality index
*IT      input    surface type
*             1 refl. (Fresnel), 2 refl. (look-up), 3 refract, 4 diffract
*EKEV    input    photon energy keV
*SRGH    input    Specific roughness (A**2 mm)
*             if -ve then rms figure gradient error radians
*FMIN    input    Minimum surface spatial frequency
*PIND    input    Surface roughness power spectrum index
*ALPHA   input    real part of diel. constant or refractive index ratio N1/N2
*GAMMA   input    imaginary part of dielectric constant
*NREF    input    number of ANGS and REFS pairs
*ANGS    input    incidence angles (degrees increasing) (QTYPE 2 and 4)
*REFS    input    reflectivity values (QTYPE 2 and 4)
*GPITCH  input    d-spacing for grating mm (QTYPE 4)
*DHUB    input    distance from ruling hub to surface reference point (QTYPE 4)
*             if <1.0 mm then d-spacing gradient across ruling (QTYPE 4)
*ORDER   input    diffraction order (QTYPE 4)
Cf2py  intent(in) IS,IT,EKEV,SRGH,FMIN,PIND,ALPHA,GAMMA
Cf2py  intent(in) NREF,ANGS,REFS,GPITCH,DHUB,ORDER
        IMPLICIT NONE
        INTEGER IS,IT,NREF
        DOUBLE PRECISION EKEV,SRGH,FMIN,PIND,ANGS(*),REFS(*)
        DOUBLE PRECISION ALPHA,GAMMA,GPITCH,DHUB,ORDER
*-Author Dick Willingale 2012-Apr-30
        INTEGER NP
        DOUBLE PRECISION AKEV,WLNG
        PARAMETER (AKEV=12.397639)
        INCLUDE 'QR_COM'
C
        IF(ISTAT.NE.0) RETURN
C
        NP=0
        IF(IT.EQ.2.OR.IT.EQ.4) THEN
                NP=NREF
        ENDIF
        WLNG=AKEV/EKEV
        CALL SRT_SETT(IS,IT,WLNG,SRGH,FMIN,PIND,GPITCH,DHUB,ORDER,
     +        ALPHA,GAMMA,NP,ANGS,REFS,ISTAT)
        END
