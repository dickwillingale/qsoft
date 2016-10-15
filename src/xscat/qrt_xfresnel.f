*+QRT_XFRESNEL Calculate X-ray reflectivity using Fresnel equations
        SUBROUTINE QRT_XFRESNEL(ALPHA,GAMMA,NANGS,ANGS,RS,RP,RUNP)
*ALPHA     input        real incremental part dielectric constant
*GAMMA     input        imaginary part of dielectric constant
*NANGS     input        number of incidence angles
*ANGS      input        incidence angles (degrees)
*RS        output       sigma reflectivity
*RP        output       pi reflectivity
*RUNP      output       unpolarized reflectivity
*If ANGS(I) out of range 0-90 degrees returns zero reflectivities
Cf2py  intent(in) ALPHA,GAMMA,NANGS,ANGS
Cf2py  intent(out) RS,RP,RUNP
        IMPLICIT NONE 
        INTEGER NANGS
        DOUBLE PRECISION ALPHA,GAMMA,ANGS(NANGS)
        DOUBLE PRECISION RS(NANGS),RP(NANGS),RUNP(NANGS)
*-Dick Willingale 2012-May-08
        DOUBLE PRECISION PI,THETAG
        PARAMETER(PI=3.1415926535898)
        INTEGER I
C
        DO I=1,NANGS
            IF(ANGS(I).GE.0.0D0.AND.ANGS(I).LE.90.0D0) THEN
                THETAG=(90.0D0-ANGS(I))*PI/180.D0
                CALL SRT_FRNL(THETAG,ALPHA,GAMMA,RS(I),RP(I))
                RUNP(I)=0.5D0*(RP(I)+RS(I))
            ELSE
                RS(I)=0.0D0
                RP(I)=0.0D0
                RUNP(I)=0.0D0
            ENDIF
        ENDDO
        END
