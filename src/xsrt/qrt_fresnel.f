*+QRT_FRESNEL calculate reflectivity using Fresnel equations
        SUBROUTINE QRT_FRESNEL(NREAL,KIMAG,NANGS,ANGS,RS,RP,RUNP)
*NREAL     input        real part of refractive index
*KIMAG     input        imaginary part of refractive index
*NANGS     input        number of incidence angles
*ANGS      input        incidence angles (degrees range 0-90)
*RS        output        sigma reflectivity
*RP        output        pi reflectivity
*RUNP      output        unpolarized reflectivity
*Reference "Handbook of Optical Constants of Solids" Ed. Edward D.Palik
*Academic Press 1985, page 70
*If ANGS(I) out of range 0-90 degrees returns zero reflectivities
Cf2py  intent(in) NREAL,KIMAG,NANGS,ANGS
Cf2py  intent(out) RS,RP,RUNP
        IMPLICIT NONE 
        INTEGER NANGS
        DOUBLE PRECISION NREAL,KIMAG,ANGS(NANGS)
        DOUBLE PRECISION RS(NANGS),RP(NANGS),RUNP(NANGS)
*-Dick Willingale 2012-May-08
        DOUBLE PRECISION N,K
        DOUBLE PRECISION A2,A,PHI,SQBIT,B2
        DOUBLE PRECISION PI,BIT1
        PARAMETER(PI=3.1415926535898)
        INTEGER I
C
        N=NREAL
        K=KIMAG
        DO I=1,NANGS
            IF(ANGS(I).GE.0.0D0.AND.ANGS(I).LE.90.0D0) THEN
                      PHI=ANGS(I)*PI/180.0D0
                SQBIT=(N**2)-(K**2)-(SIN(PHI))**2
                BIT1=(SQBIT**2)+(4.0D0)*(N**2)*(K**2)
                A2=0.5D0*(SQRT(BIT1)+SQBIT)
                B2=0.5D0*(SQRT(BIT1)-SQBIT)
                A=SQRT(A2)                       
                RS(I)=(((A-COS(PHI))**2)+B2)/(((A+COS(PHI))**2)+B2)
                RP(I)=RS(I)*(((A-SIN(PHI)*TAN(PHI))**2)+B2)/
     &                (((A+SIN(PHI)*TAN(PHI))**2)+B2)
                IF(RS(I).GT.1.0D0) RS(I)=1.0D0
                IF(RP(I).GT.1.0D0) RP(I)=1.0D0
                RUNP(I)=0.5D0*(RP(I)+RS(I))
            ELSE
                RS(I)=0.0D0
                RP(I)=0.0D0
                RUNP(I)=0.0D0
            ENDIF
        ENDDO
        END
