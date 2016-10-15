*+QRT_MLAYER Calculate X-ray reflectivities of a multilayer
        SUBROUTINE QRT_MLAYER(NG,ANGS,EKEV,NLAY,NR,NI,D,NPER,
     +  RSIG,RPI,TSIG,TPI)
*NG       input        number of incidence angles
*ANGS     input        incidence angles (degrees)
*EKEV     input        X-ray energy keV
*NLAY     input        number of layers (including substrate if there is one)
*NR       input        real part of refractive index of layers
*NI       input        imaginary part of refractive index of layers
*D        input        thickness of layers A
*NPER     input        number of periods
*RSIG     output        reflectivity sigma polarization
*RPI      output        reflectivity pi polarization
*TSIG     output        transmission sigma polarization
*TPI      output        transmission pi polarization
Cf2py  intent(in) NG,ANGS,EKEV,NLAY,NR,NI,D,NPER
Cf2py  intent(out) RSIG,RPI,TSIG,TPI
        IMPLICIT NONE
        INTEGER NG,NLAY,NPER
        DOUBLE PRECISION ANGS(NG),EKEV,NR(NLAY),NI(NLAY),D(NLAY)
        DOUBLE PRECISION RSIG(NG),RPI(NG),TSIG(NG),TPI(NG)
*-Author Dick Willingale 2012-May-08
        INCLUDE 'QR_COM'
        DOUBLE PRECISION ANG
        DOUBLE PRECISION AKEV,WL
        PARAMETER (AKEV=12.397639)
        INTEGER J,NMAX
        PARAMETER (NMAX=400)
        DOUBLE COMPLEX FR(NMAX)
        DOUBLE PRECISION PI,PIBY2
        PARAMETER (PIBY2=1.5707963267949,PI=PIBY2*2.D0)
C 
        IF(ISTAT.NE.0) RETURN
C
        IF(NMAX.LT.NLAY) THEN
                WRITE(*,*) 'QRT_MLAYER - too many layers'
                RETURN
        ENDIF
C Set complex refractive index
        DO J=1,NLAY
                FR(J)=DCMPLX(NR(J),NI(J))
        ENDDO
C Loop for incidence angles
        WL=AKEV/EKEV
        DO J=1,NG
C Convert angle to radians
             ANG=ANGS(J)*PI/180.0
             CALL SRT_MLTI(ANG,WL,NLAY,FR,D,NPER,RSIG(J),RPI(J),
     +       TSIG(J),TPI(J),ISTAT)
        ENDDO        
        END
