*+QRI_ANNULUS        Calculate average level and scatter of pixels in annulus
        SUBROUTINE QRI_ANNULUS(NEL1,NEL2,ARRAY,RMIN,RMAX,
     +        NSAM,BLEV,BVAR)
*NEL1        input        x dimension of array
*NEL2        input        y dimension of array
*ARRAY        input        data array
*RMIN        input        minimum radius
*RMAX        input        maximum radius
*NSAM        output        number of pixels in box
*BLEV        output        average level per pixel
*BVAR        output        variance on BLEV
        IMPLICIT NONE
        INTEGER NEL1,NEL2,NSAM
        DOUBLE PRECISION ARRAY(NEL1,NEL2),RMIN,RMAX,BLEV,BVAR
*-Author Dick Willingale 2012-May-14
        INCLUDE 'QRI_TRANSCOM'
        INTEGER J,K,IXP,IYP,IRAD,NXL,NXH,NYL,NYH
        DOUBLE PRECISION XP,YP,RAD,VAL,BLEVSQ,XAN,YAN,P(2),XY(2)
        INCLUDE 'QR_COM'
C 
        IF(ISTAT.NE.0) RETURN
C Set position
        XY(1)=XYNOW(1)
        XY(2)=XYNOW(2)
        CALL QRI_LTOP(XY,P)
        XAN=P(1)
        YAN=P(2)
C Initialize results
        BLEV=0.0
        BVAR=0.0
        NSAM=0
C Find pixel ranges about centre
        IXP=XAN+1.0
        IYP=YAN+1.0
        IF(IXP.LT.1.OR.IXP.GT.NEL1.OR.IYP.LT.1.OR.IYP.GT.NEL2) THEN
                RETURN
        ENDIF
        IRAD=RMAX+2.0
        NXL=MAX(MIN(IXP-IRAD,NEL1),1)
        NXH=MAX(MIN(IXP+IRAD,NEL1),1)
        NYL=MAX(MIN(IYP-IRAD,NEL2),1)
        NYH=MAX(MIN(IYP+IRAD,NEL2),1)
C
        DO J=NYL,NYH
                YP=(REAL(J)-0.5)-YAN
                DO K=NXL,NXH
                        XP=(REAL(K)-0.5)-XAN
                        RAD=SQRT(XP**2+YP**2)
                        IF(RAD.GT.RMIN.AND.RAD.LE.RMAX) THEN
                                NSAM=NSAM+1
                                VAL=ARRAY(K,J)
                                BLEV=BLEV+VAL
                                BVAR=BVAR+VAL**2
                        ENDIF
                ENDDO
        ENDDO
        IF(NSAM.GT.0) THEN
                BLEV=BLEV/REAL(NSAM)
                BLEVSQ=BLEV**2
                BVAR=BVAR/REAL(NSAM)
                BVAR=BVAR-BLEVSQ
        ENDIF
        END
