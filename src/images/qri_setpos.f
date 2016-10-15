*+QRI_SETPOS        Set current position for image field
        SUBROUTINE QRI_SETPOS(IPOS,P)
        IMPLICIT NONE
        INTEGER IPOS
        DOUBLE PRECISION P(2)
Cf2py   intent(in) IPOS,P
*IPOS     input  1 pixel 0-NCOLS, 0-NROWS
*                2 local X,Y
*                3 local azimuth,elevation degrees
*                4 Celestial RA,DEC degrees J2000
*                5 Ecliptic EA,EL degrees
*                6 Galactic LII,BII degrees
*P        input        position
*-Author Dick Willingale 2012-May-11
        INCLUDE 'QRI_TRANSCOM'
        INCLUDE 'QR_COM'
        DOUBLE PRECISION DTOR,POUT(2),CEL(2),XY(2)
C
        IF(ISTAT.NE.0) RETURN
C
        DTOR=ASIN(1.0D0)/90.D0
C
        IF(IPOS.EQ.1) THEN
                CALL QRI_PTOL(P,XY)
                XYNOW(1)=XY(1)
                XYNOW(2)=XY(2)
                CALL QRI_LTOS(XY,CEL)
                AENOW(1)=CEL(1)
                AENOW(2)=CEL(2)
        ELSEIF(IPOS.EQ.2) THEN
                XYNOW(1)=P(1)
                XYNOW(2)=P(2)
                CALL QRI_LTOS(P,CEL)
                AENOW(1)=CEL(1)
                AENOW(2)=CEL(2)
        ELSEIF(IPOS.EQ.3) THEN
                CEL(1)=P(1)*DTOR
                CEL(2)=P(2)*DTOR
                AENOW(1)=CEL(1)
                AENOW(2)=CEL(1)
                CALL QRI_STOL(CEL,XY)
                XYNOW(1)=XY(1)
                XYNOW(2)=XY(2)
        ELSEIF(IPOS.EQ.4) THEN        
                CEL(1)=P(1)*DTOR
                CEL(2)=P(2)*DTOR
                CALL QRI_CTOS(CEL,POUT)
                AENOW(1)=POUT(1)
                AENOW(2)=POUT(2)
                CALL QRI_STOL(POUT,XY)
                XYNOW(1)=XY(1)
                XYNOW(2)=XY(2)
        ELSEIF(IPOS.EQ.5) THEN        
                CEL(1)=P(1)*DTOR
                CEL(2)=P(2)*DTOR
                CALL QRI_ETOC(CEL,POUT)
                CALL QRI_CTOS(POUT,CEL)
                AENOW(1)=CEL(1)
                AENOW(2)=CEL(2)
                CALL QRI_STOL(CEL,XY)
                XYNOW(1)=XY(1)
                XYNOW(2)=XY(2)
        ELSEIF(IPOS.EQ.6) THEN        
                CEL(1)=P(1)*DTOR
                CEL(2)=P(2)*DTOR
                CALL QRI_GTOC(CEL,POUT)
                CALL QRI_CTOS(POUT,CEL)
                AENOW(1)=CEL(1)
                AENOW(2)=CEL(2)
                CALL QRI_STOL(CEL,XY)
                XYNOW(1)=XY(1)
                XYNOW(2)=XY(2)
        ENDIF
        END
