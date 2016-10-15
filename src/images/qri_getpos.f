*+QRI_GETPOS        Get current position in image field
        SUBROUTINE QRI_GETPOS(PIX,XYL,AES,EQU,ECL,GAL)
        IMPLICIT NONE
        DOUBLE PRECISION PIX(2),XYL(2),AES(2),EQU(2),ECL(2),GAL(2)
Cf2py   intent(out) PIX,XYL,AES,EQU,ECL,GAL
*pix        output        pixel position 0-NCOLS, 0-NROWS
*xyl        output        local X,Y
*aes        output        local azimuth,elevation degrees
*equ        output        Equatorial RA,DEC degrees J2000
*ecl        output        Ecliptic EA,EL degrees
*gal        output        Galactic LII,BII degrees
*-Author Dick Willingale 2012-May-11
        INCLUDE 'QRI_TRANSCOM'
        INCLUDE 'QR_COM'
        DOUBLE PRECISION DTOR,PT(2),CL(2),XY(2)
C
        IF(ISTAT.NE.0) RETURN
C
        DTOR=ASIN(1.0D0)/90.D0
C
        XY(1)=XYNOW(1)
        XY(2)=XYNOW(2)
        XYL(1)=XY(1)
        XYL(2)=XY(2)
        CALL QRI_LTOP(XY,PT)
        PIX(1)=PT(1)
        PIX(2)=PT(2)
        CALL QRI_LTOS(XY,CL)
        AES(1)=CL(1)/DTOR
        AES(2)=CL(2)/DTOR
        CALL QRI_STOC(CL,PT)
        EQU(1)=PT(1)/DTOR
        EQU(2)=PT(2)/DTOR
        CALL QRI_CTOE(PT,CL)
        ECL(1)=CL(1)/DTOR
        ECL(2)=CL(2)/DTOR
        CALL QRI_CTOG(PT,CL)
        GAL(1)=CL(1)/DTOR
        GAL(2)=CL(2)/DTOR
        END
