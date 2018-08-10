*+QRI_SETFIELD        Set image field
        SUBROUTINE QRI_SETFIELD(NX,XLEFT,XRIGHT,NY,YBOT,YTOP)
        IMPLICIT NONE
        INTEGER NX,NY
        DOUBLE PRECISION XLEFT,XRIGHT,YBOT,YTOP
Cf2py   intent(in) NX,XLEFT,XRIGHT,NY,YBOT,YTOP
*NX        input        number of columns
*XLEFT     input        local x left edge
*XRIGHT    input        local x right edge
*NY        input        number of rows
*YBOT      input        local y bottom edge
*YTOP      input        local y top edge
*-Author Dick Willingale 2012-May-11
        INCLUDE 'QRI_TRANSCOM'
        INCLUDE 'QR_COM'
C
        IF(ISTAT.NE.0) RETURN
C Set image field parameters into common
        NBYM(1)=NX
        NBYM(2)=NY
        BLHC(1)=XLEFT
        BLHC(2)=YBOT
        XYPIXEL(1)=(XRIGHT-XLEFT)/NX
        XYPIXEL(2)=(YTOP-YBOT)/NY
        XYNOW(1)=(XRIGHT+XLEFT)/2.0
        XYNOW(2)=(YTOP+YBOT)/2.0
        END
