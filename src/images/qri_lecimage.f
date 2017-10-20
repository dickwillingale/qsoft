*+QRI_LECIMAGE Create an image of the lobster eye cross beam
        SUBROUTINE QRI_LECIMAGE(S,H,B,XCEN,YCEN,NELS1,NELS2,ARRAY)
        IMPLICIT NONE
        INTEGER NELS1,NELS2
        DOUBLE PRECISION S,H,B,XCEN,YCEN,ARRAY(NELS1,NELS2)
Cf2py  intent(in) S,H,B,XCEN,YCEN,NELS1,NELS2
Cf2py  intent(out) ARRAY
*S           input        size of square area in pixels
*H           input        height of cross-arm triangle in pixels (=2d/L)
*B           input        width of cross-arm triangle in pixels
*XCEN        input        centre pixel (see coords below)
*YCEN        input        centre pixel (see coords below)
*NELS1       input        first dimension of array
*NELS2       input        second dimension of array
*ARRAY       output       image array
*Pixels assumed to run:
*        X 1 left to NELS1 right
*        Y 1 bottom to NELS2 top
*Coordinate system for XCEN,YCEN and other position is:
*        X runs from 0.0 on left to NELS1 on right
*        Y runs from 0.0 on bottom to NELS2 on top
*Centre of bottom left pixel is therefore 0.5,0.5
*Centre of top right pixel is NELS1-0.5,NELS2-0.5
*-Author Dick Willingale 2017-Sept-16
        INCLUDE 'QR_COM'
        DOUBLE PRECISION XP,YP,XR,YR,BB,BB2,BH,S2,SM,A
        INTEGER J,K
C
        IF(ISTAT.NE.0) RETURN
C Set constants
        A=1.0/H-B*0.5/H**2
        BB2=(1.0-SQRT(1.0-4.0*A*B*0.5))/A/2.0
	BB=BB2*2.0
        S2=S*0.5
        BH=BB2/H
        SM=MIN(S2,H)
C Loop around pixels
        DO J=1,NELS2
                YP=(DBLE(J)-0.5)
                YR=ABS(YP-YCEN)
                DO K=1,NELS1
                        ARRAY(K,J)=0.0
                        XP=(DBLE(K)-0.5)
                        XR=ABS(XP-XCEN)
                        IF((XR.LT.SM).AND.(YR.LT.SM)) THEN
                                IF((YR.LT.(BB2-BH*XR)).OR.
     +                          (XR.LT.(BB2-BH*YR))) THEN
                                    ARRAY(K,J)=1.0
                                ENDIF
                        ENDIF
                ENDDO
        ENDDO
        END
