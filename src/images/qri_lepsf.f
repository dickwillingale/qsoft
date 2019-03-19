*+QRI_LEPSF Create an image of the lobster eye PSF
        SUBROUTINE QRI_LEPSF(S,H,G,ETA,ALP,XCEN,YCEN,NELS1,NELS2,ARRAY)
        IMPLICIT NONE
        INTEGER NELS1,NELS2
        DOUBLE PRECISION S,H,G,ETA,ALP,XCEN,YCEN,ARRAY(NELS1,NELS2)
Cf2py  intent(in) S,H,G,ETA,ALP,XCEN,YCEN,NELS1,NELS2
Cf2py  intent(out) ARRAY
*S           input        size of square patch area in pixels
*H           input        height of cross-arm triangle in pixels (=2d/L)
*G           input        width of Lorentzian central spot pixels
*ETA         input        cross-arm to peak ratio at centre
*ALP         input        index (1 for Lorenztian)
*XCEN        input        centre of PSF (see below for coords. system)
*YCEN        input        centre of PSF (see below for coords. system)
*NELS1       input        first dimension of array
*NELS2       input        second dimension of array
*ARRAY       output       image array
*Pixels indices assumed to run:
*        X 1 left to NELS1 right
*        Y 1 bottom to NELS2 top
*Coordinate system for XCEN,YCEN and other position is:
*        X runs from 0.0 on left to NELS1 on right
*        Y runs from 0.0 on bottom to NELS2 on top
*Centre of bottom left pixel is therefore 0.5,0.5
*Centre of top right pixel is NELS1-0.5,NELS2-0.5
*-Author Dick Willingale 2017-Sept-24
        INCLUDE 'QR_COM'
        DOUBLE PRECISION XP,YP,XR,YR,G2,S2,XB,YB,SM,EG
        INTEGER J,K
C
        IF(ISTAT.NE.0) RETURN
C Set constants
        S2=S*0.5
        G2=G*0.5
        SM=MIN(S2,H)
        EG=ETA*G/H
C Loop around pixels
        DO J=1,NELS2
          YP=(DBLE(J)-0.5)
          YR=ABS(YP-YCEN)
          DO K=1,NELS1
            XP=(DBLE(K)-0.5)
            XR=ABS(XP-XCEN)
            IF((XR.LT.SM).AND.(YR.LT.SM)) THEN
              XB=1.0/((XR/G2)**2+1.0)**ALP+EG*(1.0-(XR/H)**2)
              YB=1.0/((YR/G2)**2+1.0)**ALP+EG*(1.0-(YR/H)**2)
              ARRAY(K,J)=XB*YB/(1.0+EG)**2
            ELSE
              ARRAY(K,J)=0.0
            ENDIF
          ENDDO
        ENDDO
        END
