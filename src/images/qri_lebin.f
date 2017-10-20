*+QRI_LEBIN Create an image from an event list using lobster eye psf
        SUBROUTINE QRI_LEBIN(NE,XE,YE,S,H,G,ETA,NELS1,NELS2,ARRAY)
        IMPLICIT NONE
        INTEGER NE,NELS1,NELS2
        DOUBLE PRECISION XE(NE),YE(NE),S,H,G,ETA,ARRAY(NELS1,NELS2)
Cf2py  intent(in) NE,XE,YE,S,H,G,ETA,NELS1,NELS2
Cf2py  intent(out) ARRAY
*NE          input        number of events
*XE          input        x event positions, pixels
*YE          input        y event positions, pixels
*S           input        size of square area in pixels
*H           input        height of cross-arm triangle in pixels (=2d/L)
*G           input        width of Lorentzian central spot pixels
*ETA         input        cross-arms to peak ratio at centre
*NELS1       input        first dimension of array
*NELS2       input        second dimension of array
*ARRAY       output       image array
*Event positions assumed to run:
*        X 0 left edge to NELS1 right edge
*        Y 0 bottom edge to NELS2 top edge
*Therefore centre of left bottom pixel is 0.5,0.5
*The binning is effectively a cross-correlation of the event list with
*the lobster eye cross-beam. 
*-Author Dick Willingale 2017-Sept-17
        INCLUDE 'QR_COM'
        DOUBLE PRECISION XR,YR,G2,S2,SM,XB,YB,EG
        INTEGER I,J,K
C
        IF(ISTAT.NE.0) RETURN
C Set beam constants
        S2=S*0.5
        G2=G*0.5
        SM=MIN(S2,H)
        EG=ETA*G/H
C Loop around pixels of output image
        DO J=1,NELS2
         DO K=1,NELS1
          ARRAY(K,J)=0.0
C Loop around events
          DO I=1,NE
           XR=ABS(DBLE(K)-0.5-XE(I))
           YR=ABS(DBLE(J)-0.5-YE(I))
           IF((XR.LT.SM).AND.(YR.LT.SM)) THEN
               XB=1.0/((XR/G2)**2+1.0)+EG*(1.0-(XR/H)**2)
               YB=1.0/((YR/G2)**2+1.0)+EG*(1.0-(YR/H)**2)
               ARRAY(K,J)=ARRAY(K,J)+XB*YB/(1.0+EG)**2
           ENDIF
          ENDDO
         ENDDO
        ENDDO
        END
