*+QRI_BINXY        x-y event binning to form an image
        SUBROUTINE QRI_BINXY(N,X,Y,NQ,IQ,NW,W,XLEFT,XRIGHT,YBOT,YTOP,
     +  NX,NY,A)
*N        input        number of events
*X        input        array of x positions
*Y        input        array of y positions
*NQ       input        number of quality values (1 or N)
*IQ       input        array of quality values (0 for OK)
*NW       input        number of weights (1 or N)
*W        input        array of weights
*XLEFT    input        minimum x value
*XRIGHT   input        maximum x value
*YBOT     input        minimum y value
*YTOP     input        maximum y value
*NX,NY    input        dimensions of output array
*A        output        array
Cf2py  intent(in) N,X,Y,NQ,IQ,NW,W,XLEFT,XRIGHT,YBOT,YTOP,NX,NY
Cf2py  intent(out) A
        IMPLICIT NONE
        INTEGER N,NQ,NX,NY,NW
        INTEGER IQ(NQ)
        DOUBLE PRECISION X(N),Y(N),XLEFT,XRIGHT,YBOT,YTOP
        DOUBLE PRECISION A(NX,NY),W(NW)
*-Author Dick Willingale 2012-May-10
        INTEGER J,K,I,IX,IY,II
        DOUBLE PRECISION XB,YB
C Clear destination array
        DO J=1,NY
                DO K=1,NX
                        A(K,J)=0.0
                ENDDO
        ENDDO
C Calculate bin sizes
        XB=(XRIGHT-XLEFT)/NX
        YB=(YTOP-YBOT)/NY
C Loop through events
        DO I=1,N
C Check quality
                IF(NQ.EQ.1) THEN
                        II=1
                ELSE
                        II=I
                ENDIF
                IF(IQ(II).EQ.0) THEN
C Calculate x bin
                        IX=INT((X(I)-XLEFT)/XB)
                        IX=IX+1
                        IF(IX.GT.0.AND.IX.LE.NX) THEN
C Calculate y bin
                                IY=INT((Y(I)-YBOT)/YB)
                                IY=IY+1
                                IF(IY.GT.0.AND.IY.LE.NY) THEN
                                        IF(NW.EQ.1) THEN
                                                A(IX,IY)=A(IX,IY)+W(1)
                                        ELSE
                                                A(IX,IY)=A(IX,IY)+W(I)
                                        ENDIF
                                ENDIF
                        ENDIF
                ENDIF
        ENDDO
        END
