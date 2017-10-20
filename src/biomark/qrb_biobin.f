*+QRB_BIOBIN raw spectrum binning
        SUBROUTINE QRB_BIOBIN(N,X,W,XLEFT,XRIGHT,NX,A,XA)
        IMPLICIT NONE
        INTEGER N,NX
        DOUBLE PRECISION X(N),W(N),XLEFT,XRIGHT,A(NX),XA(NX)
Cf2py    intent(in) N,X,W,XLEFT,XRIGHT,NX
Cf2py    intent(out) A,XA
*N        input        number of raw values
*X        input        array of x positions
*W        input        array of weights
*XLEFT    input        left-hand side of 1st bin
*XRIGHT   input        right-hand side of last bin
*NX       input        dimension of output array (number of bins)
*A        output       array
*XA       output       array of bin centre positions
*-Author Dick Willingale 2014-Mar-27
* qsoft version 2017-Jun-26
        INTEGER I,IX
        DOUBLE PRECISION XB
C Calculate bin size
        XB=(XRIGHT-XLEFT)/NX
C Clear destination array and set bin positions
        DO I=1,NX
                A(I)=0.0
                XA(I)=XLEFT+XB*(I-1)+XB*0.5
        ENDDO
C Loop through events
        DO I=1,N
C Calculate x bin
                IX=INT((X(I)-XLEFT)/XB)
                IX=IX+1
                IF(IX.GT.0.AND.IX.LE.NX) THEN
                        A(IX)=A(IX)+W(I)
                ENDIF
        ENDDO
        END
