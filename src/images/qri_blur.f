*+qri_blur   crosscorrelation of image array with mask
        SUBROUTINE QRI_BLUR(NARX,NARY,ARR,NMASX,NMASY,MASK,RES)
        IMPLICIT NONE
        INTEGER NARX,NARY,NMASX,NMASY
        DOUBLE PRECISION ARR(NARX,NARY),MASK(NMASX,NMASY),RES(NARX,NARY)
*NARX,Y        input        dimensions of input array
*ARR        input        2-d image array
*NMASX,Yinput        dimensions of mask array
*MASK        input        2-d mask array
*RES        output        2-d result of correlation
*This routine is intended for image processing when the mask is much
*smaller than the original array. The "centre" of the mask is assumed to be
*(1+NMAS(1)/2,1+NMAS(2)/2) using integer arithmetic. This centre is not
*allowed to go beyond the bounds of the original image thereby keeping the
*dimensions of the result the same as the input array.
*When mask elements are correlating with values outside the input bounds
*the missing image elements are assumed to be zero.
*-Author Dick Willingale 1988-Mar-31
        INTEGER NMI,NMJ,J,I,JJ,K,II,L
        NMI=1+NMASX/2
        NMJ=1+NMASY/2
        DO J=1,NARY
        DO I=1,NARX
                RES(I,J)=0.0
                DO L=1,NMASY
                        JJ=J+L-NMJ
                        IF(JJ.GT.0.AND.JJ.LE.NARY) THEN
                                DO K=1,NMASX
                                        II=I+K-NMI
                                        IF(II.GT.0.AND.II.LE.NARX) THEN
                                                RES(I,J)=RES(I,J)+
     +                                          ARR(II,JJ)*MASK(K,L)
                                        ENDIF
                                ENDDO
                        ENDIF
                ENDDO
        ENDDO
        ENDDO
        END
