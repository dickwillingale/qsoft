*+SRT_DIDI      Calculate direction and distance between positions
        SUBROUTINE SRT_DIDI(V1,V2,DIR,DIS)
        DOUBLE PRECISION V1(3),V2(3),DIR(3),DIS
*V1     input   position vector
*V2     input   position vector
*DIR    output  direction cosines of V1 to V2 (1,0,0 if V1=V2)
*DIS    output  distance from V1 to V2
*ISTAT  in/out  returned status
*-Author Dick Willingale 1996-Nov-17
        DIS=0.0
        DO J=1,3
                DIR(J)=V2(J)-V1(J)
                DIS=DIS+DIR(J)*DIR(J)
        ENDDO
        IF(DIS.GT.0.0) THEN
                DIS=SQRT(DIS)
                DO J=1,3
                        DIR(J)=DIR(J)/DIS
                ENDDO
        ELSE
                DIR(1)=1.0
                DIR(2)=0.0
                DIR(3)=0.0
        ENDIF
        END

