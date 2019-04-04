*+QRT_KINGRAND  Return samples from a modified Cauchy (King profile)
        SUBROUTINE QRT_KINGRAND(ALPHA,NR,RSAM)
        IMPLICIT NONE
        INTEGER NR
        DOUBLE PRECISION ALPHA,RSAM(NR)
*ALPHA input        index of modified Cauchy distribution
*NR    input        number of samples
*RSAM  in/out       samples returned
Cf2py  intent(in) ALPHA,NR
Cf2py  intent(out) RSAM
*-Author Dick Willingale 2019-Apr-02
        INCLUDE 'QR_COM'
        INTEGER J
        DOUBLE PRECISION SRT_KINGRAND
        EXTERNAL SRT_KINGRAND
C
        IF(ISTAT.NE.0) RETURN
        DO J=1,NR
                RSAM(J)=SRT_KINGRAND(ALPHA)
        ENDDO
        END
