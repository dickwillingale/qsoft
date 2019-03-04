*+QRT_INIT        Initialise/reset the QRT environment
        SUBROUTINE QRT_INIT()
        IMPLICIT NONE
*-Author Dick Willingale 2012-May-04
        INCLUDE 'SRT_COM'
        INTEGER J
C
        CALL QR_INIT()
C
        NPAR=0
        NSUR=0
        IDSAM=0
        ITRA=0
        NRAYS=0
        ISRC(1)=0
        ISRC(2)=0
        ISRC(3)=0
        NDET=0
        NPOS=0
        IDET=0
        IDEBUG=0
        DO J=1,MAXST
                ISQP(1,J)=0
        ENDDO
        DO J=1,MAXDF
                IDFM(1,J)=0
        ENDDO
C        write(*,*) 'qpy initialises/reset'
        END
