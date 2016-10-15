*+QRT_LIST List all current QRT parameters
        SUBROUTINE QRT_LIST()
*-Author Dick Willingale 2012-May-04
        IMPLICIT NONE
        INTEGER IU
        INCLUDE 'QR_COM'
C
        IF(ISTAT.NE.0) RETURN
C
        IU=6
        CALL SRT_LIST(IU,ISTAT)
        END
