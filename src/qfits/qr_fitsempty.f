*+QR_FITSEMPTY        set primary or image array as empty
        SUBROUTINE QR_FITSEMPTY()
        IMPLICIT NONE
*-Author: Dick Willingale 2013-Feb-09
        INCLUDE 'QR_COM'
        INTEGER IHTYPE
C
        IF(ISTAT.NE.0) RETURN
C
        IF(IFITS.EQ.0) THEN
                WRITE(*,*) 'QR_FITSEMPTY - fits file not open'
                ISTAT=1
                RETURN
        ENDIF
C Create new extension if not primary
        IF(IEXT.GT.0) THEN
                CALL FTCRHD(IFITS,ISTAT)
        ENDIF
C Set extension
        IEXT=IEXT+1
C Set header
        CALL FTPHPR(IFITS,.TRUE.,8,0,0,0,1,.TRUE.,ISTAT)
C Move to header
        CALL FTMAHD(IFITS,IEXT,IHTYPE,ISTAT)
        END
