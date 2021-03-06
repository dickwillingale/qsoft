*+QR_FITSPARRJ        put integer array as primary or image array
        SUBROUTINE QR_FITSPARRJ(NDIMS,NELS,NEL,ARR)
        IMPLICIT NONE
        INTEGER NDIMS,NELS(NDIMS),NEL,ARR(NEL)
*NDIMS        input        number of dimensions
Cf2py  intent(in) ndims
*NELS        input        dimensions
Cf2py  intent(in) nels
*NEL        input        total number of elements
Cf2py  intent(in) nel
*ARR        input        data array
Cf2py  intent(in) arr
*-Author: Dick Willingale 2013-Feb-09
        INCLUDE 'QR_COM'
        INTEGER IHTYPE,BITPIX
C
        IF(ISTAT.NE.0) RETURN
C
        IF(IFITS.EQ.0) THEN
                WRITE(*,*) 'QR_FITSPARRI - fits file not open'
                ISTAT=1
                RETURN
        ENDIF
C Create new extension if not primary
        IF(IEXT.GT.0) THEN
                CALL FTCRHD(IFITS,ISTAT)
        ENDIF
C Set extension
        IEXT=IEXT+1
        BITPIX=32
C Set header
        CALL FTPHPR(IFITS,.TRUE.,BITPIX,NDIMS,NELS,0,1,.TRUE.,ISTAT)
C Move to header
        CALL FTMAHD(IFITS,IEXT,IHTYPE,ISTAT)
C Switch off scaling
        CALL FTPSCL(IFITS,1.D0,0.D0,ISTAT)
C Put array into fits
        CALL FTPPRJ(IFITS,0,1,NEL,ARR,ISTAT)
        END
