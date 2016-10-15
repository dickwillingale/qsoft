*+QR_FITSOPEN        Open FITS data file
        SUBROUTINE QR_FITSOPEN(IN,NAME,NRW,NHDU)
        IMPLICIT NONE
        INTEGER IN,NHDU,NRW
        CHARACTER NAME*(IN)
*IN        input        length of name
Cf2py  intent(in) in
*NAME        input        file name
Cf2py  intent(in) name
*NRW        input        readonly 0, read/write 1
Cf2py  intent(in) nrw
*NDHU        output        number of HDU in file
Cf2py  intent(out) nhdu
*-Author: Dick Willingale 2012-Jul-25
        INCLUDE 'QR_COM'
        INTEGER IBLOCK
C
        IF(ISTAT.NE.0) RETURN
C
        IF(IFITS.GT.0) THEN
                CALL FTCLOS(IFITS,ISTAT)
                IFITS=0
        ENDIF
C
        CALL SYS_GETLUN(IFITS,ISTAT)
        CALL FTOPEN(IFITS,NAME,NRW,IBLOCK,ISTAT)
        IF(ISTAT.NE.0) THEN
                WRITE(*,*) 'QR_FITSOPEN failed ',NAME
                RETURN
        ENDIF
        CALL FTTHDU(IFITS,NHDU,ISTAT)
        END
