*+QR_FITSOPEN        Open FITS data file
        SUBROUTINE QR_FITSOPEN(IN,NAM,NRW,NHDU)
        IMPLICIT NONE
        INTEGER IN,NHDU,NRW
        CHARACTER NAM*(IN)
*IN        input        length of name
Cf2py  intent(in) in
*NAM        input        file name
Cf2py  intent(in) nam
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
C        write(*,*) 'qr_fitsopen',nam(1:in),nrw,ifits

        CALL FTOPEN(IFITS,NAM(1:IN),NRW,IBLOCK,ISTAT)
C        write(*,*) 'qr_fitsopen',nam(1:in),nrw,ifits
        IF(ISTAT.NE.0) THEN
              WRITE(*,*) 'QR_FITSOPEN failed ',NAM(1:IN)
              RETURN
C        ELSE
C             write(*,*) 'qr_fitsopen done'
        ENDIF
        CALL FTTHDU(IFITS,NHDU,ISTAT)
        END
