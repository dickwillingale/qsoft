*+QR_FITSCOLUNIT  Get FITS table column units keyword if it exists
        SUBROUTINE QR_FITSCOLUNIT(ICOL,IU,TUN,UI)
        IMPLICIT NONE
        INTEGER ICOL,IU,UI
        CHARACTER TUN*(IU)
*ICOL        input        column index
Cf2py  intent(in) icol
*IU        input        length of input units string
Cf2py  intent(in) iu
*TUN         output       TUNIT keyword value
Cf2py  intent(out) tun
*UI        output        length of returned units string
Cf2py  intent(out) ui
*-Author: Dick Willingale 2015-Apr-10
        INCLUDE 'QR_COM'
        INTEGER LEN_TRIM
        EXTERNAL LEN_TRIM
        INTEGER IW
        CHARACTER*(10) KWORD
        CHARACTER*(255) COM
C
        IF(ISTAT.NE.0) RETURN
C
        IF(IFITS.EQ.0) THEN
            WRITE(*,*) 'QR_FITSCOLUNIT error - file closed'
            ISTAT=1
            RETURN
        ENDIF
C Get TUNIT for column if exists
        IF(ICOL.LT.10) THEN
                IW=6
                WRITE(KWORD(6:IW),'(I1)') ICOL
        ELSEIF(ICOL.LT.100) THEN
                IW=7
                WRITE(KWORD(6:IW),'(I2)') ICOL
        ELSEIF(ICOL.LT.1000) THEN
                IW=8
                WRITE(KWORD(6:IW),'(I3)') ICOL
        ELSEIF(ICOL.LT.10000) THEN
                IW=9
                WRITE(KWORD(6:IW),'(I4)') ICOL
        ELSE
                IW=10
                WRITE(KWORD(6:IW),'(I5)') ICOL
        ENDIF
        KWORD(1:5)='TUNIT'
        CALL FTGKEY(IFITS,KWORD(1:IW),TUN,COM,ISTAT)
        IF(ISTAT.NE.0) THEN
                ISTAT=0 
                TUN=''
                UI=0
        ELSE
                UI=LEN_TRIM(TUN)
        ENDIF
        END
