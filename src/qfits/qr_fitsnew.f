*+QR_FITSNEW        Open new FITS data file
        SUBROUTINE QR_FITSNEW(IN,FNAME)
        IMPLICIT NONE
        INTEGER IN
        CHARACTER FNAME*(IN)
*IN        input        length of name
Cf2py  intent(in) in
*FNAME        input        file name
Cf2py  intent(in) fname
*-Author: Dick Willingale 2013-Feb-09
        INCLUDE 'QR_COM'
        LOGICAL LVAL
C
        IF(ISTAT.NE.0) RETURN
	write(*,*) 'start qr_fitsnew',FNAME(1:in)
C
        IF(IFITS.GT.0) THEN
                CALL FTCLOS(IFITS,ISTAT)
                IFITS=0
        ENDIF
C Check if file already exists
        INQUIRE(FILE=FNAME(1:IN),EXIST=LVAL)
        IF(LVAL) THEN
                WRITE(*,*) FNAME(1:IN),' already exists'
CALL SYS_UNLINK(FNAME(1:IN),ISTAT)
        ENDIF
C Open new file
        CALL SYS_GETLUN(IFITS,ISTAT)
        CALL FTINIT(IFITS,FNAME,0,ISTAT)
        IF(ISTAT.NE.0) THEN
                WRITE(*,*) 'QR_FITSNEW failed to open ',FNAME
                IFITS=0
                RETURN
        ENDIF
C Set extension 
        IEXT=0
        END
