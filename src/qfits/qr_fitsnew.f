*+QR_FITSNEW        Open new FITS data file
        SUBROUTINE QR_FITSNEW(II,FNAME)
        IMPLICIT NONE
        INTEGER II
        CHARACTER FNAME*(*)
Cf2py  intent(inout) fname
Cf2py  intent(in) in
*II        input        length of name
*FNAME        input        file name
*-Author: Dick Willingale 2013-Feb-09
        INCLUDE 'QR_COM'
        LOGICAL LVAL
C
        IF(ISTAT.NE.0) RETURN
C
        IF(IFITS.GT.0) THEN
                CALL FTCLOS(IFITS,ISTAT)
                IFITS=0
        ENDIF
C Check if file already exists
        INQUIRE(FILE=FNAME(1:II),EXIST=LVAL)
        IF(LVAL) THEN
                WRITE(*,*) FNAME(1:II),' already exists, overwriting'
                CALL SYS_UNLINK(FNAME(1:II),ISTAT)
        ENDIF
C Open new file
        CALL SYS_GETLUN(IFITS,ISTAT)
        CALL FTINIT(IFITS,FNAME(1:II),0,ISTAT)
        IF(ISTAT.NE.0) THEN
                WRITE(*,*) 'QR_FITSNEW failed to open ',FNAME
                IFITS=0
                RETURN
        ENDIF
C Set extension 
        IEXT=0
        END
