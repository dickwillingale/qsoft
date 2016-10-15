*+QR_FITSGCL        Get FITS column value as logical
        SUBROUTINE QR_FITSGCL(ICOL,IROW,IELE,IVAL,IUDEF)
        IMPLICIT NONE
        INTEGER ICOL,IROW,IELE,IVAL,IUDEF
*ICOL        input        column number
Cf2py  intent(in) icol
*IROW        input        row number
Cf2py  intent(in) irow
*IELE        input        element number
Cf2py  intent(in) iele
*IVAL        output        logical value returned as integer (for R)
Cf2py  intent(out) ival
*IUDEF        output        1 if not defined
Cf2py  intent(out) iudef
*-Author: Dick Willingale 2012-Dec-27
        INCLUDE 'QR_COM'
        LOGICAL LVAL,LFL,LDF
C
        IF(ISTAT.NE.0) RETURN
C
        IF(IFITS.EQ.0) THEN
                WRITE(*,*) 'QR_FITSGCL error - file closed'
                ISTAT=1
                RETURN
        ENDIF
C
        CALL FTGCL(IFITS,ICOL,IROW,IELE,1,LVAL,ISTAT)
        CALL FTGCFL(IFITS,ICOL,IROW,IELE,1,LVAL,LFL,LDF,ISTAT)
        IF(LVAL) THEN
                IVAL=1
        ELSE
                IVAL=0
        ENDIF
        IF(LDF) THEN
                IUDEF=1
        ELSE
                IUDEF=0
        ENDIF
        END
