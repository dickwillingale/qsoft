*+QR_FITSGCX        Get FITS bit row value as character string
        SUBROUTINE QR_FITSGCX(ICOL,IROW,NBIT,SVAL)
        IMPLICIT NONE
        INTEGER ICOL,IROW,NBIT
        CHARACTER SVAL*(NBIT)
*ICOL        input        column number
Cf2py  intent(in) icol
*IROW        input        row number
Cf2py  intent(in) irow
*NBIT        input        number of bits per row
Cf2py  intent(in) nbit
*SVAL        output        row bit values returned as character string
Cf2py  intent(out) sval
*                characters 0 or 1
*-Author: Dick Willingale 2012-Dec-27
        INCLUDE 'QR_COM'
        INTEGER I,K
        LOGICAL LVAL(8)
C
        IF(ISTAT.NE.0) RETURN
C
        IF(IFITS.EQ.0) THEN
                WRITE(*,*) 'QR_FITSGCX error - file closed'
                ISTAT=1
                RETURN
        ENDIF
C
        DO I=1,NBIT,8
                CALL FTGCX(IFITS,ICOL,IROW,I,8,LVAL,ISTAT)
                DO K=I,MIN(I+7,NBIT)
                        IF(LVAL(K-I+1)) THEN
                                SVAL(K:K)='1'
                        ELSE
                                SVAL(K:K)='0'
                        ENDIF
                ENDDO
        ENDDO
        END
