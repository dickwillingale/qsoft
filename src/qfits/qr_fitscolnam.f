*+QR_FITSCOLNAM   Get FITS table column name and row repeat count if variable
        SUBROUTINE QR_FITSCOLNAM(ICOL,IR,NROWS,IS,COLNAM,SI,REP)
        IMPLICIT NONE
        INTEGER ICOL,IR,NROWS,IS,SI,REP(NROWS)
        CHARACTER COLNAM*(IS)
*ICOL        input        column index
Cf2py  intent(in) icol
*IR        input        0 if variable column width
Cf2py  intent(in) ir 
*NROWS        input        number of rows
Cf2py  intent(in) nrows
*IS        input        length of name string
Cf2py  intent(in) is
*COLNAM        output        column name
Cf2py  intent(out) colnam
*SI        output        length of returned name
Cf2py  intent(out) si
*REP        output        variable repeat count for rows if IR=0
Cf2py  intent(out) rep
*-Author: Dick Willingale 2012-Dec-29
        INCLUDE 'QR_COM'
        INTEGER LEN_TRIM
        EXTERNAL LEN_TRIM
        INTEGER OFF,I
C
        IF(ISTAT.NE.0) RETURN
C
        IF(IFITS.EQ.0) THEN
                WRITE(*,*) 'QR_FITSCOLNAM error - file closed'
                ISTAT=1
                RETURN
        ENDIF
        COLNAM=TTYPE(ICOL)
        SI=LEN_TRIM(COLNAM)
        DO I=1,NROWS        
                IF(IR.EQ.0) THEN
                        CALL FTGDES(IFITS,ICOL,I,REP(I),OFF,ISTAT)
                ELSE
                        REP(I)=IR
                ENDIF
        ENDDO
        END
