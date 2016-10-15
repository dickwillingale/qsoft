*+QR_FITSTYPES        Get FITS header data types
        SUBROUTINE QR_FITSTYPES(HDUTYPE,NCOLS,CTYPE,REP)
        IMPLICIT NONE
        INTEGER HDUTYPE,NCOLS,CTYPE(NCOLS)
        INTEGER REP(NCOLS)
*HDUTYPEinput        type of HDU
Cf2py  intent(in) hdytype
*NCOLS        input        number of columns (1 if HDUTYPE=0)
Cf2py  intent(in) ncols
*CTYPE        output        column types (single type if HDUTYPE=0)
Cf2py  intent(out) ctype
*REP        output        column repeat count (0 if variable)
Cf2py  intent(out) rep
* ctype= 1 integer, 2 number (real), 3 logical, 4 string
* ctype= 5 complex, 6 double complex, 7 byte, 8 bit
*-Author: Dick Willingale 2012-Jul-25
        INCLUDE 'QR_COM'
        INTEGER BITPIX,NAXES(2)
        INTEGER J,LFORM,LUNIT,LVAL,LCOM
        PARAMETER (LVAL=400,LCOM=80,LFORM=8,LUNIT=8)
        INTEGER ROWLEN,WIDTH
        INTEGER TBCOL(MAXL),VARIDAT
        CHARACTER TUNIT(MAXL)*(LUNIT),EXTNAME*(FITS_LNAM)
        CHARACTER TFORM(MAXL)*(LFORM)
        INTEGER CTYP
        INTEGER LEN_TRIM
        EXTERNAL LEN_TRIM
C
        IF(ISTAT.NE.0) RETURN
C
        IF(IFITS.EQ.0) THEN
                WRITE(*,*) 'QR_FITSHDU error - file closed'
                ISTAT=1
                RETURN
        ENDIF
        IF(HDUTYPE.EQ.0) THEN
                CALL FTGIET(IFITS,BITPIX,ISTAT)
        ELSEIF(HDUTYPE.EQ.1) THEN
                CALL FTGHTB(IFITS,MAXL,ROWLEN,NAXES(1),NAXES(2),TTYPE,
     +                TBCOL,TFORM,TUNIT,EXTNAME,ISTAT)
        ELSEIF(HDUTYPE.EQ.2) THEN
                CALL FTGHBN(IFITS,MAXL,NAXES(1),NAXES(2),TTYPE,
     +                TFORM,TUNIT,EXTNAME,VARIDAT,ISTAT)
        ELSE
                WRITE(*,*) 'QR_FITSHDU error unknown HDU type ',HDUTYPE
                ISTAT=1
                RETURN
        ENDIF
C Decode data/column types
        IF(HDUTYPE.EQ.0) THEN
                IF(BITPIX.EQ.8) THEN
                        CTYPE(1)=1
                ELSEIF(BITPIX.EQ.16) THEN
                        CTYPE(1)=1
                ELSEIF(BITPIX.EQ.20) THEN
                        CTYPE(1)=1
                ELSEIF(BITPIX.EQ.32) THEN
                        CTYPE(1)=1
                ELSEIF(BITPIX.EQ.-32) THEN
                        CTYPE(1)=2
                ELSEIF(BITPIX.EQ.40) THEN
                        CTYPE(1)=2
                ELSEIF(BITPIX.EQ.-64) THEN
                        CTYPE(1)=2
                ENDIF
        ELSEIF(HDUTYPE.EQ.1) THEN
                DO J=1,NCOLS
                        IF(TFORM(J)(1:1).EQ.'I') THEN
C integer
                                CTYPE(J)=1
                        ELSEIF(TFORM(J)(1:1).EQ.'F'.OR.
     +                        TFORM(J)(1:1).EQ.'E'.OR.
     +                        TFORM(J)(1:1).EQ.'D') THEN
C number
                                CTYPE(J)=2
                        ELSEIF(TFORM(J)(1:1).EQ.'A') THEN
C string
                                CTYPE(J)=4
                        ENDIF
                        REP(J)=1
                ENDDO
        ELSE
                DO J=1,NCOLS
                        CALL FTBNFM(TFORM(J),CTYP,REP(J),WIDTH,ISTAT)
C                        write(*,*) j,tform(j),ctyp,rep(j),width
                        IF(CTYP.LT.0) THEN
C trap variable width column
                                CTYP=-CTYP
                                REP(J)=0
                        ENDIF
                        IF(CTYP.EQ.1) THEN
C bit, X --> integer
                                CTYPE(J)=8
                        ELSEIF(CTYP.EQ.11) THEN
C byte, B --> raw byte
                                CTYPE(J)=7
                        ELSEIF(CTYP.EQ.14) THEN
C logical, L --> integer
                                CTYPE(J)=3
                        ELSEIF(CTYP.EQ.16) THEN
C character, A --> string
                                CTYPE(J)=4
                        ELSEIF(CTYP.EQ.21) THEN
C short integer, I --> integer
                                CTYPE(J)=1
                        ELSEIF(CTYP.EQ.41) THEN
C integer, J --> integer
                                CTYPE(J)=1
                        ELSEIF(CTYP.EQ.42) THEN
C real, E --> number
                                CTYPE(J)=2
                        ELSEIF(CTYP.EQ.81) THEN
C integer, K --> integer
                                CTYPE(J)=1
                        ELSEIF(CTYP.EQ.82) THEN
C double precision, D --> number
                                CTYPE(J)=2
                        ELSEIF(CTYP.EQ.83) THEN
C complex --> complex
                                CTYPE(J)=5
                        ELSEIF(CTYP.EQ.163) THEN
C double complex --> complex
                                CTYPE(J)=6
                        ENDIF
                ENDDO
        ENDIF
        END
