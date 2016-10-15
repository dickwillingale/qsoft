*+QR_FITSHDU        Get FITS header dimension info
        SUBROUTINE QR_FITSHDU(IHDU,MAXD,HDUTYPE,NAXIS,NAXES,NKEYS)
        IMPLICIT NONE
        INTEGER IHDU,MAXD,HDUTYPE,NAXIS,NAXES(MAXD),NKEYS
*IHDU        input        HDU index
Cf2py  intent(in) ihdu
*MAXD        input        maximum number of dimensions
Cf2py intent(in) maxd
*HDUTYPEoutput        type of HDU
Cf2py intent(out) hdutype
*NAXIS        output        number of dimensions
Cf2py  intent(out) naxis
*NAXES        output        size of dimensions (NAXES(1)=nrows, NAXES(2)=ncols)
Cf2py  intent(out) naxes
*NKEYS        output        number of keywords
Cf2py  intent(out) nkeys
*-Author: Dick Willingale 2012-Jul-25
        INCLUDE 'QR_COM'
        INTEGER NMORE
C
        IF(ISTAT.NE.0) RETURN
C
        IF(IFITS.EQ.0) THEN
                WRITE(*,*) 'QR_FITSHDU error - file closed'
                ISTAT=1
                RETURN
        ENDIF
C
        CALL FTMAHD(IFITS,IHDU,HDUTYPE,ISTAT)
        IF(HDUTYPE.EQ.0) THEN
                CALL FTGIDM(IFITS,NAXIS,ISTAT)
                IF(NAXIS.GT.0) THEN
                        CALL FTGISZ(IFITS,MAXD,NAXES,ISTAT)
                ENDIF
        ELSEIF(HDUTYPE.EQ.1.OR.HDUTYPE.EQ.2) THEN
                NAXIS=2
                CALL FTGNRW(IFITS,NAXES(1),ISTAT)
                CALL FTGNCL(IFITS,NAXES(2),ISTAT)
        ELSE
                WRITE(*,*) 'QR_FITSHDU error unknown HDU type ',HDUTYPE
                ISTAT=1
                RETURN
        ENDIF
C Get number of keywords
        CALL FTGHSP(IFITS,NKEYS,NMORE,ISTAT)
        END
