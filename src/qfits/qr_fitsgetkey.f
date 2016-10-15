*+QR_FITSGETKEY        Get FITS keyword 
        SUBROUTINE QR_FITSGETKEY(IKEY,KEY,KI,SVAL,SI,
     +        JVAL,DVAL,LVAL,KTYPE,COMMENT,CI)
        IMPLICIT NONE
        INTEGER IKEY,SI,KI,JVAL,LVAL,KTYPE,CI
        CHARACTER SVAL*(68),KEY*(8),COMMENT*(80)
        DOUBLE PRECISION DVAL
*IKEY        input        keyword index
Cf2py  intent(in) ikey
*KEY        output        keyword string returned
Cf2py  intent(out) key
*KI        output        number of characters set in keyword string
Cf2py  intent(out) ki
*SVAL        output        string value
Cf2py  intent(out) sval
*SI        output        number of characters set in string
Cf2py  intent(out) si
*JVAL        output        integer value
Cf2py  intent(out) jval
*LVAL        output        logical value
Cf2py  intent(out) lval
*DVAL        output        real value
Cf2py  intent(out) dval
*KTYPE        output        type of value returned
*                1 integer in JVAL
*                2 real in DVAL
*                3 logical in JVAL
*                4 string in SVAL
Cf2py  intent(out) ktype
*COMMENToutput        comment string
Cf2py  intent(out) comment
*CI        output        number of characters in comment string
Cf2py  intent(out) ci
*-Author: Dick Willingale 2012-Dec-27
        INCLUDE 'QR_COM'
        CHARACTER DTYPE*(1)
        LOGICAL LV
        INTEGER LEN_TRIM
        EXTERNAL LEN_TRIM
C
        IF(ISTAT.NE.0) RETURN
C
        IF(IFITS.EQ.0) THEN
                WRITE(*,*) 'QR_FITSGETKEY error - file closed'
                ISTAT=1
                RETURN
        ENDIF
C
        CALL FTGKYN(IFITS,IKEY,KEY,SVAL,COMMENT,ISTAT)
        KI=LEN_TRIM(KEY)
        SI=LEN_TRIM(SVAL)
        CI=LEN_TRIM(COMMENT)
        IF(SI.GT.0) THEN
                CALL FTDTYP(SVAL(1:SI),DTYPE,ISTAT)
        ELSE
                DTYPE='C'
        ENDIF
        IF(ISTAT.EQ.204) THEN
C Trap value string blank
                ISTAT=0
                DTYPE='C'
        ENDIF
        IF(DTYPE.EQ.'L') THEN
                KTYPE=3
                CALL FTGKYL(IFITS,KEY(1:KI),LV,COMMENT,ISTAT)
                IF(LV) THEN
                        LVAL=1
                ELSE
                        LVAL=0
                ENDIF
                CI=LEN_TRIM(COMMENT)
        ELSEIF(DTYPE.EQ.'I') THEN
                KTYPE=1
                CALL FTGKYJ(IFITS,KEY(1:KI),JVAL,COMMENT,ISTAT)
                CI=LEN_TRIM(COMMENT)
        ELSEIF(DTYPE.EQ.'F') THEN
                KTYPE=2
                CALL FTGKYD(IFITS,KEY(1:KI),DVAL,COMMENT,ISTAT)
                CI=LEN_TRIM(COMMENT)
        ELSEIF(DTYPE.EQ.'C') THEN
                KTYPE=4
                CALL FTGKYS(IFITS,KEY(1:KI),SVAL,COMMENT,ISTAT)
                IF(ISTAT.EQ.204) THEN
                        ISTAT=0
                        SI=0
                ELSE
                        SI=LEN_TRIM(SVAL)
                ENDIF
                CI=LEN_TRIM(COMMENT)
        ENDIF
        END
