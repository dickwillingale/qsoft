*+QR_FITSEXISTS        See if key exists in CHDU over range of records
        SUBROUTINE QR_FITSEXISTS(IUF,IKL,IKH,KEY,IKEY,IFSTAT)
        INTEGER IUF,IKL,IKH,IKEY,IFSTAT
        CHARACTER KEY*(*)
*-Author Dick Willingale 1994-Sep-23
        CHARACTER WORD*20,VALUE*80,COM*80
C
        IF(IFSTAT.NE.0) RETURN
C
        IKEY=0
        J=IKL
        DO WHILE(IKEY.EQ.0.AND.J.LE.IKH)
                CALL FTGKYN(IUF,J,WORD,VALUE,COM,IFSTAT)
                IF(WORD.EQ.KEY) THEN
                        IKEY=J
                ENDIF
                J=J+1
        ENDDO
        END
