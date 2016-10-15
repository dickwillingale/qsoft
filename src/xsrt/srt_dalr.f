*+SRT_DALR        Calculate distance along ray to point
        SUBROUTINE SRT_DALR(DIR,XA,P,DIS,ISTAT)
        DOUBLE PRECISION DIR(3),XA(3),P(3),DIS
        INTEGER ISTAT
*DIR    input   direction cosines of ray
*XA     input   origin of ray
*P      input   point on ray
*DIS    output  distance to point (-ve for before origin)
*ISTAT  in/out  returned status
*-Author Dick Willingale 1996-Nov-17
        DOUBLE PRECISION DR(3),CA
        IF(ISTAT.NE.0) RETURN
C        
        CALL SRT_DIDI(XA,P,DR,DIS)
        IF(DIS.GT.0.0) THEN
C Check direction cosines
                CALL SRT_VDOT(DIR,DR,CA)
                IF(CA.LT.-0.9999) THEN
                        DIS=-DIS
                ELSEIF(CA.LT.0.9999) THEN
                        ISTAT=1
                        WRITE(*,*) 'SRT_DALR error - point not on ray'
                ENDIF
        ENDIF
        END

