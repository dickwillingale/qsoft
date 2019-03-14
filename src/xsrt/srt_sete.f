*+SRT_SETE	Set reflection efficiency for surface quality
	SUBROUTINE SRT_SETE(IS,EFF,ISTAT)
	IMPLICIT NONE
	INTEGER IS,ISTAT
	DOUBLE PRECISION EFF
*IS	input	surface index
*EFF	input	efficiency wrt theory for reflection
*ISTAT	in/out	returned status
*-Author Dick Willingale 2019-Mar-14
	INCLUDE 'SRT_COM'
        INTEGER NPI,IT
C
	IF(ISTAT.NE.0) RETURN
C Check surface index
	IF(IS.LT.0.OR.IS.GT.MAXST) THEN
		WRITE(*,*) 'SRT_SETT error - surface index out of range'
		ISTAT=1
		RETURN
	ENDIF
C Get type and index
	IT=ISQP(1,IS)
	NPI=ISQP(2,IS)
C Set parameter
        IF(IT.EQ.1.OR.IT.EQ.2) THEN
                PAR(NPI+4)=EFF
	ELSE
		WRITE(*,*) 'SRT_SETE error - surface type not 1 or 2'
		ISTAT=1
		RETURN
	ENDIF
	END
