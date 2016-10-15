*+SRT_NSUR	Get next free surface number
	SUBROUTINE SRT_NSUR(KSUR,ISTAT)
	INTEGER KSUR,ISTAT
*KSUR	output	next free surface number
*ISTAT	in/out	returned status
*-Author Dick Willingale 1996-Dec-16
	INCLUDE 'SRT_COM'
C
	IF(ISTAT.NE.0) RETURN
C New surface
	IF(NSUR.EQ.MAXSUR) THEN
		WRITE(*,*) 'SRT_SETF error - surface list full'
		ISTAT=1
		RETURN
	ENDIF
	KSUR=NSUR+1
	END	
