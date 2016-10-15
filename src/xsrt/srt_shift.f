*+SRT_SHIFT	Linear shift surface element
	SUBROUTINE SRT_SHIFT(IS,VSH,ISTAT)
	IMPLICIT NONE
	INTEGER IS,ISTAT
	DOUBLE PRECISION VSH(3)
*IS	input	surface element number
*VSH	input	vector shift of position
*ISTAT	in/out	returned status
*-Author Dick Willingale 2003-Sep-11
	INCLUDE 'SRT_COM'
	INTEGER I
C
	IF(ISTAT.NE.0) RETURN
C
	IF(IS.LT.1.OR.IS.GT.NSUR) THEN
		WRITE(*,*) 'SRT_SHIFT error - element ',IS,' not set'
		ISTAT=1
		RETURN
	ENDIF
C shift position of element
	I=IPAR(IS)+6
	PAR(I)=PAR(I)+VSH(1)
	PAR(I+1)=PAR(I+1)+VSH(2)
	PAR(I+2)=PAR(I+2)+VSH(3)
	END
