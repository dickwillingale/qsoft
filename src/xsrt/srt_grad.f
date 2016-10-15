*+SRT_GRAD	Calculate gradient
	SUBROUTINE SRT_GRAD(NP,XP,YP,NSMOOTH,GR)
	IMPLICIT NONE
	INTEGER NP,NSMOOTH
	DOUBLE PRECISION XP(NP),YP(NP),GR(NP)
*NP	input	number of points
*XP	input	x positions
*YP	input	y positions
*NSMOOTHinput	range over which to interpolate for each gradient estimate
*GR	output	gradient vector
*-Author Dick Willingale 1995-Dec-11
	DOUBLE PRECISION SX,SY,SXX,SXY,RN,DEN
	INTEGER NW,K,J
C
	NW=MAX(NSMOOTH/2,3)
	DO J=1,NP
C Find gradient by least squares fit over a range
		SX=0.0
		SY=0.0
		SXX=0.0
		SXY=0.0
		RN=0.0
		DO K=J-NW,J+NW
			IF(K.GE.1.AND.K.LE.NP) THEN
				RN=RN+1.0
				SX=SX+XP(K)
				SY=SY+YP(K)
				SXX=SXX+XP(K)*XP(K)
				SXY=SXY+XP(K)*YP(K)
			ENDIF
		ENDDO
		DEN=RN*SXX-SX*SX
		IF(DEN.NE.0.0) THEN
			GR(J)=(RN*SXY-SX*SY)/DEN
		ELSE
			GR(J)=0.0
		ENDIF
	ENDDO
	END
