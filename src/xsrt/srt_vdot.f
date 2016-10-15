*+SRT_VDOT  Take dot product of 2 vectors
	SUBROUTINE SRT_VDOT(VA,VB,DP)
	DOUBLE PRECISION VA(3),VB(3),DP
*-Author Dick Willingale 1996-Nov-16
	DP=VA(1)*VB(1)+VA(2)*VB(2)+VA(3)*VB(3)
	END
