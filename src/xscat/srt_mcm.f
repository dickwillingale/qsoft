*+SRT_MCM	Multiplication of complex matrices, result in place
	SUBROUTINE SRT_MCM(AR,BR,TF)
	IMPLICIT NONE
	DOUBLE COMPLEX AR(2,2),BR(2,2)
	LOGICAL TF
*AR	in/out	left hand matrix
*BR	input	right hand matrix
*TF	in/out	overflow flag (true if not overflowed)
* Calculate A*B and returns the result in A
*If overflow then A left as is
*-Author Dick Willingale 1997-Feb-17
	DOUBLE COMPLEX CR(2,2)
	INTEGER K,J,I
C Check overflow from element (1,1)
	TF=TF.AND.(ABS(AR(1,1)).LT.1.E15.AND.ABS(BR(1,1)).LT.1.E15)
	IF(TF) THEN
		DO K=1,2
			DO J=1,2
				CR(K,J)=CMPLX(0.0,0.0)
				DO I=1,2
					CR(K,J)=CR(K,J)+AR(K,I)*BR(I,J)
				ENDDO
			ENDDO
		ENDDO
		DO K=1,2
			DO J=1,2
				AR(K,J)=CR(K,J)
			ENDDO
		ENDDO
	ENDIF
	END
