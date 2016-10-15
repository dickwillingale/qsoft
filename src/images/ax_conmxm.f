*+AX_CONMXM	Performs matrix multiplication on two 3x3 matrices [C]=[A].[B]
	SUBROUTINE AX_CONMXM(A,B,C)
	REAL A(3,3), B(3,3), C(3,3)
*A	Input 	Matrix
*B	Input 	Matrix
*C	Output	Returns matrix product: [C] = [A].[B]
*-Author	Clive Page	1988 Feb 22 (After AX_DONMXM of Mike Watson)
	DO J = 1,3
	    DO I = 1,3
		C(I,J) = 0.0
		DO K = 1,3
		    C(I,J) = C(I,J) + A(I,K) * B(K,J)
		END DO
	    END DO
	END DO
	END
