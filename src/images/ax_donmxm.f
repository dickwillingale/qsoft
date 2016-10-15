*+AX_DONMXM	Performs matrix multiplication on two 3x3 matrices [C]=[A].[B]
	SUBROUTINE AX_DONMXM(A,B,C)
	DOUBLE PRECISION A(3,3), B(3,3), C(3,3)
*A	Input 	Matrix
*B	Input 	Matrix
*C	Output	Returns matrix product: [C] = [A].[B]
*-
*Author Mike Watson, 1984 Jan.
	DO J = 1,3
	    DO I = 1,3
		C(I,J) = 0.D0
		DO K = 1,3
		    C(I,J) = C(I,J) + A(I,K)*B(K,J)
		END DO
	    END DO
	END DO
	END
