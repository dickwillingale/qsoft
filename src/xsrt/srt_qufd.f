*+SRT_QUFD      Evaluate derivitive of quadratic function
	DOUBLE PRECISION FUNCTION SRT_QUFD(P,X)
	DOUBLE PRECISION P(3),X
*-Author Dick Willingale 1996-Nov-20
	SRT_QUFD=2.0*P(1)*X+P(2)
	END
