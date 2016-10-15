*+SRT_QUFU	Evaluate quadratic function
	DOUBLE PRECISION FUNCTION SRT_QUFU(P,X)
	DOUBLE PRECISION P(3),X
*-Author Dick Willingale 1996-Nov-20
	SRT_QUFU=P(1)*X**2+P(2)*X+P(3)
	END
