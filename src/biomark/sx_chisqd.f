*+SX_CHISQD	Computes chi-squared deviate for given probability and freedom.
	REAL FUNCTION SX_CHISQD(PROB,NFREE)
	REAL PROB
	INTEGER NFREE
*PROB	input	Probability level, e.g. 0.01 for 1 %.
*NFREE	input	Number of degrees of freedom, range 1 upwards.
*SX_CHISQD output	Returns deviate corresponding to the given upper
*			tail area (i.e. probability).
*Accuracy: better than 3 digits over entire range.
*-
*Author	Clive Page	1983-JULY-5
*Standard Fortran-77 version 1984 Sept 25, CGP.
*Algorithm of R.B.Goldstein, Comm. Ass. Comp. Mach. 16,483 (1973)

	PARAMETER (C1=1.565326E-3, C2=1.1060438E-3, C3=-6.950356E-3,
     &   C4 =-1.323293E-2,  C5 = 2.277679E-2,  C6 =-8.986007E-3,
     &   C7 =-1.51390E-2,   C8 = 2.530010E-3,  C9 =-1.450117E-3,
     &   C10= 5.169654E-3,  C11=-1.153761E-2,  C12= 1.128186E-2,
     &   C13= 2.607083E-2,  C14=-0.2237368,    C15= 9.780499E-5,
     &   C16=-8.426812E-4,  C17= 3.125580E-3,  C18=-8.553069E-3,
     &   C19= 1.348028E-4,  C20= 0.4713941,    C21= 1.0000886)
	PARAMETER (A1=1.264616E-2, A2=-1.425296E-2, A3=1.400483E-2,
     &   A4 =-5.886090E-3, A5 =-1.091214E-2, A6 =-2.304527E-2,
     &   A7 = 3.135411E-3, A8 =-2.728484E-4, A9 =-9.9699681E-3,
     &   A10= 1.316872E-2, A11= 2.618914E-2, A12=-0.2222222,
     &   A13= 5.406674E-5, A14= 3.483789E-5, A15=-7.274761E-4,
     &   A16= 3.292181E-3, A17=-8.729713E-3, A18= 0.4714045, A19=1.0)

	IF(NFREE .EQ. 1) THEN
	    SX_CHISQD = (SX_GAUSSD(0.5 * PROB))**2

	ELSE IF(NFREE .EQ. 2) THEN
	    SX_CHISQD = -2.0 * LOG(PROB)

	ELSE
	    F1 = 1.0 / NFREE
	    T  = SX_GAUSSD(PROB)
	    F2 = SQRT(F1) * T
	    IF(NFREE .LT. (2+INT(4.0*ABS(T)))) THEN
		SX_CHISQD =
     &   (((((((C1*F2+C2)*F2+C3)*F2+C4)*F2+C5)*F2+C6)*F2+C7)*F1+
     &   ((((((C8+C9*F2)*F2+C10)*F2+C11)*F2+C12)*F2+C13)*F2+C14))*F1+
     &   (((((C15*F2+C16)*F2+C17)*F2+C18)*F2+C19)*F2+C20)*F2+C21
	    ELSE
		SX_CHISQD =
     &   (((A1+A2*F2)*F1+(((A3+A4*F2)*F2+A5)*F2+A6))*F1+
     &   (((((A7+A8*F2)*F2+A9)*F2+A10)*F2+A11)*F2+A12))*F1+
     &   (((((A13*F2+A14)*F2+A15)*F2+A16)*F2+A17)*F2*F2+A18)*F2+A19
	    END IF
	    SX_CHISQD = NFREE * SX_CHISQD**3
	END IF
	END
