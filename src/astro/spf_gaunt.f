*+SPF_GAUNT	Born approximation to Gaunt factor
	FUNCTION SPF_GAUNT(E,T)
	REAL SPF_GAUNT,E,T
*E		Energy keV
*T		Temperature keV
*SPF_GAUNT	Gaunt factor at energy E
*-Author Dick Willingale 1986-Sep-4
C From Gordon Stewart's code May 1984
	X=(E/T)*0.5
	IF(X.LT.2.) THEN
		T2=(X/3.75)**2
		AIO=(((((.0045813*T2+.0360768)*T2+.2659732)*T2+1.2067492)*T2+
     +		3.0899424)*T2+3.5156229)*T2+1.0
		Y=X*0.5
		Z=ALOG(Y)
		Y2=Y**2
		BO=-AIO*Z+(((((.0000074*Y2+.0001075)*Y2+.00262698)*Y2+
     +		.0348859)*Y2+.23069756)*Y2+.4227842)*Y2-.57721566
		SPF_GAUNT=.551329*EXP(X)*BO
	ELSE
		Y=2./X
		BO=(((((.00053208*Y-.0025154)*Y+.00587872)*Y-.01062446)*Y+
     +		.02189568)*Y-.07832358)*Y+1.25331414
		FAC=BO/SQRT(X)
		SPF_GAUNT=.551329*FAC
	ENDIF
	END
