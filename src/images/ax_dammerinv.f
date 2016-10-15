*+AX_DAMMERINV	Inverse Hammer (Aitoff) Projection, double precision
	SUBROUTINE AX_DAMMERINV(X,Y,AZ,EL,STATUS)
	IMPLICIT NONE
	DOUBLE PRECISION X,Y,AZ,EL
	INTEGER STATUS
*X	input	projected azimuth range -PI to PI
*Y	input	projected elevation range -PI/2 to PI/2
*AZ	output	azimuth radians
*EL	output	elevation radians
*STATUS	in/out	returned status, 0 if OK, 1 if out of range
*-Author Dick Willingale 1991-Mar-4
	DOUBLE PRECISION PI,PIBY2,XN,YN,R,RR,A,B,P,ZZ,Q,HAZ
C
	PIBY2=ASIN(1.D0)
	PI=PIBY2*2.D0
C Normalize to range -1 to 1
	XN=X/PI
	IF(XN.GT.1.0) THEN
		XN=XN-2.0D0
	ELSEIF(XN.LT.-1.0) THEN
		XN=XN+2.0D0
	ENDIF
	YN=Y/PIBY2
C
	ZZ=XN*XN+YN*YN
C Check in range
	IF(ZZ.GT.1.0D0) THEN
		STATUS=1
	ELSE
		STATUS=0
		IF(YN.NE.0.0D0) THEN
			R = XN/YN
			RR = R*R
			A = (1.0D0-ZZ)**2 +RR
			B = 1.0D0 + RR
			P = SQRT(A/B)
			EL = ACOS(P)*SIGN(1.0D0,YN)
			IF(P.NE.0.0D0) THEN
				Q=(1.0D0-ZZ)/P
			ELSE
				Q=0.0D0
			ENDIF
			HAZ = ACOS(Q)*SIGN(1.0D0,XN)
			AZ = 2.0D0*HAZ
		ELSE
			EL = 0.0D0
			Q = (1.0D0-ZZ)
			HAZ = ACOS(Q)*SIGN(1.0D0,XN)
			AZ = 2.0D0*HAZ
		ENDIF
	ENDIF
	END
