*+SRT_BDIPOLE	B field of a magnetic dipole
	SUBROUTINE SRT_BDIPOLE(D,P,GHAT,F,B,R,ISTAT)
	DOUBLE PRECISION D,P(3),GHAT(3),F(3),B(3),R
	INTEGER ISTAT
*D	input	dipole moment (Gauss cm3)
*P	input	position of dipole (cm)
*GHAT	input	direction cosines of dipole moment
*F	input	position for field (cm)
*B	output	magnetic field at position F (Gauss)
*R	output	distance from dipole
*ISTAT	in/out	returned status
*-Author Dick Willingale 1994-Sep-14
	DOUBLE PRECISION RHAT(3),AHAT(3),THAT(3),CTH,STH
C
	if(ISTAT.ne.0) return
C Calculate distance R from dipole position to field position
	RHAT(1)=F(1)-P(1)
	RHAT(2)=F(2)-P(2)
	RHAT(3)=F(3)-P(3)
	R=SQRT(RHAT(1)**2+RHAT(2)**2+RHAT(3)**2)
C Calculate radial unit vector RHAT
	IF(R.GT.0.1D0) THEN
		RHAT(1)=RHAT(1)/R
		RHAT(2)=RHAT(2)/R
		RHAT(3)=RHAT(3)/R
C Calculate elevation sine and cosine of position wrt dipole direction
		CALL SRT_VDOT(RHAT,GHAT,CTH)
		STH=SQRT(ABS(1.D0-CTH*CTH))
		IF(STH.GT.0.D0) THEN
C Calculate azimuthal unit vector AHAT
			CALL SRT_VCRS(RHAT,GHAT,AHAT)
			CALL SRT_VNRM(AHAT,ISTAT)
C Calculate the tangential unit vector THAT
			CALL SRT_VCRS(RHAT,AHAT,THAT)
			CALL SRT_VNRM(THAT,ISTAT)
		ELSE
			THAT(1)=0.D0
			THAT(2)=0.D0
			THAT(3)=0.D0
		ENDIF
C Sum radial and tangential components of field
		RR=D/R**3
		B(1)=RR*(2.D0*CTH*RHAT(1)+STH*THAT(1))
		B(2)=RR*(2.D0*CTH*RHAT(2)+STH*THAT(2))
		B(3)=RR*(2.D0*CTH*RHAT(3)+STH*THAT(3))
	ELSE
C Too close so set to zero
		B(1)=0.D0
		B(2)=0.D0
		B(3)=0.D0
	ENDIF
	END
