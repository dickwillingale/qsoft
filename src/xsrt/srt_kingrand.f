*+SRT_KINGRAND	return a random number from a King profile
	FUNCTION SRT_KINGRAND(ALPHA)
	DOUBLE PRECISION SRT_KINGRAND,ALPHA
* Modified Cauchy or King profile random sample
* ALPHA    index of the modified Cauchy distribution 1/(1+x^2)^alpha
*-Author Dick Willingale 2019_Mar_31
	INTEGER IDUM
	COMMON/RSEED/IDUM
        INTEGER NA,I,N2,J
        DOUBLE PRECISION STEP,S,HW
	PARAMETER (STEP=0.01)
	PARAMETER (NA=4001)
        DOUBLE PRECISION RIND,UP(NA),X(NA),Y(NA),XX
        COMMON/KING/RIND,UP
	DOUBLE PRECISION SYS_DRAN0
	EXTERNAL SYS_DRAN0
	SYS_DRAND=SYS_DRAN0(IDUM)
	IF(RIND.NE.ALPHA) THEN
		RIND=ALPHA
		HW=(NA-1)*STEP*0.5
		Y(1)=0.0
		X(1)=-HW
	        DO I=2,NA
			X(I)=X(I-1)+STEP
			Y(I)=1.0/((1.0+X(I)**2)**RIND)
			Y(I)=Y(I-1)+Y(I)
		ENDDO
C Normalise to maximum probability of 1
		S=Y(NA)
	        DO I=1,NA
			Y(I)=Y(I)/S
		ENDDO
C Interpolate to set up lookup table of deviates
		J=2
		DO I=1,NA
			YY=DBLE(I-1)/DBLE(NA-1)
			DO WHILE(Y(J).LT.YY.AND.J.LT.NA)
				J=J+1
			ENDDO
			UP(I)=(YY-Y(J-1))/(Y(J)-Y(J-1))*(X(J)-X(J-1))
			UP(I)=UP(I)+X(J-1)
		ENDDO
	ENDIF
	S=SYS_DRAN0(IDUM)*(NA-1)
	I=INT(S)+1
	SRT_KINGRAND=UP(I)+(UP(I+1)-UP(I))*(S-(I-1))
	END
