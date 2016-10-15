*+SRT_VNRM Normalisation of vector to unit length (in place)
	SUBROUTINE SRT_VNRM(V,ISTAT)
	INTEGER ISTAT
	DOUBLE PRECISION V(3)
*-Author Dick Willingale 1996-Nov-16
	DOUBLE PRECISION RN
	IF(ISTAT.NE.0) RETURN
	RN=V(1)**2+V(2)**2+V(3)**2
	IF(RN.LE.0.0) THEN
		ISTAT=1
		WRITE(*,*) 'SRT_VNRM error - null vector'
	ELSE
		RN=SQRT(RN)
		V(1)=V(1)/RN
		V(2)=V(2)/RN
		V(3)=V(3)/RN
	ENDIF
	END
