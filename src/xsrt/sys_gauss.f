*+SYS_GAUSS	Internal routine to perform vector initialisation
	SUBROUTINE SYS_GAUSS(NEL,V,MEAN,SIGMA,ISTAT)
	INTEGER NEL,ISTAT
	DOUBLE PRECISION V(NEL),MEAN,SIGMA
*NEL	input	number of vector elements
*V	output	vector for result
*MEAN	input	mean of Gaussian distribution
*SIGMA	input	width of Gaussian distribution
*ISTAT	output	returned status
*-Author Dick Willingale 1990-Nov-9
* Use algorithm from Numerical Recipes page 279-280
C	INCLUDE 'Q_COMMON'
	DOUBLE PRECISION V1,V2,RSQ,FAC,SYS_DRAND
C
	IF(ISTAT.NE.0) RETURN
C Generate sequence
	DO J=1,NEL,2
		RSQ=0.0
		DO WHILE(RSQ.GE.1.D0.OR.RSQ.EQ.0.0)
			V1=SYS_DRAND()*2.D0-1.D0
			V2=SYS_DRAND()*2.D0-1.D0
			RSQ=V1**2+V2**2
		ENDDO
		FAC=SQRT(-2.D0*LOG(RSQ)/RSQ)
		V(J)=FAC*V1*SIGMA+MEAN
		IF(J.LT.NEL) V(J+1)=FAC*V2*SIGMA+MEAN
	ENDDO
	END
