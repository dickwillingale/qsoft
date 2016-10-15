	SUBROUTINE SRT_KBFINDMODULE(INRAY,NROFELTS,MAXLEN,
     +	RLIST,PHILIST,ELTSIZE,INRAD,OUTRAD,MODHIT)

* this file by FS october 2008


* this routine determines in which module the incoming ray falls

	IMPLICIT NONE

	INTEGER NROFELTS,MODHIT,MAXLEN
	DOUBLE PRECISION RLIST(MAXLEN),PHILIST(MAXLEN),
     +	C,INRAD,OUTRAD,ELTSIZE,INRAY(2)

	DOUBLE PRECISION
     *	RRAY,PHIRAY,RLOW,RHIGH,PHILOW,PHIHIGH,PHIGOLD,PI,ANGLE,
     +	HALFRIB,A(4),B(4),R,PHI,MODPHI,
     *	ALLLIST(MAXLEN,3),CENTRE(2),CORNER(4,2)
	INTEGER INDEX,COUNTER,NMODS,INDEXOFF,
     *	INDEXBEG,INDEXEND,SINDBEG,SINDEND,MODLIST(MAXLEN),
     *	UNITNR,READSTATUS,MODINDEX,QUAD
	LOGICAL FOPEN,REFILE

* initialize:
	PI=ASIN(1.D0)*2.D0
	PHIGOLD=2.3999632297286533
	HALFRIB=SQRT(2.D0)*ELTSIZE/2
	C=ELTSIZE/(PI/SQRT(2.D0))

	HALFRIB=HALFRIB/2.0


* first a shortlist of modules is drafted:
* the r-bounderies are calculated, the phi-bounderies have to be looped

* cartesian to polar:
	RRAY=SQRT(INRAY(1)**2+INRAY(2)**2)
	PHIRAY=ATAN2(INRAY(2),INRAY(1))
* the bounds for defining a neighbour:
	RLOW=RRAY-SQRT(2.)*ELTSIZE/2
	RHIGH=RRAY+SQRT(2.)*ELTSIZE/2
	PHILOW=PHIRAY-ATAN2(SQRT(2.)*ELTSIZE/2,RLOW)
	PHIHIGH=PHIRAY+ATAN2(SQRT(2.)*ELTSIZE/2,RLOW)
* calculate the index range from the r values:
	INDEXBEG=INT(((RLOW/C)**2)/PHIGOLD)+1
	INDEXEND=INT(((RHIGH/C)**2)/PHIGOLD)
* this way of counting starts at phi=0, but the lists start from 
* the inner radius, so an offset is necessary:
	INDEXOFF=INT((((INRAD+SQRT(2.)*ELTSIZE/2)/C)**2)/PHIGOLD)+1
* set the search index range:
	SINDBEG=MAX(INDEXBEG-INDEXOFF+1,1)
	SINDEND=INDEXEND-INDEXOFF+1
* now loop over the phi values:
	COUNTER=0
	DO INDEX=SINDBEG,SINDEND
* note the range issue: normally -pi->pi, but Mod: 0->2pi
	   IF(MOD(PHILIST(INDEX),2*PI).LE.PI)THEN
	      MODPHI=MOD(PHILIST(INDEX),2*PI)
	   ELSE
	      MODPHI=MOD(PHILIST(INDEX),2*PI)-2*PI
	   ENDIF
	   IF(PHILOW.LT.MODPHI.AND.
     +	   MODPHI.LT.PHIHIGH)THEN
	      COUNTER=COUNTER+1
	      MODLIST(COUNTER)=INDEX
	   ENDIF
	ENDDO
	NMODS=COUNTER
*	WRITE(*,*) (MODLIST(INDEX),INDEX=1,NMODS)


* secondly, find out which module of the shortlist is the right one (if any): 
	IF(NMODS.EQ.0)THEN
	   MODHIT=0
	ENDIF
	IF(NMODS.GT.1) THEN
	   DO MODINDEX=1,NMODS
* the centre of the module:
	      CENTRE(1)=RLIST(MODLIST(MODINDEX))*
     +		COS(PHILIST(MODLIST(MODINDEX)))
	      CENTRE(2)=RLIST(MODLIST(MODINDEX))*
     +		SIN(PHILIST(MODLIST(MODINDEX)))
* the angle and the quadrant of the ray position wrt the centre of the module:
	      ANGLE=ATAN2(CENTRE(2),CENTRE(1))
	      IF (ANGLE.GE.-PI/4.AND.ANGLE.LT.PI/4) QUAD=1
	      IF (ANGLE.GE.PI/4.AND.ANGLE.LT.3*PI/4) QUAD=2
	      IF (ANGLE.GE.3*PI/4.OR.ANGLE.LT.-3*PI/4) QUAD=3
	      IF (ANGLE.GE.-3*PI/4.AND.ANGLE.LT.-PI/4) QUAD=4
* the corners of the module:
* corner 1
	      R=RLIST(MODLIST(MODINDEX))-HALFRIB
	      PHI=PHILIST(MODLIST(MODINDEX))
	      CORNER(1,1)=R*COS(PHI)
	      CORNER(1,2)=R*SIN(PHI)
* corner 2
	      R=SQRT(RLIST(MODLIST(MODINDEX))**2+HALFRIB**2)
	      PHI=PHILIST(MODLIST(MODINDEX))+
     +	      ATAN2(HALFRIB,RLIST(MODLIST(MODINDEX)))
	      CORNER(2,1)=R*COS(PHI)
	      CORNER(2,2)=R*SIN(PHI)
* corner 3
	      R=RLIST(MODLIST(MODINDEX))+HALFRIB
	      PHI=PHILIST(MODLIST(MODINDEX))
	      CORNER(3,1)=R*COS(PHI)
	      CORNER(3,2)=R*SIN(PHI)
* corner 4
	      R=SQRT(RLIST(MODLIST(MODINDEX))**2+HALFRIB**2)
	      PHI=PHILIST(MODLIST(MODINDEX))-
     +	      ATAN2(HALFRIB,RLIST(MODLIST(MODINDEX)))
	      CORNER(4,1)=R*COS(PHI)
	      CORNER(4,2)=R*SIN(PHI)

* the sides of the module, described by  y=a*x+b:
	      A(1)=(CORNER(2,2)-CORNER(1,2))/(CORNER(2,1)-CORNER(1,1))
	      B(1)=CORNER(1,2)-A(1)*CORNER(1,1)
	      A(2)=(CORNER(3,2)-CORNER(2,2))/(CORNER(3,1)-CORNER(2,1))
	      B(2)=CORNER(2,2)-A(2)*CORNER(2,1)
	      A(3)=(CORNER(4,2)-CORNER(3,2))/(CORNER(4,1)-CORNER(3,1))
	      B(3)=CORNER(3,2)-A(3)*CORNER(3,1)
	      A(4)=(CORNER(1,2)-CORNER(4,2))/(CORNER(1,1)-CORNER(4,1))
	      B(4)=CORNER(4,2)-A(4)*CORNER(4,1)

* see if the ray is within the sides of the module;
* different criteria are necessary for different quadrants:
	      MODHIT=0
	      IF (QUAD.EQ.1.AND.
     +	      INRAY(2).LT.A(1)*INRAY(1)+B(1).AND.
     +	      INRAY(2).LT.A(2)*INRAY(1)+B(2).AND.
     +	      INRAY(2).GT.A(3)*INRAY(1)+B(3).AND.
     +	      INRAY(2).GT.A(4)*INRAY(1)+B(4)
     +	      )THEN
	      MODHIT=MODLIST(MODINDEX)
	      ENDIF

	      IF (QUAD.EQ.2 .AND.
     +	      INRAY(2).GT.A(1)*INRAY(1)+B(1).AND.
     +	      INRAY(2).LT.A(2)*INRAY(1)+B(2).AND.
     +	      INRAY(2).LT.A(3)*INRAY(1)+B(3).AND.
     +	      INRAY(2).GT.A(4)*INRAY(1)+B(4)
     +	      )THEN
	      MODHIT=MODLIST(MODINDEX)
	      ENDIF

	      IF (QUAD.EQ.3 .AND.
     +	      INRAY(2).GT.A(1)*INRAY(1)+B(1).AND.
     +	      INRAY(2).GT.A(2)*INRAY(1)+B(2).AND.
     +	      INRAY(2).LT.A(3)*INRAY(1)+B(3).AND.
     +	      INRAY(2).LT.A(4)*INRAY(1)+B(4)
     +	      )THEN
	      MODHIT=MODLIST(MODINDEX)
	      ENDIF

	      IF (QUAD.EQ.4 .AND.
     +	      INRAY(2).LT.A(1)*INRAY(1)+B(1).AND.
     +	      INRAY(2).GT.A(2)*INRAY(1)+B(2).AND.
     +	      INRAY(2).GT.A(3)*INRAY(1)+B(3).AND.
     +	      INRAY(2).LT.A(4)*INRAY(1)+B(4)
     +	      )THEN
	      MODHIT=MODLIST(MODINDEX)
	      ENDIF

	   ENDDO
	ELSE
	   MODHIT=MODLIST(1)
	ENDIF



* end of code:
	RETURN
	END
