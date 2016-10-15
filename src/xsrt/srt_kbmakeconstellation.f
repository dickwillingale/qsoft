	SUBROUTINE srt_KBmakeconstellation(IPACK,INRAD,OUTRAD,ELTSIZE,
     *	LISTLEN,RLIST,PHILIST,THLIST,NROFELTS)

* this file by FS, october 2008

* RLIST = LIST OF RADIUS VALUES
* RLISTLEN = LENGTH OF RLIST

* initialise:
	IMPLICIT NONE
	INTEGER LISTLEN,INDEX,COUNTER,NROFELTS,IPACK,NN,NX,NY,NDO
	DOUBLE PRECISION RLIST(LISTLEN),PHILIST(LISTLEN),THLIST(LISTLEN),
     *	R,PHI,HALFRIB,INRAD,OUTRAD,ELTSIZE,C,PHIGOLD,PI
        DOUBLE PRECISION X,Y,B

* initialize:
	PI=ASIN(1.D0)*2.D0
	PHIGOLD=2.3999632297286533
	HALFRIB=ELTSIZE/SQRT(2.D0)
	C=ELTSIZE/(PI/SQRT(2.D0))
	B=ELTSIZE*1.1
* IPACK determines type of constellation
	if(ipack.eq.1) then
* Sunflower packing
		NDO=LISTLEN
	elseif(ipack.eq.2.or.ipack.eq.3) then
* Cartesian packing
		NN=int(OUTRAD/ELTSIZE)*2
		NDO=MIN(NN*NN,LISTLEN)
	else
		HALFRIB=0.0
		NDO=1
	endif
* loop over the elements:
	NROFELTS=0
	DO INDEX=1,NDO
	  if(IPACK.eq.1) then
* this is the sunflower constellation rule:
	   	PHI=(INDEX-1)*PHIGOLD
	  	R=C*SQRT(PHI)
	  elseif(ipack.eq.2.or.ipack.eq.3) then
* Cartesian packing
		NX=mod(INDEX-1,NN)
		NY=(INDEX-1)/NN
		X=NX*B-NN/2*B+B*0.5
		Y=NY*B-NN/2*B+B*0.5
		PHI=ATAN2(Y,X)
		R=SQRT(X**2+Y**2)
	   else
* Single module
		X=(OUTRAD+INRAD)*0.5/SQRT(2.0)
		Y=X
		PHI=ATAN2(Y,X)
		R=SQRT(X**2+Y**2)
	   endif
* check for the element being inside the mirror bounds:
	   IF (R.GT.INRAD+HALFRIB.AND.R.LT.OUTRAD-HALFRIB) THEN
	      NROFELTS=NROFELTS+1
	      if(IPACK.eq.1) THEN
	              RLIST(NROFELTS)=R
	              PHILIST(NROFELTS)=MOD(PHI,PI*2.0D0)
		      THLIST(NROFELTS)=-PI*0.25D0
	      elseif(IPACK.eq.2) then
	              RLIST(NROFELTS)=R
	              PHILIST(NROFELTS)=PHI
		      IF(PHI.LT.-PI*0.5) THEN
		      		THLIST(NROFELTS)=-PHI-PI
		      ELSEIF(PHI.LT.0.0) THEN
		      		THLIST(nrOFELTS)=-PHI-PI*0.5
		      ELSEIF(PHI.LT.PI*0.5) THEN
		      		THLIST(NROFELTS)=-PHI
		      ELSE
		      		THLIST(NROFELTS)=-PHI+PI*0.5
		      ENDIF
	      elseif(IPACK.eq.3) then
	              RLIST(NROFELTS)=R
	              PHILIST(NROFELTS)=PHI
		      THLIST(NROFELTS)=-PHI
	      endif
	   ENDIF
	ENDDO
	RETURN
	END
