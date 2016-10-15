*+SRT_SETKBSCOM    Set values in constellation common block
	SUBROUTINE SRT_SETKBSCOM(IPACK,FL,PLMIN,PLMAX,SLOT,
     +	RMIN,RMAX,CSIZE,NMAX,RC,PC,TC,AC,NSET)
*IPACK	input	module packing 0 1 module, 1 sunflower, 2 cartesian
*FL	input	focal length
*PLMIN	input	minimum module length
*PLMAX	input	maximum module length
*SLOT	input	slot aperture width
*RMIN	input	minimum radius of aperture
*RMAX	input	maximum radius of aperture
*CSIZE	input	size of each module mm
*NMAX	input	maximum number of element positions returned
*RC,PC	output	radius and azimuth of each module
*TC	output	rotation angle of module
*AC	output	axial length of each module
*NSET	output	number of elements set
* IF IPACK=0 then single module with centre at x=y=(rmax+rmin)/2/SQRT(2)
      	IMPLICIT NONE
	INTEGER IPACK,NMAX,NSET
	DOUBLE PRECISION FL,PLMIN,PLMAX,SLOT
	DOUBLE PRECISION RMIN,RMAX,CSIZE,RC(NMAX),PC(NMAX),TC(NMAX)
	DOUBLE PRECISION AC(NMAX)
*-Dick Willingale 2012-Jun-28
	INCLUDE 'SRT_COM'
	INTEGER J,NDO
	DOUBLE PRECISION THETAG
C
	inrad=RMIN
	outrad=RMAX
	NDO=MIN(maxlist,NMAX)
	CALL srt_kbmakeconstellation(IPACK,inrad,outrad,CSIZE,
     +          NDO,rlist,philist,thlist,nrofelts)
	iptype=1
	DO J=1,nrofelts
		wlist(J)=CSIZE
		hlist(J)=CSIZE
		THETAG=ATAN2(rlist(J)/SQRT(2.0),ABS(FL))/2.0
	        llist(J)=MIN(MAX(PLMIN,SLOT/THETAG),PLMAX)
		IF(FL.LT.0.0) THEN
C If 2nd stack then decrease radius to correct position for top of module
			rlist(j)=rlist(j)*(FL*2.0)/(FL*2.0+PLMAX)
		ENDIF
	ENDDO
C 
	NSET=MIN(NMAX,nrofelts)
	DO J=1,NSET
		IF(J.LE.NMAX) THEN
			RC(J)=rlist(J)
			PC(J)=philist(J)
			TC(J)=thlist(J)
			AC(J)=llist(J)
		ENDIF
	ENDDO
	END
