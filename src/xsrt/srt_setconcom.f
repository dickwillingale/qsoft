*+SRT_SETCONCOM    Set values in constellation common block
	SUBROUTINE SRT_SETCONCOM(ITMOD,RIN,ROUT,NMOD,RC,PC,TC,WC,HC,AC,
     +	CC,GC,ISTAT)
*ITMOD	input	module aperture type 1 rectangle, 2 sector
*RIN	input	inner radius of full aperture
*ROUT	input	outer radius of full aperture
*NMOD	input	number of modules
*RC,PC	input	radius and azimuth of each module
*TC	input	rotation angle of module
*WC	input	width of each module (x if ITMOD=1, radius if ITMOD=2)
*HC	input	height of each module (y if ITMOD=1, azimuth radians if ITMOD=2)
*		if HC=0 then set ITMOD=3 to specify circular aperture
*AC	input	axial length of each module
*CC	input	curvature signature of each module
*GC	input	grazing angle ratio for each module
*ISTAT	in/out	status 0 OK
      	IMPLICIT NONE
	INTEGER ITMOD,NMOD,ISTAT
	DOUBLE PRECISION RIN,ROUT
	DOUBLE PRECISION RC(NMOD),PC(NMOD),TC(NMOD)
	DOUBLE PRECISION WC(NMOD),HC(NMOD),AC(NMOD)
	DOUBLE PRECISION CC(NMOD),GC(NMOD)
*-Dick Willingale 2012-Jun-28
	INCLUDE 'SRT_COM'
	INTEGER J
C
	IF(ISTAT.NE.0) RETURN
C
	IF(NMOD.GT.maxlist) THEN
		WRITE(*,*) 'SRT_CONCOM error -# modules >',maxlist
		ISTAT=1
		RETURN
	ENDIF
	iptype=ITMOD
	if(HC(1).eq.0) then
		iptype=3
	endif
	inrad=RIN
	outrad=ROUT
	nrofelts=NMOD
	DO J=1,nrofelts
		rlist(J)=RC(J)
		philist(J)=PC(J)
		thlist(J)=TC(J)
		wlist(J)=WC(J)
		if(HC(J).gt.0) then
			hlist(J)=HC(J)
		else
			hlist(J)=WC(J)
		endif
	        llist(J)=AC(J)
	        clist(J)=CC(J)
	        glist(J)=GC(J)
	ENDDO
	END
