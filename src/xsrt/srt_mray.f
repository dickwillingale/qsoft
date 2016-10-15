*+SRT_MRAY	Make a ray
	SUBROUTINE SRT_MRAY(ITYPE,IDEF,SD,SP,SA,AP,AN,AR,PL,RPOS,RDIR,RQU,
     +	IQU,ISTAT)
	IMPLICIT NONE
	INTEGER ITYPE,IDEF,IQU,ISTAT
	DOUBLE PRECISION SD(3),SP(3),SA,AP(3),AN(3),AR(3),PL(4)
	DOUBLE PRECISION RPOS(3),RDIR(3),RQU(3)
*ITYPE	input	type of source
*IDEF	input	deformation index
*SD	input	direction of source
*SP	input	position of source
*SA	input	aperture area per ray
*AP	input	reference position in source aperture
*AN	input	normal to source aperture
*AR	input	reference axis in aperture
*PL	input	limits of aperture and any further parameters
*RPOS	output	position of ray in aperture
*RDIR	output	direction of ray
*RQU	output	quality factors of ray (wavelength, area, etc.)
*IQU	output	history of ray (integer for each point along ray)
*ISTAT	in/out	returned status
* ITYPE		Configuration/position	Aperture limits
*  1		infinity		radial	
*  2		infinity		cartesian
*  3		finite			radial
*  4		finite			cartesian
*  5		diffuse			radial
*  6		diffuse			cartesian
* If aperture limits cartesian then PL(1)=xmin, PL(2)=ymin, PL(3)=xmax,
* PL(4)=ymax
* If aperture limits radial then PL(1)=rmin, PL(2)=rmax
*-Author Dick Willingale 1996-Dec-4
	INCLUDE 'SRT_COM'
	DOUBLE PRECISION AY(3),X,Y,R,T,P,COSA,DX,DY
        DOUBLE PRECISION SYS_DRAND
	INTEGER J,ID
C
	IF(ISTAT.NE.0) RETURN
C Generate ray positions in aperture
	CALL SRT_VCRS(AN,AR,AY)
	IF(MOD(ITYPE,2).EQ.0) THEN
C Cartesian limits
		X=PL(1)+(PL(3)-PL(1))*SYS_DRAND()
		Y=PL(2)+(PL(4)-PL(2))*SYS_DRAND()
	ELSE
C Radial limits
		R=SQRT(PL(1)**2+(PL(2)**2-PL(1)**2)*SYS_DRAND())
		T=PI*2.0*SYS_DRAND()
		X=R*COS(T)
		Y=R*SIN(T)
	ENDIF
	DO J=1,3
		RPOS(J)=AP(J)+AR(J)*X+AY(J)*Y
	ENDDO
C Get tweaks from deformation array
	CALL SRT_NXTSHFT(IDEF,DX,DY,ISTAT)
C Generate ray direction
	ID=(ITYPE+1)/2
	IF(ID.EQ.1) THEN
C Source at infinity
		DO J=1,3
			RDIR(J)=SD(J)+AR(J)*DX+AY(J)*DY
		ENDDO
		CALL SRT_VNRM(RDIR,ISTAT)
	ELSEIF(ID.EQ.2) THEN
C Source at finite distance
		DO J=1,3
			RDIR(J)=RPOS(J)-(SP(J)+AR(J)*DX+AY(J)*DY)
		ENDDO
		CALL SRT_VNRM(RDIR,ISTAT)
	ELSE
C Diffuse source (isotropic)
		T=ACOS(1.0-SYS_DRAND())
		P=2.0*PI*SYS_DRAND()
		RDIR(1)=COS(T)
		RDIR(2)=SIN(T)*COS(P)
		RDIR(3)=SIN(T)*SIN(P)
	ENDIF
C Set quality and history
C Quality 1 will pick up a wavelength from reflection etc.
	RQU(1)=0.0
C Quality 2 is area carried by ray
	CALL SRT_VDOT(RDIR,AN,COSA)
	RQU(2)=SA*ABS(COSA)
C History flag -2 for source
	IQU=-2
	END
