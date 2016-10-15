*+SRT_CONE      Calculate parameters for cone
	SUBROUTINE SRT_CONE(XMIN,RMIN,XMAX,RMAX,AX,AR,RP,P)
	IMPLICIT NONE
	DOUBLE PRECISION XMIN,RMIN,XMAX,RMAX,AX(3),AR(3),RP(3),P(16)
*XMIN	input	axial distance of first plane
*RMIN	input	radius at XMIN
*XMAX	input	axial distance at second plane
*RMAX	input	radius at XMAX
*AX	input	axis
*AR	input	reference axis
*RP	input	reference position on axis corresponding to origin 
*P	output  parameters
* r**2=P(10).x**2+P(11).x+P(12) where x is axial distance from vertex
* P(10)=tan(theta)**2, P(11)=P(12)=0.0
* xmin=P(13), rmin=P(14), xmax=P(15), rmax=P(16)
*-Author Dick Willingale 1996-Nov-20
	DOUBLE PRECISION TANT
	INTEGER J
C Calculate tan of half angle
	TANT=(RMAX-RMIN)/(XMAX-XMIN)
C Set conic parameters and limits
	P(10)=TANT**2
	P(11)=(RMIN-XMIN*TANT)*2.0*TANT
	P(12)=(RMIN-XMIN*TANT)**2
	P(13)=XMIN
	P(14)=RMIN
	P(15)=XMAX
	P(16)=RMAX
C Set vectors
	DO J=1,3
		P(J)=AX(J)
	ENDDO
	DO J=1,3
		P(J+3)=AR(J)
	ENDDO
	DO J=1,3
		P(J+6)=RP(J)
	ENDDO
	END
