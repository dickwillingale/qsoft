*+SRT_SPH1      Calculate parameters for spherical approx. of Wolter Type I
	SUBROUTINE SRT_SPH1(RWIDTH,PP,PY,PPS,PYS,ISTAT)
	IMPLICIT NONE
	INTEGER ISTAT
	DOUBLE PRECISION RWIDTH,PP(16),PY(16),PPS(14),PYS(14)
*RWIDTH	input	width of spherical surface along tangential reference
*PP     input	parameters for the paraboloid
*PY     input	parameters for the hyperboloid
*PPS    output	parameters for the sphere 1
*PYS    output	parameters for the sphere 2
*ISTAT	in/out	returned status
*-Author Dick Willingale 1997-Oct-30
	DOUBLE PRECISION AXP(3),ARP(3),FFP(3),AXH(3),ARH(3),FFH(3)
	DOUBLE PRECISION SXP(3),SRP(3),SFP(3),SXH(3),SRH(3),SFH(3)
	DOUBLE PRECISION AP,BP,CP,AH,BH,CH,XP,XH,RP,RH
	DOUBLE PRECISION XP0,XH0,RP0,RH0,RADP,RADH,FU,FQ
	DOUBLE PRECISION TTT(3)
	INTEGER J
C Get vectors that define coordinate system of Wolter I
	DO J=1,3
		AXP(J)=PP(J)
		AXH(J)=PY(J)
	ENDDO
	DO J=1,3
		ARP(J)=PP(J+3)
		ARH(J)=PY(J+3)
	ENDDO
	DO J=1,3
		FFP(J)=PP(J+6)
		FFH(J)=PY(J+6)
	ENDDO
C Get coefficients of conic sections
	AP=PP(10)
	BP=PP(11)
	CP=PP(12)
	AH=PY(10)
	BH=PY(11)
	CH=PY(12)
C Calculate datum for parabola and hyperbola
	XP=(PP(15)+PP(13))*0.5
	XH=(PY(15)+PY(13))*0.5
	RP=SQRT(AP*XP**2+BP*XP+CP)
	RH=SQRT(AH*XH**2+BH*XH+CH)
	write(*,*) 'xp,rp',xp,rp
	write(*,*) 'xh,rh',xh,rh
C Calculate parameters for spheres
	RP0=-4.0*RP**3/BP**2
	XP0=BP*0.5+2.0*RP**2/BP+XP
	write(*,*) 'rp0,xp0',rp0,xp0
	FU=AH/RH-AH**2*XH**2/RH**3-BH**2/(4.0*RH**3)-BH*AH*XH/RH**3
	FQ=AH*XH/RH+BH/(2.0*RH)
	RH0=RH-(FQ**2-1)/FU
	XH0=(RH-RH0)*FQ+XH
	write(*,*) 'rh0,xh0',rh0,xh0
	RADP=SQRT((XP-XP0)**2+(RP-RP0)**2)
	RADH=SQRT((XH-XH0)**2+(RH-RH0)**2)
C Calculate axes of spheres
	DO J=1,3
		SXP(J)=AXP(J)*(XP-XP0)+ARP(J)*(RP-RP0)
		SXH(J)=AXH(J)*(XH-XH0)+ARH(J)*(RH-RH0)
	ENDDO
	CALL SRT_VNRM(SXP,ISTAT)
	CALL SRT_VNRM(SXH,ISTAT)
C Calculate tangential reference axes of spheres
	CALL SRT_VCRS(AXP,ARP,TTT)
	CALL SRT_VCRS(SXP,TTT,SRP)
	CALL SRT_VCRS(AXH,ARH,TTT)
	CALL SRT_VCRS(SXH,TTT,SRH)
C Calculate centres of spheres
	DO J=1,3
		SFP(J)=FFP(J)+AXP(J)*XP0+ARP(J)*RP0
		SFH(J)=FFH(J)+AXH(J)*XH0+ARH(J)*RH0
	ENDDO
C Pack parameters into arrays
	DO J=1,3
		PPS(J)=SXP(J)
		PYS(J)=SXH(J)
		PPS(J+3)=SRP(J)
		PYS(J+3)=SRH(J)
		PPS(J+6)=SFP(J)
		PYS(J+6)=SFH(J)
	ENDDO
	PPS(10)=RADP
	PYS(10)=RADH
C Set surface limits parabola sphere
	PPS(11)=-(PP(15)-PP(13))*0.5
	PPS(12)=-RWIDTH*0.5
	PPS(13)=(PP(15)-PP(13))*0.5
	PPS(14)=RWIDTH*0.5
C Set surface limits hyperbola sphere
	PYS(11)=-(PY(15)-PY(13))*0.5
	PYS(12)=-RWIDTH*0.5
	PYS(13)=(PY(15)-PY(13))*0.5
	PYS(14)=RWIDTH*0.5
	END
