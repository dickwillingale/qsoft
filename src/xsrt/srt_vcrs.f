*+SRT_VCRS        Vector cross product
	SUBROUTINE SRT_VCRS(VA,VB,VC)
	DOUBLE PRECISION VA(3),VB(3),VC(3)
*-Author Dick Willingale 1996-Nov-16
	VC(1)=VA(2)*VB(3)-VA(3)*VB(2)
	VC(2)=VA(3)*VB(1)-VA(1)*VB(3)
	VC(3)=VA(1)*VB(2)-VA(2)*VB(1)
	END
