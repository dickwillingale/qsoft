*+AX_DONVCP	Conversion 3-vector cross-product
	SUBROUTINE AX_DONVCP(VA,VB,VC)
	DOUBLE PRECISION VA(3), VB(3), VC(3)
*VA	Input 	3-vector
*VB	Input 	3-vector
*VC	Output 	3-vector, VC = VA x VB
*-
*Author: Gavin Eadie, 1972;  Fortran-77 version by Clive Page, 1984 November.
	DOUBLE PRECISION VNORM
	VC(1) = VA(2)*VB(3) - VA(3)*VB(2)
	VC(2) = VA(3)*VB(1) - VA(1)*VB(3)
	VC(3) = VA(1)*VB(2) - VA(2)*VB(1)
	VNORM = SQRT(VC(1)**2 + VC(2)**2 + VC(3)**2)
	DO I = 1,3
	  VC(I) = VC(I) / VNORM
	END DO
	END
