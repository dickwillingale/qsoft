*+QR_XYZTRANS	3 axis transformation
	SUBROUTINE QR_XYZTRANS(N,X,Y,Z,CEN,AER,XP,YP,ZP)
	INTEGER N
	DOUBLE PRECISION X(N),Y(N),Z(N),CEN(3),AER(3),XP(N),YP(N),ZP(N)
*N	input	number of points
*X	input	x-positions
*Y	input	y-positions
*Z	input	z-positions
*CEN	input	new origin
*AER	input	azimuth, elevation and roll of new x-axis radians
*		roll angle from Z-axis to +ve elev. +ve clockwise
*XP	ouput	new x
*YP	output	new y
*ZP	output	new z
*-Author: Dick Willingale 2012-Jun-28
	DOUBLE PRECISION XC,YC,ZC
	DOUBLE PRECISION MAT(3,3),BACK(3,3)
	INTEGER I
	CALL AX_DMAT(AER,AER(3),MAT,BACK)
	DO I=1,N
		XC=X(I)-CEN(1)
		YC=Y(I)-CEN(2)
		ZC=Z(I)-CEN(3)
	  	XP(I) = XC*MAT(1,1)+YC*MAT(2,1)+ZC*MAT(3,1)
	  	YP(I) = XC*MAT(1,2)+YC*MAT(2,2)+ZC*MAT(3,2)
	  	ZP(I) = XC*MAT(1,3)+YC*MAT(2,3)+ZC*MAT(3,3)
	ENDDO
	END
