*+QRT_BDIPOLE B field of a magnetic dipole
        SUBROUTINE QRT_BDIPOLE(D,PX,PY,PZ,GX,GY,GZ,FX,FY,FZ,BX,BY,BZ,R)
        IMPLICIT NONE
        DOUBLE PRECISION D,PX,PY,PZ,GX,GY,GZ,FX,FY,FZ,BX,BY,BZ,R
*D        input        dipole moment (Gauss cm3)
*PX       input        x position of dipole (cm)
*PY       input        y position of dipole (cm)
*PZ       input        z position of dipole (cm)
*GX       input        x direction cosine of dipole moment
*GY       input        y direction cosine of dipole moment
*GZ       input        z direction cosine of dipole moment
*FX       input        x position for field (cm)
*FY       input        y position for field (cm)
*FZ       input        z position for field (cm)
*BX       output       x magnetic field at position FX,FY,FZ (Gauss)
*BY       output       y magnetic field at position FX,FY,FZ (Gauss)
*BZ       output       z magnetic field at position FX,FY,FZ (Gauss)
*R        output       distance from dipole
Cf2py  intent(in) D,PX,PY,PZ,GX,GY,GZ,FX,FY,FZ
Cf2py  intent(out) BX,BY,BZ,R
*-Author Dick Willingale 1994-Sep-14
        DOUBLE PRECISION GHAT(3),RHAT(3),AHAT(3),THAT(3),CTH,STH,RR
        INCLUDE 'QR_COM'
C
        IF(ISTAT.NE.0) return
C Calculate distance R from dipole position to field position
        RHAT(1)=FX-PX
        RHAT(2)=FY-PY
        RHAT(3)=FZ-PZ
        GHAT(1)=GX
        GHAT(2)=GY
        GHAT(3)=GZ
        CALL SRT_VNRM(GHAT,ISTAT)
        R=SQRT(RHAT(1)**2+RHAT(2)**2+RHAT(3)**2)
C Calculate radial unit vector RHAT
        IF(R.GT.0.01D0) THEN
                RHAT(1)=RHAT(1)/R
                RHAT(2)=RHAT(2)/R
                RHAT(3)=RHAT(3)/R
        ELSE
C If too close then fix at 0.01 cm on axis
                R=0.01D0
                RHAT(1)=GHAT(1)
                RHAT(2)=GHAT(2)
                RHAT(3)=GHAT(3)
        ENDIF
C Calculate elevation sine and cosine of position wrt dipole direction
        CALL SRT_VDOT(RHAT,GHAT,CTH)
        IF(ABS(CTH).LT.0.9999D0) THEN
C Calculate azimuthal unit vector AHAT
                CALL SRT_VCRS(RHAT,GHAT,AHAT)
                STH=SQRT(AHAT(1)**2+AHAT(2)**2+AHAT(3)**2)
                CALL SRT_VNRM(AHAT,ISTAT)
C Calculate the tangential unit vector THAT
                CALL SRT_VCRS(RHAT,AHAT,THAT)
                CALL SRT_VNRM(THAT,ISTAT)
        ELSE
                CTH=DSIGN(1.0D0,CTH)
                STH=0.0D0
                THAT(1)=0.D0
                THAT(2)=0.D0
                THAT(3)=0.D0
        ENDIF
C Sum radial and tangential components of field
        RR=D/R**3
        BX=RR*(2.D0*CTH*RHAT(1)+STH*THAT(1))
        BY=RR*(2.D0*CTH*RHAT(2)+STH*THAT(2))
        BZ=RR*(2.D0*CTH*RHAT(3)+STH*THAT(3))
        END
