*+QRT_BFIELD  Do calculation for BFIELD for array of dipoles and positions
        SUBROUTINE QRT_BFIELD(NPOLE,DM,PDX,PDY,PDZ,DDX,DDY,DDZ,
     +  NPOS,PX,PY,PZ,BX,BY,BZ,RMIN)
        IMPLICIT NONE
        INTEGER NPOLE,NPOS
        DOUBLE PRECISION DM(NPOLE),PDX(NPOLE),PDY(NPOLE),PDZ(NPOLE)
        DOUBLE PRECISION DDX(NPOLE),DDY(NPOLE),DDZ(NPOLE)
        DOUBLE PRECISION PX(NPOS),PY(NPOS),PZ(NPOS)
        DOUBLE PRECISION BX(NPOS),BY(NPOS),BZ(NPOS),RMIN
*NPOLE    input        number of dipoles
*DM       input        dipole moments (Gauss cm3)
*POX      input        x positions of dipoles (cm)
*POY      input        y positions of dipoles (cm)
*POZ      input        z positions of dipoles (cm)
*DDX      input        x direction cosines of dipole moments
*DDY      input        y direction cosines of dipole moments
*DDZ      input        z direction cosines of dipole moments
*NPOS     input        number of positions for field calculation
*PX       input        x positions for calculation
*PY       input        y positions for calculation
*PZ       input        z positions for calculation
*BX       output       x magnetic field at positions PX,PY,PZ (Gauss)
*BY       output       y magnetic field at positions PX,PY,PZ (Gauss)
*BZ       output       z magnetic field at positions PX,PY,PZ (Gauss)
*RMIN     output       minimum distance from dipoles
Cf2py  intent(in) NPOLE,DM,PDX,PDY,PDZ,DDX,DDY,DDZ,NPOS,PX,PY,PZ
Cf2py  intent(out) BX,BY,BZ,RMIN
*-Author Dick Willingale 1994-Sep-16
        DOUBLE PRECISION BBX,BBY,BBZ,R
        INTEGER J,I
        INCLUDE 'QR_COM'
C
        IF(ISTAT.NE.0) RETURN
C
        RMIN=-1.0
        DO J=1,NPOS
                BX(J)=0.0
                BY(J)=0.0
                BZ(J)=0.0
                DO I=1,NPOLE
                        CALL QRT_BDIPOLE(DM(I),PDX(I),PDY(I),PDZ(I),
     +                  DDX(I),DDY(I),DDZ(I),
     +                  PX(J),PY(J),PZ(J),BBX,BBY,BBZ,R)
                        BX(J)=BX(J)+BBX
                        BY(J)=BY(J)+BBY
                        BZ(J)=BZ(J)+BBZ
                        IF(RMIN.LT.0.0) THEN
                                RMIN=R
                        ELSE
                                RMIN=MIN(RMIN,R)
                        ENDIF
                ENDDO
        ENDDO
        END
