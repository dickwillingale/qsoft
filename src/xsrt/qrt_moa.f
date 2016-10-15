*+QRT_MOA Set up a cylindrial MOA
        SUBROUTINE QRT_MOA(PCEN,PNO,RAX,RCUR,XYAP,PITCH,WALL,PLEN,
     +  IDF,IQ)
*PCEN       input        centre of MOA
*PNO        input        normal at centre of MOA
*RAX        input   cylindrical axis of MOA
*RCUR       input        radius of curvature of cylinder
*XYAP       input        half width of MOA aperture
*PITCH      input   pitch of slots
*WALL       input        wall thickness between slots
*PLEN       input   depth of slots (thickness MOA)
*IDF        input   deformation index
*IQ         input   surface quality index
Cf2py  intent(in) PCEN,PNO,RAX,RCUR,XYAP,PITCH,WALL,PLEN,IDF,IQ
        IMPLICIT NONE
        DOUBLE PRECISION PCEN(3),PNO(3),RAX(3),RCUR,XYAP
        DOUBLE PRECISION PITCH,WALL,PLEN
        INTEGER IDF,IQ
*-Author Dick Willingale 2012-Jun-28
        INTEGER IDEF(2),J,NTHIS,NP
        DOUBLE PRECISION OTHER(3),SCEN(3),PP(50)
        DOUBLE PRECISION PNORM(3),RAXIS(3)
        INCLUDE 'QR_COM'
C
        IF(ISTAT.NE.0) RETURN
C Get surface deformation and quality indices
        IDEF(1)=IDF
        IDEF(2)=2
C normalise reference vectors
        DO J=1,3
                PNORM(J)=PNO(J)
                RAXIS(J)=RAX(J)
        ENDDO
        CALL SRT_VNRM(PNORM,ISTAT)
        CALL SRT_VNRM(RAXIS,ISTAT)
C Find centre of cylinder
        SCEN(1)=PCEN(1)-PNORM(1)*RCUR
        SCEN(2)=PCEN(2)-PNORM(2)*RCUR
        SCEN(3)=PCEN(3)-PNORM(3)*RCUR
C force reference axis to be at 90 degrees to normal
        CALL SRT_VCRS(PNORM,RAXIS,OTHER)
        CALL SRT_VNRM(OTHER,ISTAT)
        CALL SRT_VCRS(OTHER,PNORM,RAXIS)
        CALL SRT_VNRM(RAXIS,ISTAT)
C generate other axis
        CALL SRT_VCRS(PNORM,RAXIS,OTHER)
C Pack parameters into single array
        DO J=1,3
                PP(J)=RAXIS(J)
                PP(J+3)=OTHER(J)
                PP(J+6)=SCEN(J)
        ENDDO
        PP(10)=0.0
        PP(11)=0.0
        PP(12)=RCUR**2
        PP(13)=-XYAP
        PP(14)=-XYAP
        PP(15)=XYAP
        PP(16)=XYAP
        PP(17)=PITCH
        PP(18)=WALL
        PP(19)=PLEN
C index of surface and surface quality
        CALL SRT_NSUR(NTHIS,ISTAT)
        PP(20)=NTHIS
        PP(21)=IQ
        NP=21
C Set parameters in common for cylinder with axial/circum limits type 23
        CALL SRT_SETF(0,23,NP,PP,IDEF,0,0,-1,ISTAT)
C Set place marker for square aperture of pore
        CALL SRT_SETF(0,1,13,PP,0,0,0,-1,ISTAT)
C Set place markers for 4 sides of pore
        CALL SRT_SETF(0,5,13,PP,0,IQ,NTHIS+3,-1,ISTAT)
        CALL SRT_SETF(0,5,13,PP,0,IQ,NTHIS+4,-1,ISTAT)
        CALL SRT_SETF(0,5,13,PP,0,IQ,NTHIS+5,-1,ISTAT)
        CALL SRT_SETF(0,5,13,PP,0,IQ,NTHIS+2,-1,ISTAT)
        END
