*+QRT_SIPORE Set up Si pore optics
        SUBROUTINE QRT_SIPORE(PCEN,PNORM,RAXIS,FLEN,RPITCH,APITCH,
     +  WALL,NMOD,RM,PM,TM,WM,HM,AM,CM,GM,WFR,A2J,IDF,IQ)
*PCEN      input   position of centre of aperture (above join plane)
*RNORM     input   normal at centre of aperture
*RAXIS     input   reference axis at centre of aperture
*FLEN      input        focal length
*RPITCH    input        pore radial pitch
*APITCH    input        pore azimuthal pitch
*WALL      input        wall thickness
*NMOD      input        number of modules (rectangular apertures)
*RM        input        array of module radii
*PM        input        array of module azimuths (radians)
*TM        input        array of module rotations (radians normally 0)
*WM        input        array of module widths (radial mm)
*HM        input        array of module heights (azimuthal mm)
*AM        input        array of module lengths (axial pore length mm)
*CM        input        array of module curvature signatures
*                0 conical, 1 Wolter, 2 curve-plane, 3 constant...
*GM        input        array of module grazing angle ratios
*WRF       input        module frame width (surrounding module aperture)
*A2J       input        aperture to join plane axial distance
*IDF       input        deformation index
*IQ        input        reflecting surface quality
Cf2py  intent(in) PCEN,PNORM,RAXIS,FLEN,RPITCH,APITCH
Cf2py  intent(in) WALL,NMOD,RM,PM,TM,WM,HM,AM,CM,GM,WFR,A2J,IDF,IQ
        IMPLICIT NONE
        DOUBLE PRECISION PCEN(3),PNORM(3),RAXIS(3)
        DOUBLE PRECISION FLEN,RPITCH,APITCH,WALL,WFR,A2J
        INTEGER IQ,IDF,NMOD
        DOUBLE PRECISION RM(NMOD),PM(NMOD),TM(NMOD)
        DOUBLE PRECISION WM(NMOD),HM(NMOD),AM(NMOD)
        DOUBLE PRECISION CM(NMOD),GM(NMOD)
*-Author Dick Willingale 2012-Jul-12
        INTEGER IDEF(2),J,NTHIS,NP,MXP
        PARAMETER (MXP=500)
        DOUBLE PRECISION OTHER(3)
        DOUBLE PRECISION PP(MXP)
        INCLUDE 'QR_COM'
C
        IF(ISTAT.NE.0) RETURN
C
        DO J=1,3
                PP(J)=PCEN(J)
                PP(J+3)=PNORM(J)
                PP(J+6)=RAXIS(J)
        ENDDO
C RMIN and RMAX
        PP(10)=0.0
        PP(11)=0.0
        DO J=1,NMOD
                PP(11)=MAX(PP(11),RM(J)+WM(J)*2.0)
        ENDDO
        PP(12)=FLEN
        PP(13)=RPITCH
        PP(14)=APITCH
        PP(15)=WALL
        PP(16)=IQ
C The next 2 are place markers for the x,y of the ray intersection in locals
        PP(17)=0.0
        PP(18)=0.0
C index of surface and surface quality
        CALL SRT_NSUR(NTHIS,ISTAT)
        PP(19)=NTHIS
C Axial curvature of pores (now set per module using CMOD)
        PP(20)=0.0
C Module frame width
        PP(21)=WFR
C Axial distance from aperture to join plane
        PP(22)=A2J
        NP=22
C surface deformation
        IDEF(1)=IDF
        IDEF(2)=2
C normalise reference vectors
        CALL SRT_VNRM(PP(4),ISTAT)
        CALL SRT_VNRM(PP(7),ISTAT)
C force reference axis to be at 90 degrees to normal
        CALL SRT_VCRS(PP(4),PP(7),OTHER)
        CALL SRT_VNRM(OTHER,ISTAT)
        CALL SRT_VCRS(OTHER,PP(4),PP(7))
        CALL SRT_VNRM(PP(7),ISTAT)
C Set constellation common
        CALL SRT_SETCONCOM(1,PP(10),PP(11),NMOD,RM,PM,TM,WM,HM,AM,
     +  CM,GM,ISTAT)
C Set parameters in common for plane with radial/cart+return limits, type=24
        CALL SRT_SETF(0,24,NP,PP,IDEF,0,0,-1,ISTAT)
C Set place markers for grid aperture in front of module
        CALL SRT_SETF(0,1,13,PP,0,0,0,-1,ISTAT)
        CALL SRT_SETF(0,1,13,PP,0,0,0,-1,ISTAT)
C Set place marker for aperture of 1st pore
        CALL SRT_SETF(0,1,13,PP,0,0,0,-1,ISTAT)
C Set place markers for 4 sides of pore
        CALL SRT_SETF(0,5,16,PP,0,IQ,NTHIS+5,-1,ISTAT)
        CALL SRT_SETF(0,5,16,PP,0,IQ,NTHIS+6,-1,ISTAT)
        CALL SRT_SETF(0,5,16,PP,0,IQ,NTHIS+7,-1,ISTAT)
        CALL SRT_SETF(0,5,16,PP,0,IQ,NTHIS+4,-1,ISTAT)
C Set place marker for square aperture of 2nd pore
        CALL SRT_SETF(0,1,13,PP,0,0,0,-1,ISTAT)
C Set place markers for 4 sides of pore
        CALL SRT_SETF(0,5,16,PP,0,IQ,NTHIS+10,-1,ISTAT)
        CALL SRT_SETF(0,5,16,PP,0,IQ,NTHIS+11,-1,ISTAT)
        CALL SRT_SETF(0,5,16,PP,0,IQ,NTHIS+12,-1,ISTAT)
        CALL SRT_SETF(0,5,16,PP,0,IQ,NTHIS+9,-1,ISTAT)
        END
