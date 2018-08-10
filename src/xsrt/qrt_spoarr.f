*+QRT_SPOARR Set up array of  Si pore optics modules
        SUBROUTINE QRT_SPOARR(PCEN,PNORM,RAXIS,FLEN,A2J,NMOD,
     +  RM,PM,TM,WM,HM,AM,CM,GM,RPITCH,WALL,APITCH,RWI,WFR,SIQ,IDF)
*PCEN      input   position of centre of aperture (above join plane)
*PNORM     input   normal at centre of aperture
*RAXIS     input   reference axis at centre of aperture
*FLEN      input   focal length of array
*A2J       input   aperture to join plane (axial mm)
*NMOD      input   number of modules (rectangular apertures)
*RM        input   array of module radial positions (mm)
*PM        input   array of module azimuthal positions (radians)
*TM        input   array of module rotations (radians normally 0)
*WM        input   array of module widths (radial mm)
*HM        input   array of module heights (azimuthal mm)
*AM        input   array of module lengths (axial pore length mm)
*CM        input   array of module curvature indices
*GM        input   array of module grazing angle ratios
*RPITCH    input   array of module pore radial pitch (mm)
*WALL      input   array of module radial wall thickness (mm)
*APITCH    input   array of module pore azimuthal pitch (rib spacing mm)
*RWI       input   array of module pore rib thickness (mm)
*WFR       input   array of module frame widths (mm)
*SIQ       input   array of module reflecting surface quality indices
*IDF       input   deformation index
Cf2py  intent(in) PCEN,PNORM,RAXIS,FLEN,NMOD
Cf2py  intent(in) RM,PM,TM,WM,HM,AM,CM,GM
Cf2py  intent(in) RPITCH,WALL,APITCH,RWI,SIQ,WFR,A2J,IDF
        IMPLICIT NONE
        DOUBLE PRECISION PCEN(3),PNORM(3),RAXIS(3)
        DOUBLE PRECISION FLEN,A2J
        INTEGER IDF,NMOD
        DOUBLE PRECISION RM(NMOD),PM(NMOD),TM(NMOD)
        DOUBLE PRECISION WM(NMOD),HM(NMOD),AM(NMOD)
        DOUBLE PRECISION CM(NMOD),GM(NMOD)
        DOUBLE PRECISION RPITCH(NMOD),APITCH(NMOD),WALL(NMOD)
        DOUBLE PRECISION RWI(NMOD),SIQ(NMOD),WFR(NMOD)
*-Author Dick Willingale 2018-Mar-26
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
        PP(13)=A2J
        PP(14)=0.0
        PP(15)=0.0
        PP(16)=0.0
C The next 2 are place markers for the x,y of the ray intersection in locals
        PP(17)=0.0
        PP(18)=0.0
C index of surface
        CALL SRT_NSUR(NTHIS,ISTAT)
        PP(19)=NTHIS
        NP=19
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
        CALL SRT_SETSPOCOM(NMOD,RM,PM,TM,WM,HM,AM,
     +  CM,GM,RPITCH,WALL,APITCH,RWI,WFR,SIQ,ISTAT)
C Set parameters in common for plane with radial/cart+return limits,
C type=30
        CALL SRT_SETF(0,30,NP,PP,IDEF,0,0,-1,ISTAT)
C Set place markers for grid aperture in front of module
        CALL SRT_SETF(0,1,13,PP,0,0,0,-1,ISTAT)
        CALL SRT_SETF(0,1,13,PP,0,0,0,-1,ISTAT)
C Set place marker for aperture of 1st pore
        CALL SRT_SETF(0,1,13,PP,0,0,0,-1,ISTAT)
C Set place markers for 4 sides of pore
        CALL SRT_SETF(0,5,16,PP,0,0,NTHIS+5,-1,ISTAT)
        CALL SRT_SETF(0,5,16,PP,0,0,NTHIS+6,-1,ISTAT)
        CALL SRT_SETF(0,5,16,PP,0,0,NTHIS+7,-1,ISTAT)
        CALL SRT_SETF(0,5,16,PP,0,0,NTHIS+4,-1,ISTAT)
C Set place marker for square aperture of 2nd pore
        CALL SRT_SETF(0,1,13,PP,0,0,0,-1,ISTAT)
C Set place markers for 4 sides of pore
        CALL SRT_SETF(0,5,16,PP,0,0,NTHIS+10,-1,ISTAT)
        CALL SRT_SETF(0,5,16,PP,0,0,NTHIS+11,-1,ISTAT)
        CALL SRT_SETF(0,5,16,PP,0,0,NTHIS+12,-1,ISTAT)
        CALL SRT_SETF(0,5,16,PP,0,0,NTHIS+9,-1,ISTAT)
        END
