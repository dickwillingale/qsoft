*+QRT_SLE Set up Schmidt Lobster Eye (derived from QRT_KBS)
        SUBROUTINE QRT_SLE(PCEN,PNOR,RAXI,FLEN,CWIDTH,CHEIGHT,
     +  PITCH,WALL,PL,NMIR,IDF,IQ,NMAX,RC,PC,TC,AC,NSET)
*PCEN    input      centre of telescope aperture
*PNOR    input      normal to aperture
*RAXI    input      reference axis at centre of aperture
*FLEN    input      focal length (-ve for second stack) = focus to front aperture
*CWIDTH  input      cell width = mirror width (its active part)
*CHEIGHT input      cell height
*PITCH   input      pitch of mirrors thickness included = pdd+wall
*WALL    input      mirror thickness
*PL      input      mirror axial length (depth)
*IDF     input      deformation index
*IQ      input      surface quality index
*NMAX    input      maximum number of module coordinate returned
*RC      output     radius of each module
*PC      output     azimuth of each module
*TC      output     rotation of each module
*AC      output     axial length of each module
*NSET    output     number of module coordinates returned
*                   currently, a single module is created always
Cf2py  intent(in)  PCEN,PNOR,RAXI,FLEN,CWIDTH,CHEIGHT,
Cf2py  intent(in)  PITCH,WALL,PL,NMIR,IDF,IQ,NMAX
Cf2py  intent(out) RC,PC,TC,AC,NSET
c
c Derived from QRT_KBS by Vladimir Tichy
c
c Last modify 24 Jan 2018 
C Changed surface types to 30 and 31 Dick Willingale 2018-Nov-21
c
        IMPLICIT NONE
        DOUBLE PRECISION PCEN(3),PNOR(3),RAXI(3),FLEN,RAPMIN,RAPMAX
        DOUBLE PRECISION CWIDTH, CHEIGHT,PL
        INTEGER NMIR
        DOUBLE PRECISION PITCH,WALL,PLMIN,PLMAX
        INTEGER IPACK,IDF,IQ,NMAX,NSET
        DOUBLE PRECISION RC(NMAX),PC(NMAX),TC(NMAX),AC(NMAX)
*
        INTEGER IDEF(2),J,NTHIS,NP
        DOUBLE PRECISION PP(50),OTHER(3),PNORM(3),RAXIS(3),SCEN(3)
        DOUBLE PRECISION RSPHERE
        DOUBLE PRECISION RMID
        INCLUDE 'QR_COM'
C
        IF(ISTAT.NE.0) RETURN
C Get surface deformation and quality indices
        IDEF(1)=IDF
        IDEF(2)=2
C
C >outer< radius of SLE stack
        RSPHERE=ABS(FLEN)*2.0-PL/2
        RMID=ABS(FLEN)*2.0-PL
C
C set up legacy parameters for right behaviour
CPLMIN   input      minimum axial length of K-B slots
CPLMAX   input      maximum axial length of K-B slots
cIPACK   input      packing 0 1 module, 1 sunflower, 2 cartesian, 3 wide
C                    field cartesian
CRAPMIN  input      minimum radius for aperture of constellation
CRAPMAX  input      maximum radius for aperture of constellation
        PLMIN=PL
        PLMAX=PL
        RAPMIN=0.0
        RAPMAX=RSPHERE*2.0
        IPACK=0
C normalise reference vectors
C Find centre of curvature of spherical aperture
        DO J=1,3
                PNORM(J)=PNOR(J)
                RAXIS(J)=RAXI(J)
                SCEN(J)=PCEN(J)-PNOR(J)*RSPHERE
        ENDDO
        CALL SRT_VNRM(PNORM,ISTAT)
        CALL SRT_VNRM(RAXIS,ISTAT)
C force reference axis to be at 90 degrees to normal
        CALL SRT_VCRS(PNORM,RAXIS,OTHER)
        CALL SRT_VNRM(OTHER,ISTAT)
        CALL SRT_VCRS(OTHER,PNORM,RAXIS)
        CALL SRT_VNRM(RAXIS,ISTAT)
C Pack parameters into single array
        DO J=1,3
                PP(J)=PNORM(J)
                PP(J+3)=RAXIS(J)
                PP(J+6)=SCEN(J)
        ENDDO
        PP(10)=RSPHERE
        PP(11)=RAPMAX
        PP(12)=FLEN
        PP(13)=PITCH
        PP(14)=WALL
        PP(15)=PLMIN
        PP(16)=PLMAX
C The next 2 are place markers for the x,y of the ray intersection in locals
        PP(17)=0.0
        PP(18)=0.0
C index of surface and surface quality
        CALL SRT_NSUR(NTHIS,ISTAT)
        PP(19)=NTHIS
        PP(20)=IQ
C number of mirrors
        PP(21)=NMIR
C        PP(21)=2
C angle distance between adjacent mirrors
        PP(22)=2.0*ASIN(PITCH/(2.0*RMID))
C angle position of first mirror
        PP(23)=-PP(22)*(NMIR-1)/2.0
        NP=23
C Set up constellation aperture
        CALL SRT_SETSLECOM(IPACK,FLEN,PLMIN,PLMAX,PITCH-WALL,
     +        RAPMIN,RAPMAX,CWIDTH,CHEIGHT,NMAX,RC,PC,TC,AC,NSET)
C Set parameters in common for sphere with radial/cart+return limits,
C type=30 (was 25)
C        CALL SRT_SETF(0,25,NP,PP,IDEF,0,0,-1,ISTAT)
        CALL SRT_SETF(0,31,NP,PP,IDEF,0,0,-1,ISTAT)
C Set place marker for rectangular aperture of K-B slot
        CALL SRT_SETF(0,1,13,PP,0,0,0,-1,ISTAT)
C Set place markers for 4 sides of K-B slot
        CALL SRT_SETF(0,32,16,PP,0,IQ,NTHIS+3,-1,ISTAT)
        CALL SRT_SETF(0,32,16,PP,0,IQ,NTHIS+4,-1,ISTAT)
        CALL SRT_SETF(0,32,16,PP,0,IQ,NTHIS+5,-1,ISTAT)
        CALL SRT_SETF(0,32,16,PP,0,IQ,NTHIS+2,-1,ISTAT)
C
c        WRITE(*,*) '*** QRT_LIST FROM QRT_SLE ***'
c        CALL QRT_LIST()
c        WRITE(*,*) '*****************************'
C
        END
