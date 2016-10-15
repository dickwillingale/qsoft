*+QRT_KBS Set up Si Kirkpatrick-Baez stack
        SUBROUTINE QRT_KBS(PCEN,PNOR,RAXI,IPACK,RAPMIN,RAPMAX,FLEN,
     +  CSIZE,PITCH,WALL,PLMIN,PLMAX,IDF,IQ,NMAX,RC,PC,TC,AC,NSET)
*PCEN    input      centre of telescope aperture
*PNOR    input      normal to aperture
*RAXI    input      reference axis at centre of aperture
*IPACK   input      packing 0 1 module, 1 sunflower, 2 cartesian, 3 wide
*                    field cartesian
*RAPMIN  input      minimum radius for aperture of constellation
*RAPMAX  input      maximum radius for aperture of constellation
*FLEN    input      focal length (-ve for 2nd stack)
*CSIZE   input      size of each module in constellation
*PITCH   input      pitch of K-B slots
*WALL    input      wall thickness of K-B slots
*PLMIN   input      minimum axial length of K-B slots
*PLMAX   input      maximum axial length of K-B slots
*IDF     input      deformation index
*IQ      input      surface quality index
*NMAX    input      maximum number of module coordinate returned
*RC      output     radius of each module
*PC      output     azimuth of each module
*TC      output     rotation of each module
*AC      output     axial length of each module
*NSET    output     number of module coordinates returned
* IF IPACK=0 then single module with centre at x=y=(RAPMIN+RAPMAX)/2/SQRT(2)
Cf2py  intent(in) PCEN,PNOR,RAXI,IPACK,RAPMIN,RAPMAX,FLEN
Cf2py  intent(in) CSIZE,PITCH,WALL,PLMIN,PLMAX,IDF,IQ,NMAX
Cf2py  intent(out) RC,PC,TC,AC,NSET
        IMPLICIT NONE
        DOUBLE PRECISION PCEN(3),PNOR(3),RAXI(3),FLEN,RAPMIN,RAPMAX
        DOUBLE PRECISION CSIZE,PITCH,WALL,PLMIN,PLMAX
        INTEGER IPACK,IDF,IQ,NMAX,NSET
        DOUBLE PRECISION RC(NMAX),PC(NMAX),TC(NMAX),AC(NMAX)
*-Author Dick Willingale 2012-Jun-28
        INTEGER IDEF(2),J,NTHIS,NP
        DOUBLE PRECISION PP(50),OTHER(3),PNORM(3),RAXIS(3),SCEN(3)
        DOUBLE PRECISION RSPHERE
        INCLUDE 'QR_COM'
C
        IF(ISTAT.NE.0) RETURN
C Get surface deformation and quality indices
        IDEF(1)=IDF
        IDEF(2)=2
C Find radius of curvature of spherical aperture
        IF(FLEN.LT.0.0) THEN
                RSPHERE=ABS(FLEN)*2.0
        ELSE
                RSPHERE=ABS(FLEN)*2.0+PLMAX
        ENDIF
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
        NP=20
C Set up constellation aperture
        CALL SRT_SETKBSCOM(IPACK,FLEN,PLMIN,PLMAX,PITCH-WALL,
     +        RAPMIN,RAPMAX,CSIZE,NMAX,RC,PC,TC,AC,NSET)
C Set parameters in common for sphere with radial/cart+return limits, type=25
        CALL SRT_SETF(0,25,NP,PP,IDEF,0,0,-1,ISTAT)
C Set place marker for rectangular aperture of K-B slot
        CALL SRT_SETF(0,1,13,PP,0,0,0,-1,ISTAT)
C Set place markers for 4 sides of K-B slot
        CALL SRT_SETF(0,5,16,PP,0,IQ,NTHIS+3,-1,ISTAT)
        CALL SRT_SETF(0,5,16,PP,0,IQ,NTHIS+4,-1,ISTAT)
        CALL SRT_SETF(0,5,16,PP,0,IQ,NTHIS+5,-1,ISTAT)
        CALL SRT_SETF(0,5,16,PP,0,IQ,NTHIS+2,-1,ISTAT)
        END
