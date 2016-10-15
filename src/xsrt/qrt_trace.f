*+QRT_TRACE Perform ray tracing
        SUBROUTINE QRT_TRACE(IDEB,RIRIS,IOPT,AREA,DSHFT,YBAR,ZBAR,RMS)
*IDEB        input        debugging level (0 none)
*RIRIS       input        radius about centre of detector for analysis
*                if 0.0 then no analysis of detected distribution
*IOPT        input  <0 save traced.dat and detected.dat files
*                   >0 adjust focus and save traced.dat and detected.dat
*                   only rays with IOPT reflections are used in adjustment
*                   =0 don't save or adjust focus
*AREA        output        detected area within RIRIS
*DSHFT       output        axial shift to optimum focus (0.0 if IOPT<=0)
*YBAR        output        y centroid of detected distribution
*ZBAR        output        z centroid of detected distribution
*RMS         output        rms radius of detected distribution
Cf2py  intent(in) IDEB,RIRIS,IOPT
Cf2py  intent(out) AREA,DSHFT,YBAR,ZBAR,RMS
        IMPLICIT NONE
        INTEGER IDEB,IOPT
        DOUBLE PRECISION RIRIS,AREA,DSHFT,YBAR,ZBAR,RMS
*-Author Dick Willingale 2012-May-2
        INCLUDE 'SRT_COM'
        INCLUDE 'QR_COM'
C
        IF(ISTAT.NE.0) RETURN
C
        IDEBUG=IDEB
C
        CALL SRT_TRACE(RIRIS,IOPT,AREA,DSHFT,YBAR,ZBAR,RMS,ISTAT)
        END
