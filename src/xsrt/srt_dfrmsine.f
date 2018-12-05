*-SRT_DFRMSINE sine deformation
        SUBROUTINE SRT_DFRMSINE(X,Y,NS,NSMAX,NX,NY,XPAR,YPAR,
     +        Z,DZDX,DZDY,ISTAT)
        IMPLICIT NONE
        INTEGER NS,NSMAX,NX,NY
        DOUBLE PRECISION X,Y,XPAR(NX,NSMAX),YPAR(NY,NSMAX)
        DOUBLE PRECISION Z,DZDX,DZDY
        INTEGER ISTAT
*X        input        X position on surface (axial for conic)
*Y        input        Y position in surface (azimuth for conic)
*NS        input        sub-surface number (shell for Wolter I nest)
*NSMAX        input        maximum number of sub-surfaces held
*NX        input        number of x sine parameters (0 or 3)
*NY        input        number of y sine parameters (0 or 3)
*XPAR        input        x sine parameters: ampl., period, phase
*YPAR        input        y sine parameters: ampl., period, phase
*Z        output        deformation at X,Y
*DZDX        output        gradient of deformation wrt x at X,Y
*DZDY        output        gradient of deformation wrt y at X,Y
*ISTAT        in/out        returned status
*-Author Vladimir Tichy (2017)
*
        DOUBLE PRECISION DX,DY
        DOUBLE PRECISION Z0,Z1
        DOUBLE PRECISION AMP,LAM,PHI
C
        DOUBLE PRECISION SRTPI
        SRTPI=4.D0*DATAN(1.D0)
CC
        IF(ISTAT.NE.0) RETURN
C
        Z0=0
        Z1=0
        DX=0
        DY=0
c
c        WRITE(*,*) 'DEBUG SRT_DFRMSINE: NS = ',NS
c        WRITE(*,*) 'DEBUG SRT_DFRMSINE: PI = ',SRTPI
C    
        IF(NS.GT.0.AND.NS.LE.NSMAX) THEN
C
C        WRITE (*,*) 'CAU VOLE'
c
                IF(NX.EQ.3) THEN
                        AMP=XPAR(1,NS)
                        LAM=XPAR(2,NS)
                        PHI=XPAR(3,NS)
                        Z0=AMP*SIN(2*SRTPI*(X-PHI)/LAM)
                        DX=AMP*COS(2*SRTPI*(X-PHI)/LAM)*2*SRTPI/LAM
                ENDIF
C
                IF(NY.EQ.3) THEN
                        AMP=YPAR(1,NS)
                        LAM=YPAR(2,NS)
                        PHI=YPAR(3,NS)
                        Z1=AMP*SIN(2*SRTPI*(Y-PHI)/LAM)
                        DY=AMP*COS(2*SRTPI*(Y-PHI)/LAM)*2*SRTPI/LAM
                ENDIF
C
        ELSE
                WRITE(*,*) 'SRT_DFRMSINE error - internal index',NS,'  out of range'
                WRITE(*,*) 'It must be > 0 and <=',NSMAX
                ISTAT=1
                RETURN
        ENDIF
C
        Z=Z0+Z1
        DZDX=DX
        DZDY=DY
        END
