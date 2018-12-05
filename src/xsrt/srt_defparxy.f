*+SRT_DEFPARXY        Set x,y deformation parameters
        SUBROUTINE SRT_DEFPARXY(NSM,NS,NX,NY,X,Y,XD,YD,ISTAT)
        IMPLICIT NONE
        INTEGER NSM,NS,NX,NY,ISTAT
        DOUBLE PRECISION X(NX),Y(NY)
        DOUBLE PRECISION XD(NX,NSM),YD(NY,NSM)
*NSM        input        total number of sub-matrices
*NS        input        index of sub-surface
*NX        input        number of x parameters
*NY        input        number of y parameters
*X        input        x values
*Y        input        y values
*XD        in/out        x values
*YD        in/out        y values
*ISTAT        in/out        returned status
*-Author Vladimir Tichy (2017)
        INTEGER J
C
        IF(ISTAT.NE.0) RETURN
C
        DO J=1,NX
                XD(J,NS)=X(J)
        ENDDO
        DO J=1,NY
                YD(J,NS)=Y(J)
        ENDDO
        END
