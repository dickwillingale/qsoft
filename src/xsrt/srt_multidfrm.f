*-SRT_MULTIDFRM      Deformation of surfaces - sum of more matrices
        SUBROUTINE SRT_MULTIDFRM(IDEF,X,Y,Z,DZDX,DZDY,ISTAT)
        IMPLICIT NONE
        INTEGER IDEF(2),ISTAT
        INTEGER NDI
        DOUBLE PRECISION X,Y,Z,DZDX,DZDY
        DOUBLE PRECISION AZ,ADZDX,ADZDY
        INTEGER I
*IDEF   input   IDEF(1)=main deformation index, IDEF(2)=internal index
*X,Y    input   position on surface (in surface coordinates)
*Z      output  deformation normal to surface 
*DZDX   output  gradient of deformation wrt X
*DZDY   output  gradient of deformation wrt Y
*ISTAT  in/out  returned status
*Note for a surface of revolution X is axial direction and Y is azimuth.
*-Author Dick Willingale 1996-Nov-21
* extended by Vladimir Tichy 
        INCLUDE 'SRT_COM'
        INTEGER ID
        INTEGER DFT
        INTEGER FC
C
        IF(ISTAT.NE.0) RETURN
C Format of deformation data in common IDFM(3,MAXDF) and IDFP(5,MAXDF)
C Deformation parameter values held in array DSAM
C IDFM(1,IDEF(1))       maximum number of sub-matrices
C IDFM(2,IDEF(1))       number of X sample positions
C IDFM(3,IDEF(1))       number of Y sample positions
C IDFP(1,IDEF(1))       index to X samples
C IDFP(2,IDEF(1))       index to Y samples
C IDFP(3,IDEF(1))       index to Z matrices
C IDFP(4,IDEF(1))       index to X gradient matrices
C IDFP(5,IDEF(1))       index to Y gradient matcices
C
        ID=IDEF(1)
        Z=0
        DZDX=0
        DZDY=0
        FC=1
C        WRITE(*,*) 'DEBUG: SRT_MULTIDFRM ID = 'ID
C        IF (NDI.LE.0) RETURN
c        WRITE(*,*) 'DEBUG: SRT_DFRM ID = ',ID
c        WRITE(*,*) 'DEBUG: SRT_DFRM IDFM(1,1) = ',IDFM(1,1)
C        IF(ID.LT.1.OR.ID.GT.MAXDF) THEN
C                    WRITE(*,*) 'SRT_MULTIDFRM error - starting index ',ID,' out of range'
C                ISTAT=1
C                RETURN
C        ENDIF
c
        DO WHILE (FC.GE.0)
            IF(ID.LT.1.OR.ID.GT.MAXDF) THEN
                WRITE(*,*) 'SRT_MULTIDFRM error - index ',ID,'run out of range'
                ISTAT=1
                RETURN
            ENDIF
            DFT=IDFT(ID)
            IF (DFT.GE.0) THEN
                FC=-1
            ELSE
                DFT=-DFT
            ENDIF
            IF(DFT.EQ.0) THEN
                AZ=0
                ADZDX=0
                ADZDY=0
            ELSEIF(DFT.EQ.1) THEN
c                WRITE(*,*) 'SRT_MULTIDFRM DEBUG: CALLING IN2D'
                    CALL SRT_IN2D(X,Y,IDEF(2),IDFM(1,ID),IDFM(2,ID),IDFM(3,ID),
     +                  DSAM(IDFP(1,ID)),DSAM(IDFP(2,ID)),DSAM(IDFP(3,ID)),
     +                  DSAM(IDFP(4,ID)),DSAM(IDFP(5,ID)),AZ,ADZDX,ADZDY,ISTAT)
                Z=Z+AZ
                DZDX=DZDX+ADZDX
                DZDY=DZDY+ADZDY
            ELSEIF(DFT.EQ.2) THEN
c                WRITE(*,*) 'SRT_MULTIDFRM DEBUG: CALLING DFRMSINE'
                    CALL SRT_DFRMSINE(X,Y,IDEF(2),IDFM(1,ID),IDFM(2,ID),IDFM(3,ID),
     +                  DSAM(IDFP(1,ID)),DSAM(IDFP(2,ID)),AZ,ADZDX,ADZDY,ISTAT)
                Z=Z+AZ
                DZDX=DZDX+ADZDX
                DZDY=DZDY+ADZDY
                ELSE
                    WRITE(*,*) 'SRT_MULTIDFRM error - unknown deformation type.', DFT
                ENDIF
            IF(ISTAT.NE.0) RETURN
            ID=ID+1
        ENDDO         
C        WRITE(*,*) 'SRT_MULTIDFRM DEBUG: returned :',Z,DZDX,DZDY
        END
