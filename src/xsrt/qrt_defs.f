*+QRT_DEFS Set surface deformation parameters
        SUBROUTINE QRT_DEFS(ID,IT,NM,NX,NY)
        IMPLICIT NONE
*ID        input        deformation index
*IT        input        deformation type (1 matrix)
*NM        input        number of sub-matrices
*NX        input        number of x samples
*NY        input        number of y samples
Cf2py  intent(in) ID,IT,NM,NX,NY
        INTEGER ID,IT,NM,NX,NY
*-Author Dick Willingale 2012-May-15
C Modified by Vladimir Tichy (2017)
        INCLUDE 'SRT_COM'
        INTEGER IX,IY,IMZ,IMX,IMY,NXT,AIT
        INCLUDE 'QR_COM'
C
        IF(ISTAT.NE.0) RETURN
C 
C Check deformation index
        IF(ID.LT.0.OR.ID.GT.MAXDF) THEN
          WRITE(*,*) 'QRT_DEFS error - deformation index out of range'
          ISTAT=1
          RETURN
        ENDIF
        IDFT(ID)=IT
c        WRITE(*,*) 'DEBUG QRT_DEFS: IDD           = ',ID
c        WRITE(*,*) 'DEBUG QRT_DEFS: IT            = ',IT
C        WRITE(*,*) 'DEBUG QRT_DEFS: IT (backread) = ',IDFT(ID)
c        WRITE(*,*) 'DEBUG QRT_DEFS: NM            = ',NM
c        WRITE(*,*) 'DEBUG QRT_DEFS: NX            = ',NX
c        WRITE(*,*) 'DEBUG QRT_DEFS: NY            = ',NY
c        WRITE(*,*) 'DEBUG QRT_DEFS: IT=',IT
c 0 = no deformation
c 1 = deformation matrix
c 2 = sine wave
c negative value = scan also for next deformation and sum them together
        AIT=ABS(IT)
        IF(AIT.EQ.1.OR.AIT.EQ.222) THEN
C Set indicies
                NXT=IDSAM+NX*NM+NY*NM+NX*NY*NM*3
                IF(NXT.LE.MAXDSAM) THEN
                    IX=IDSAM+1
                    IDSAM=IDSAM+NX*NM
                    IY=IDSAM+1
                    IDSAM=IDSAM+NY*NM
                    IMZ=IDSAM+1
                    IDSAM=IDSAM+NX*NY*NM
                    IMX=IDSAM+1
                    IDSAM=IDSAM+NX*NY*NM
                    IMY=IDSAM+1
                    IDSAM=IDSAM+NX*NY*NM
C                      Set parameters in common
                    CALL SRT_SETD(ID,NM,NX,NY,IX,IY,IMZ,IMX,IMY,ISTAT)
                ELSE
                    WRITE(*,*) 'QRT_DEFS - insufficient storage'
                    ISTAT=1
                    RETURN
                ENDIF
        ELSEIF (AIT.EQ.2) THEN
C Sine deformation
                NXT=IDSAM+NX*NM+NY*NM
                IF(NXT.LE.MAXDSAM) THEN
                        IX=IDSAM+1
                        IDSAM=IDSAM+NX*NM
                        IY=IDSAM+1
                        IDSAM=IDSAM+NY*NM
C                         Parameters are saved as "x-samples"
C                       and "y-samples"
                        CALL SRT_SETD(ID,NM,NX,NY,IX,IY,0,0,0,ISTAT)
                ELSE
                        WRITE(*,*) 'QRT_DEFS - insufficient storage'
                        ISTAT=1
                        RETURN
                ENDIF
        ELSEIF (AIT.NE.0) THEN
C 0 = no deformations (empty deformation)
                WRITE(*,*) 
     +  'QRT_DEFS error - unknown deformation type',AIT
                ISTAT=1
                RETURN
        ENDIF
C        IDFT(ID)=IT
c        WRITE(*,*) 'DEBUG QRT_DEFS 2 : IDD           = ',ID
c        WRITE(*,*) 'DEBUG QRT_DEFS 2 : IT            = ',IT
c        WRITE(*,*) 'DEBUG QRT_DEFS 2 : IT (backread) = ',IDFT(ID)
c        WRITE(*,*) 'DEBUG QRT_DEFS 2 : NM            = ',NM
        END
        
