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
        INCLUDE 'SRT_COM'
        INTEGER IX,IY,IMZ,IMX,IMY,NXT
        INCLUDE 'QR_COM'
C
        IF(ISTAT.NE.0) RETURN
C 
        IF(IT.EQ.1) THEN
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
                ELSE
                        WRITE(*,*) 'QRT_DEFS - insufficient storage'
                        ISTAT=1
                        RETURN
                ENDIF
        ELSE
                WRITE(*,*) 'QRT_DEFS error - unknown deformation type'
                ISTAT=1
                RETURN
        ENDIF
C Set parameters in common
        CALL SRT_SETD(ID,NM,NX,NY,IX,IY,IMZ,IMX,IMY,ISTAT)
        END
