*+SRT_PRAY      Put ray to store
        SUBROUTINE SRT_PRAY(IU,MAXNPR,NPR,RPOS,RDIR,RQUA,IQU,ISTAT)
        IMPLICIT NONE
        INTEGER IU,MAXNPR,NPR,IQU(MAXNPR),ISTAT
        DOUBLE PRECISION RPOS(3,MAXNPR),RDIR(3),RQUA(3)
*IU     input   unit number or pointer of store
*MAXNPR input   maximum number of positions for ray
*NPR    input  number of positions set
*RPOS   input  positions along ray
*RDIR   input  direction of ray
*RQUA   input  quality of ray
*IQU    input  history of positions along ray
*ISTAT  in/out  returned status
*-Author Dick Willingale 1996-Dec-3
	INTEGER J,I
C
        IF(ISTAT.NE.0) RETURN
C      
        IF(IU.GT.0) THEN
C Put ray onto file
                WRITE(IU,IOSTAT=ISTAT) 
     +          NPR,((RPOS(I,J),I=1,3),J=1,NPR),RDIR,RQUA,
     +		(IQU(J),J=1,NPR)
                IF(ISTAT.NE.0) THEN
                        WRITE(*,*) 'SRT_GRAY error - IOSTAT',ISTAT
                ENDIF
	ELSE
		WRITE(*,*) 'SRT_PRAY error - illegal channel number'
		ISTAT=1
        ENDIF
        END
