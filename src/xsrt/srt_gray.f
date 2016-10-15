*+SRT_GRAY      Get next ray from store
        SUBROUTINE SRT_GRAY(IU,MAXNPR,NPR,RPOS,RDIR,RQUA,IQU,ISTAT)
        IMPLICIT NONE
        INTEGER IU,MAXNPR,NPR,IQU(MAXNPR),ISTAT
        DOUBLE PRECISION RPOS(3,MAXNPR),RDIR(3),RQUA(3)
*IU     input   unit number or pointer
*MAXNPR input   maximum number of positions for ray
*NPR    output  number of positions set
*RPOS   output  positions along ray
*RDIR   output  direction of ray
*RQUA   output  quality of ray
*IQU    output  history of positions along ray
*ISTAT  in/out  returned status
*-Author Dick Willingale 1996-Dec-3
	INCLUDE 'SRT_COM'
	INTEGER J,I,IP
C
        IF(ISTAT.NE.0) RETURN
C      
        IF(IU.GT.0) THEN
C Get ray from file
                READ(IU,IOSTAT=ISTAT) 
     +          NPR,((RPOS(I,J),I=1,3),J=1,NPR),RDIR,RQUA,
     +		(IQU(J),J=1,NPR)
                IF(ISTAT.NE.0) THEN
                        WRITE(*,*) 'SRT_GRAY error - no ray on file'
                ENDIF
        ELSE
C Generate ray
                NPR=1
		IP=ISRC(2)
		CALL SRT_MRAY(ISRC(1),ISRC(3),PAR(IP),PAR(IP+3),PAR(IP+6),
     +		PAR(IP+13),PAR(IP+7),PAR(IP+10),PAR(IP+16),
     +		RPOS,RDIR,RQUA,IQU,ISTAT)
        ENDIF
        END
