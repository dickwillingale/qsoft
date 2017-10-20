*+QRB_BIOMATCH        search for peak position matches in two lists
        SUBROUTINE QRB_BIOMATCH(N1,N2,PL1,PL2,MAXDIF,NMATCH,IMATCH)
        INTEGER N1,N2,NMATCH,IMATCH(N1)
        DOUBLE PRECISION PL1(N1),PL2(N2),MAXDIF
Cf2py    intent(in) N1,N2,PL1,PL2,MAXDIF
Cf2py    intent(out) NMATCH,IMATCH
*N1        input        number of peaks in list 1
*N2        input        number of peaks in list 2
*PL1       input        list of peak positions
*PL2       input        list of peak positions
*MAXDIF    input        maximum difference for a match
*NMATCH    output        number of matches found
*IMATCH    output        index of positions from list 2 which match list 1        
* QR version RW 2014-Feb-16
*-Author Dick Willingale 2005-Apr-28
        INCLUDE 'QR_COM'
        INTEGER K,J
C Check status
        IF(ISTAT.NE.0) RETURN
C
        NMATCH=0
        DO J=1,N1
                DO K=1,N2
                        IF(ABS(PL1(J)-PL2(K)).LT.MAXDIF) THEN
                                NMATCH=NMATCH+1
                                IMATCH(J)=K
                        ENDIF
                ENDDO
        ENDDO
        END
