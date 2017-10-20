*+QRB_BIOPEAKS        Extract peaks from mass spectra
        SUBROUTINE QRB_BIOPEAKS(NB,NS,NP,IPEAK,IHWID,FLX,VAR,BAK,
     +        MOVERZ,CFLX,PPOS,PFLX,PVAR,PBAK)
        IMPLICIT NONE
        INTEGER NB,NS,NP,IPEAK(NP),IHWID(NB)
        DOUBLE PRECISION FLX(NB,NS),VAR(NB,NS),BAK(NB,NS)
        DOUBLE PRECISION MOVERZ(NB),CFLX(NP)
        DOUBLE PRECISION PPOS(NP),PFLX(NP,NS),PVAR(NP,NS),PBAK(NP,NS)
Cf2py    intent(in) NB,NS,NP,IPEAK,IHWID,FLX,VAR,BAK,MOVERZ
Cf2py    intent(out) CFLX,PPOS,PFLX,PVAR,PBAK
*NB        input        number of mass spec bins
*NS        input        number of individuals
*NP        input        number of peaks
*IPEAK     input        peak indices (peaks centres in original bins)
*IHWID     input        array of half widths (in bins)
*FLX       input        array of mass spec fluxes
*VAR       input        array of mass spec variances
*BAK       input        array of mass spec background
*MOVERZ    input        array of bin m/z positions in spectrum
*CFLX      output        combined mean flux for each peak
*PPOS      output        peak positions
*PFLX      output        peak fluxes
*PVAR      output        peak variances
*PBAK      output        peak background
* QR version RW 2014-Feb-16
* qsoft version RW 2017-Jun-26
*-Author Dick Willingale 2005-Apr-21
        INCLUDE 'QR_COM'
        INTEGER J,I,N1,N2,K
        DOUBLE PRECISION TOTAL
C Check status
        IF(ISTAT.NE.0) RETURN
C Extract peak data using indices in original binning IPEAK()
C PPOS is array of peak positions in units moverz
        DO I=1,NP
                PPOS(I)=0.0
                N1=MAX(1,IPEAK(I)-IHWID(I))
                N2=MIN(NB,IPEAK(I)+IHWID(I))
C Loop for instances
                TOTAL=0.0
                DO J=1,NS
                        PFLX(I,J)=0.0
                        PVAR(I,J)=0.0
                        PBAK(I,J)=0.0
                        DO K=N1,N2
                                PFLX(I,J)=PFLX(I,J)+FLX(K,J)
                                PVAR(I,J)=PVAR(I,J)+VAR(K,J)
                                PBAK(I,J)=PBAK(I,J)+BAK(K,J)
                                PPOS(I)=PPOS(I)+FLX(K,J)*MOVERZ(K)
                        ENDDO
                        TOTAL=TOTAL+PFLX(I,J)
                ENDDO
                PPOS(I)=PPOS(I)/TOTAL
        ENDDO
C accumulate average spectrum CFLX
        DO I=1,NP
                CFLX(I)=0.0
                DO J=1,NS
                        CFLX(I)=CFLX(I)+PFLX(I,J)
                ENDDO
                CFLX(I)=CFLX(I)/NS
        ENDDO
        END
