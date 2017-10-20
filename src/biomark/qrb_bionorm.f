*+QRB_BIONORM        Normalise peaks using index
        SUBROUTINE QRB_BIONORM(NS,NP,CFLX,NPN,INORM,RNORM,
     +  PFLX,PVAR,PBAK)
        INTEGER NS,NP,NPN,INORM(NP)
        DOUBLE PRECISION CFLX(NP),RNORM(NS)
        DOUBLE PRECISION PFLX(NP,NS),PVAR(NP,NS),PBAK(NP,NS)
Cf2py    intent(in) NS,NP,CFLX,NPN,INORM
Cf2py    intent(out) RNORM,PFLX,PVAR,PBAK
*NS        input        number of individuals
*NP        input        number of peaks
*CFLX      inout        average peak intensities
*NPN       input        number of normalisation peaks
*INORM     input        index of normalisation peaks
*RNORM     output        normalisation for each instance
*PFLX      output        peak fluxes
*PVAR      output        peak variances
*PBAK      output        peak background
*-Author Dick Willingale 2005-Aug-30
        INTEGER J,I,II,NN
        DOUBLE PRECISION RN,RNN
C loop for instances
        DO J=1,NS
C find total in normalisation peaks and calculate normalisation factor
C We normalise to minimize the S/N weighted fractional error over
C the chosen normalisation peaks
                RN=0.0
                RNN=0.0
                NN=0
                DO I=1,NPN
                        II=INORM(I)
                        IF(CFLX(II).GT.0.0) THEN
C Poisson S/N weight
C                                RN=RN+PFLX(II,J)/SQRT(CFLX(II))
c                                RNN=RNN+PFLX(II,J)**2/(CFLX(II))**(3.0/2.0)
C Fractional weight
                                RN=RN+PFLX(II,J)/CFLX(II)
                                RNN=RNN+PFLX(II,J)**2/CFLX(II)**2
C total count weight
C                                RN=RN+PFLX(II,J)
C                                RNN=RNN+PFLX(II,J)**2/CFLX(II)
                                NN=NN+1
                        ENDIF
                ENDDO
                RN=RN/RNN
                RNN=RN*RN
                DO I=1,NP
                        PFLX(I,J)=PFLX(I,J)*RN
                        PVAR(I,J)=PVAR(I,J)*RNN
                        PBAK(I,J)=PBAK(I,J)*RN
                ENDDO
                RNORM(J)=RN
        ENDDO
        END
