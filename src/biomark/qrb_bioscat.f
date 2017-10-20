*+QRB_BIOSCAT        Calculate mean and scatter of peaks
        SUBROUTINE QRB_BIOSCAT(NS,NP,CFLX,PFLX,PCTS,PSIG,IMAX,PND)
        INTEGER NS,NP,IMAX(NP)
        DOUBLE PRECISION CFLX(NP),PFLX(NP,NS)
        DOUBLE PRECISION PCTS(NP),PSIG(NP),PND(NP,NS)
Cf2py    intent(in) NS,NP,CFLX,PFLX
Cf2py    intent(out) PCTS,PSIG,IMAX,PND
*NS        input        number of individuals
*NP        input        number of peaks
*CFLX      input        average spectrum
*PFLX      input        peak fluxes
*PCTS      output       mean flux in peaks
*PSIG      output       fractional rms variation in peaks
*IMAX      output       instance of maximum for peak
*PND       utput        difference between normalised peaks and average spectrum
* QR version RW 2014-Feb-16
* qsoft version RW 2017-Jun-26
*-Author Dick Willingale 2005-Aug-30
        INTEGER J,I
        DOUBLE PRECISION SS,PMAX
C
        DO I=1,NP
                PCTS(I)=0.0
                SS=0.0
                IMAX(I)=0
                PMAX=0.0
                DO J=1,NS
                        PCTS(I)=PCTS(I)+PFLX(I,J)
                        SS=SS+PFLX(I,J)**2
                        IF(PFLX(I,J).GT.PMAX) THEN
                                IMAX(I)=J
                                PMAX=PFLX(I,J)
                        ENDIF
                        PND(I,J)=(PFLX(I,J)-CFLX(I))
                ENDDO
                PCTS(I)=PCTS(I)/NS
                SS=SS/NS
                PSIG(I)=SQRT(SS/PCTS(I)**2-1.0)
        ENDDO
        END
