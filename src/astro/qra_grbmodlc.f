*+QRA_GRBMODLC Generate GRB model light curve
        SUBROUTINE QRA_GRBMODLC(ND,NT,TIME,ELO,EHI,IPUL,RATES,FLUX)
        IMPLICIT NONE
        INTEGER ND,NT,IPUL
        DOUBLE PRECISION TIME(NT),ELO,EHI,RATES(ND,NT),FLUX(NT)
*ND        input        number of detector rates channels (should be 8)
*NT        input        number of time samples
*TIME        input        time samples (wrt to trigger)
*ELO        input        low energy of flux band keV
*EHI        input        high energy of flux band keV
*IPUL        input        0 afterglow, >0 individual pulse, -1 all, -2 prompt
*RATES        output        count rates (1-4 XRT, 5-8 BAT)
*FLUX        output        flux keV cm-2 s-1 in band elo-ehi
*-Author Dick Willingale 2013-Jan-24
        include 'QR_COM'
        include 'GRBS_COM'
C
        integer nr
        PARAMETER (nr=8)
        double precision rats(nr)
        integer i,j
C
        IF(ISTAT.NE.0) RETURN
C
        IF(ND.ne.nr) then
                write(*,*) 'qra_grblc - must be 8 rates channels',ND
                istat=1
                return
        ENDIF
        IF(npars.eq.0) then
                write(*,*) 'qra_grblc - model parameters not set'
                istat=1
                return
        ENDIF
C Loop for all time samples
        do i=1,NT
                call qra_grbrates(TIME(i),npars,pars,IPUL,ELO,EHI,
     +                nr,rats,FLUX(i))
                do j=1,nr
                        RATES(j,i)=rats(j)
                enddo
        enddo
        END
