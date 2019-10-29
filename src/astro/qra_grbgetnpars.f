*+QRA_GRBGETNPARS get total number of GRB light curve parameters set
        SUBROUTINE QRA_GRBGETNPARS(NP)
        IMPLICIT NONE
        INTEGER NP
*NP        output        total number of fit parameters 
*-Author Dick Willingale 2013-Jan-24
        include 'QR_COM'
        include 'GRBS_COM'
C
        IF(ISTAT.NE.0) RETURN
C
        NP=npars
        END
