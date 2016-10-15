*+QRA_BAND        D.Band GRB continuum spectrum
        FUNCTION QRA_BAND(E,ALPHA,BETA,EC)
        DOUBLE PRECISION QRA_BAND,E,ALPHA,BETA,EC
*E                Energy keV
*ALPHA                first photon index
*BETA                second photon index
*EC                cut energy keV
*QRA_BAND        continuum at E (photons/keV)
*-Author Dick Willingale 2013-Mar-22
        DOUBLE PRECISION DIFF,EP,R
        DIFF=ALPHA-BETA
        EP=-DIFF*EC
        IF(E.LT.EP) THEN
                QRA_BAND=(E**(-ALPHA))*EXP(-E/EC)
        ELSE
                QRA_BAND=(E**(-BETA))*EXP(DIFF)*(EP**(-DIFF))
        ENDIF
        END
