*+QRA_HABS absorption by Hydrogen column
        SUBROUTINE qra_habs(cd,ne,ekev,fact)
        integer ne
        real cd,ekev(ne),fact(ne)
Cf2py  intent(in) cd,ne,ekev
Cf2py  intent(out) fact
*cd        input        hydrogen column 10**21 cm-2
*ne        input        number of energy values
*ekev        input        array of energies keV
*fact        output        array of absorption factors
*-Author Dick Willingale 2014-Feb-13
        INCLUDE 'QR_COM'
C
        IF(ISTAT.NE.0) RETURN
C Generate absorption factors
        CALL SPX_HYABS(cd,ekev,ne,fact)
        END
