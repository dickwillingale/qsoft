*+SPF_FEABS	Iron absorption
	FUNCTION SPF_FEABS(E,IEDGE,CD)
	REAL SPF_FEABS,E,CD
	INTEGER IEDGE
*E		Energy keV
*IEDGE		Edge, 1-cold, 2-He like, 3-H like
*CD		Column density 10**16.52 cm-2 so that values equi. to hydrogen
*		column in units 10**21 cm-2 (assumes abundance ratio Fe/H of
*		10**(7.52-12)
*SPF_FEABS	Absorption factor at energy E
*-Author Dick Willingale 1986-Sep-4
C I have calculated parameters using Tucker stuff page 243-244. For cold
C Fe they agree with Morrison and McCammon and the Henke tabulation
C they claim to have used. I can't see why the higher ionization states
C should be far out. OK for fitting to low resolution pulse height spectra.
C Cold edge energy corrected to 7.11 1988-Oct-10 RW
C He-like and H-like corrected to 8.85 and 9.29  from
C Makishima 1986 (Tenerife proceedings) 1988-Oct-10 KNA.
C IRCO corrected to include L-shell cross-sections
C power law fit to the data of Reilman and Manson Ap.J supp ser 40 used.
C 1988-Dec-12 KNA.
	IF(IEDGE.EQ.1) THEN
		TAUL=0.043*CD*(0.73/E)**2.2
		TAUK=0.0012*CD*(7.11/E)**3.11
		IF ((E.GE.0.73).AND.(E.LT.7.11)) THEN
C			PEABS = TAUL
			PEABS=0.
		ELSE IF(E.GE.7.11) THEN
C			PEABS = TAUK+TAUL
			PEABS = TAUK
		ELSE
			PEABS=0.
		ENDIF
	ELSEIF(IEDGE.EQ.2) THEN
		IF(E.GE.8.85) THEN
			PEABS = 0.0011*CD*(8.85/E)**3.15
		ELSE
			PEABS=0.
		ENDIF
	ELSEIF(IEDGE.EQ.3) THEN
		IF(E.GE.9.29) THEN
			PEABS = 0.0011*CD*(9.29/E)**3.16
		ELSE
			PEABS=0.
		ENDIF
	ENDIF
	SPF_FEABS=EXP(-PEABS)
	END

