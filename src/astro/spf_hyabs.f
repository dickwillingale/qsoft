*+SPF_HYABS	Photoelectric absorption by interstellar medium
	FUNCTION SPF_HYABS(E,CD)
	REAL SPF_HYABS,E,CD
*E		Energy keV
*CD		Column density 10**21 cm-2
*SPF_HYABS	Absorption factor at energy E
*Uses piecewise polynomial fit of Morrison and McCammon Ap.J. 270, 119
*for range 0.03 to 10 keV. Below 0.03 keV uses power law fit to hydrogen and
*helium edge profiles interpolated/extrapolated from Henke data (1982) also
*used by Morrisom and McCammon.
*Above 10 keV crude "eyeball" fit provided by Gordon Stewart!?
*
*1992-Nov-18 RW. I've improved the fit to Henke data (1982)
*to match exactly with Morrison and McCammon at 30.5 eV.
*Originally I used the Cromer and Liberman crosssection for H
*for E<24.6 eV.
*Note. There is a considerable discrepency between the Cromer and Liberman and
*Henke data for H and HE in the range 13.6 to 30.5 eV. The Henke absorption
*edges are markedly smaller. This routine now (1992-Nov-18) uses Henke alone
*which makes it entirely consistent with Morrison and McCammon.
*-Author Dick Willingale 1986-Sep-4
	REAL C0(14),C1(14),C2(14),ET(15)
	DATA C0/17.3,34.6,78.1,71.4,95.5,308.9,120.6,141.3,
     +		202.7,342.7,352.2,433.9,629.0,701.2/
	DATA C1/608.1,267.9,18.8,66.8,145.8,-380.6,169.3,146.8,
     +		104.7,18.7,18.7,-2.4,30.9,25.2/
	DATA C2/-2150.,-476.1,4.3,-51.4,-61.1,294.,-47.7,-31.5,
     +		-17.,0.,0.,0.75,0.,0./
	DATA ET/.03,.1,.284,.4,.532,.707,.867,1.303,1.84,2.471,
     + 		3.21,4.038,7.111,8.331,10./
	IF(E.LT.0.0136) THEN
		TAU=0.
	ELSEIF(E.GE.0.0136.AND.E.LT.0.0246) THEN
		TAU=CD*21.05E-3/(E**2.976)
	ELSEIF(E.GE.0.0246.AND.E.LT.0.03) THEN
		TAU=CD*(21.05E-3/(E**2.976)+320.1E-3/(E**2.113))
	ELSEIF(E.GE.10.) THEN
		TAU=(0.30*CD)/(E*E*SQRT(E))
	ELSE
		DO K=1,14
			IF(E.GE.ET(K).AND.E.LT.ET(K+1)) THEN
			     TAU=(C0(K)+C1(K)*E+C2(K)*E*E)/E/E/E
			     TAU=0.001*TAU*CD
			ENDIF
		ENDDO
	ENDIF
	SPF_HYABS=EXP(-TAU)
	END
