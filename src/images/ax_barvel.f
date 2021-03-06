*+AX_BARVEL	Returns barycentric/heliocentric velocity cmpts of Earth.
	SUBROUTINE AX_BARVEL(DMJDE,DEQ,DVELH,DVELB)
	IMPLICIT DOUBLE PRECISION(D)
	DIMENSION DVELH(3),DVELB(3)
*DMJDE	Input	Specifies Modified Julian Date (ephemeris time).
*DEQ	Input	Specifies epoch of mean equator and equinox to use,
*		e.g. 1950.0D0; if DEQ = 0D0 then the mean equator
*		and equinox of DJE will be used.
*DVELH	Output	Heliocentric velocity components
*		for dx/dt, dy/dt, dz/dt in A.U. per second.
*DVELB	Output	Barycentric components as above.
*-
*    The COMMON block /BARXYZ/ contains those intermediate results of
*    BARVEL which can be used to compute the heliocentric and
*    barycentric coordinates (subroutine BARCOR).
* Author	P. Stumpff, 1979 July 15,	Astron & Astrophys
*  Converted to standard Fortran by C.G.Page.
*Modified CGP 1985 Jan 24: takes JD or MJD of epoch.
*
	DIMENSION SN(4),DCFEL(3,8),DCEPS(3),CCSEL(3,17),DCARGS(2,15),
     &    CCAMPS(5,15),CCSEC(3,4),DCARGM(2,3),CCAMPM(4,3),CCPAMV(4)
	EQUIVALENCE (SORBEL(1),E),(FORBEL(1),G)

	COMMON /BARXYZ/ DPREMA(3,3), DPSI, D1PDRO, DSINLS, DCOSLS,
     &   DSINEP, DCOSEP, FORBEL(7), SORBEL(17), SINLP(4), COSLP(4),
     &   SINLM, COSLM, SIGMA, IDEQ

	PARAMETER (DC2PI = 6.2831853071796D0, CC2PI = 6.283185,
     &	 DC1 = 1.0D0, DCT0 = 2415020.0D0, DCJUL = 36525.0D0,
     &   DCBES = 0.313D0, DCTROP = 365.24219572D0,
     &   DC1900 = 1900.0D0)

* SIDEREAL RATE DCSLD IN LONGITUDE, RATE CCSGD IN MEAN ANOMALY
	PARAMETER (DCSLD = 1.990987D-07, CCSGD = 1.990969E-07)

* SOME CONSTANTS USED IN THE CALCULATION OF THE LUNAR CONTRIBUTION
	PARAMETER(CCKM = 3.122140E-05, CCMLD = 2.661699E-06,
     &   CCFDI = 2.399485E-07)

* CONSTANTS DCFEL(I,K) OF FAST CHANGING ELEMENTS
*		       I = 1		I = 2		I = 3
	DATA DCFEL/ 1.7400353D+00, 6.2833195099091D+02, 5.2796D-06,
     &		    6.2565836D+00, 6.2830194572674D+02,-2.6180D-06,
     &              4.7199666D+00, 8.3997091449254D+03,-1.9780D-05,
     &              1.9636505D-01, 8.4334662911720D+03,-5.6044D-05,
     &              4.1547339D+00, 5.2993466764997D+01, 5.8845D-06,
     &              4.6524223D+00, 2.1354275911213D+01, 5.6797D-06,
     &              4.2620486D+00, 7.5025342197656D+00, 5.5317D-06,
     &              1.4740694D+00, 3.8377331909193D+00, 5.6093D-06/

* CONSTANTS DCEPS AND CCSEL(I,K) OF SLOWLY CHANGING ELEMENTS
*		      I = 1           I = 2           I = 3
	DATA DCEPS/ 4.093198D-01,-2.271110D-04,-2.860401D-08/
	DATA CCSEL/ 1.675104E-02,-4.179579E-05,-1.260516E-07,
     &		    2.220221E-01, 2.809917E-02, 1.852532E-05,
     &              1.589963E+00, 3.418075E-02, 1.430200E-05,
     &              2.994089E+00, 2.590824E-02, 4.155840E-06,
     &              8.155457E-01, 2.486352E-02, 6.836840E-06,
     &              1.735614E+00, 1.763719E-02, 6.370440E-06,
     &		    1.968564E+00, 1.524020E-02,-2.517152E-06,
     &              1.282417E+00, 8.703393E-03, 2.289292E-05,
     &		    2.280820E+00, 1.918010E-02, 4.484520E-06,
     &		    4.833473E-02, 1.641773E-04,-4.654200E-07,
     &		    5.589232E-02,-3.455092E-04,-7.388560E-07,
     &		    4.634443E-02,-2.658234E-05, 7.757000E-08,
     &		    8.997041E-03, 6.329728E-06,-1.939256E-09,
     &		    2.284178E-02,-9.941590E-05, 6.787400E-08,
     &		    4.350267E-02,-6.839749E-05,-2.714956E-07,
     &		    1.348204E-02, 1.091504E-05, 6.903760E-07,
     &		    3.106570E-02,-1.665665E-04,-1.590188E-07/

*   CONSTANTS OF THE ARGUMENTS OF THE SHORT-PERIOD PERTURBATIONS
*   BY THE PLANETS:   DCARGS(I,K)
*			I = 1		 I = 2
	DATA DCARGS/ 5.0974222D+00,-7.8604195454652D+02,
     &               3.9584962D+00,-5.7533848094674D+02,
     &		     1.6338070D+00,-1.1506769618935D+03,
     &		     2.5487111D+00,-3.9302097727326D+02,
     &		     4.9255514D+00,-5.8849265665348D+02,
     &               1.3363463D+00,-5.5076098609303D+02,
     &		     1.6072053D+00,-5.2237501616674D+02,
     &		     1.3629480D+00,-1.1790629318198D+03,
     &   	     5.5657014D+00,-1.0977134971135D+03,
     &		     5.0708205D+00,-1.5774000881978D+02,
     &		     3.9318944D+00, 5.2963464780000D+01,
     &		     4.8989497D+00, 3.9809289073258D+01,
     &		     1.3097446D+00, 7.7540959633708D+01,
     &		     3.5147141D+00, 7.9618578146517D+01,
     & 		     3.5413158D+00,-5.4868336758022D+02/

*   AMPLITUDES CCAMPS(N,K) OF THE SHORT-PERIOD PERTURBATIONS
*	      N = 1         N = 2          N = 3         N = 4          N = 5
	DATA CCAMPS/
     & -2.279594E-5, 1.407414E-5, 8.273188E-6, 1.340565E-5,-2.490817E-7,
     & -3.494537E-5, 2.860401E-7, 1.289448E-7, 1.627237E-5,-1.823138E-7,
     &  6.593466E-7, 1.322572E-5, 9.258695E-6,-4.674248E-7,-3.646275E-7,
     &  1.140767E-5,-2.049792E-5,-4.747930E-6,-2.638763E-6,-1.245408E-7,
     &  9.516893E-6,-2.748894E-6,-1.319381E-6,-4.549908E-6,-1.864821E-7,
     &  7.310990E-6,-1.924710E-6,-8.772849E-7,-3.334143E-6,-1.745256E-7,
     & -2.603449E-6, 7.359472E-6, 3.168357E-6, 1.119056E-6,-1.655307E-7,
     & -3.228859E-6, 1.308997E-7, 1.013137E-7, 2.403899E-6,-3.736225E-7,
     &  3.442177E-7, 2.671323E-6, 1.832858E-6,-2.394688E-7,-3.478444E-7,
     &  8.702406E-6,-8.421214E-6,-1.372341E-6,-1.455234E-6,-4.998479E-8,
     & -1.488378E-6,-1.251789E-5, 5.226868E-7,-2.049301E-7, 0.0E0,
     & -8.043059E-6,-2.991300E-6, 1.473654E-7,-3.154542E-7, 0.0E0,
     &  3.699128E-6,-3.316126E-6, 2.901257E-7, 3.407826E-7, 0.0E0,
     &  2.550120E-6,-1.241123E-6, 9.901116E-8, 2.210482E-7, 0.0E0,
     & -6.351059E-7, 2.341650E-6, 1.061492E-6, 2.878231E-7, 0.0E0/

*   CONSTANTS OF THE SECULAR PERTURBATIONS IN LONGITUDE
*   CCSEC3 AND CCSEC(N,K)
*		       N = 1	    N = 2		  N = 3
	DATA CCSEC3/-7.757020E-08/
	DATA CCSEC/ 1.289600E-06, 5.550147E-01, 2.076942E+00,
     &              3.102810E-05, 4.035027E+00, 3.525565E-01,
     &              9.124190E-06, 9.990265E-01, 2.622706E+00,
     &              9.793240E-07, 5.508259E+00, 1.559103E+01/

*   CONSTANTS DCARGM(I,K) OF THE ARGUMENTS OF THE PERTURBATIONS
*   OF THE MOTION OF THE MOON
*			I = 1			I = 2
	DATA DCARGM/  5.1679830D+00, 8.3286911095275D+03,
     &                5.4913150D+00,-7.2140632838100D+03,
     &		      5.9598530D+00, 1.5542754389685D+04/

*     AMPLITUDES  CCAMPM(N,K) OF THE PERTURBATIONS OF THE MOON
*	     N = 1     N = 2	   N = 3           N = 4
	DATA CCAMPM/
     & 1.097594E-01, 2.896773E-07, 5.450474E-02, 1.438491E-07,
     &-2.223581E-02, 5.083103E-08, 1.002548E-02,-2.291823E-08,
     & 1.148966E-02, 5.658888E-08, 8.249439E-03, 4.063015E-08/

*     CCPAMV(K) = A*M*DL/DT (PLANETS), DC1MME = 1-MASS(EARTH+MOON)

	DATA CCPAMV/8.326827E-11, 1.843484E-11, 1.988712E-12,
     &   1.881276E-12/,  DC1MME/0.99999696D0/

*     EXECUTION
*     CONTROL-PARAMETER IDEQ, AND TIME-ARGUMENTS

	IDEQ = NINT(DEQ)
*Test whether original or Modified Julian Date provided.
	IF(DMJDE .LT. 100000.0D0) THEN
	   DJE = DMJDE + 2400000.5D0
	ELSE
	   DJE = DMJDE
	END IF
	DT = (DJE-DCT0)/DCJUL
	T = DT
	DTSQ = DT*DT
	TSQ = DTSQ

*     VALUES OF ALL ELEMENTS FOR THE INSTANT DJE

	DO 100, K = 1,8
	  DLOCAL = MOD(DCFEL(1,K)+DT*DCFEL(2,K)+DTSQ*DCFEL(3,K), DC2PI)
	  IF(K .EQ. 1) THEN
		DML = DLOCAL
	  ELSE
		FORBEL(K-1) = DLOCAL
	  END IF
100 	CONTINUE

	DEPS = MOD(DCEPS(1)+DT*DCEPS(2)+DTSQ*DCEPS(3), DC2PI)

	DO 200, K = 1,17
	  SORBEL(K) = MOD(CCSEL(1,K)+T*CCSEL(2,K)+TSQ*CCSEL(3,K), CC2PI)
200 	CONTINUE

*     SECULAR PERTURABTIONS IN LONGITUDE

	DO 300, K = 1,4
	  A = MOD(CCSEC(2,K)+T*CCSEC(3,K), CC2PI)
	  SN(K) = SIN(A)
300 	CONTINUE

*     PERIODIC PERTURBATIONS OF THE EMB (EARTH-MOON BARYCENTER)

	PERTL = CCSEC(1,1) * SN(1) +CCSEC(1,2)*SN(2) +
     &         (CCSEC(1,3)+T*CCSEC3) * SN(3) + CCSEC(1,4)*SN(4)
	PERTLD = 0.0
	PERTR = 0.0
	PERTRD = 0.0
	DO 400, K = 1,15
	  A = MOD(DCARGS(1,K)+DT*DCARGS(2,K), DC2PI)
	  COSA = COS(A)
	  SINA = SIN(A)
	  PERTL = PERTL + CCAMPS(1,K)*COSA+CCAMPS(2,K)*SINA
	  PERTR = PERTR + CCAMPS(3,K)*COSA+CCAMPS(4,K)*SINA
	  IF (K .LE. 10)  THEN
	    PERTLD = PERTLD+(CCAMPS(2,K)*COSA-CCAMPS(1,K)*SINA)*
     &		CCAMPS(5,K)
	    PERTRD = PERTRD+(CCAMPS(4,K)*COSA-CCAMPS(3,K)*SINA)*
     &		CCAMPS(5,K)
	  END IF
400 	CONTINUE

*     ELLIPTIC PART OF THE MOTION OF THE EMB

	ESQ = E*E
	DPARAM = DC1-ESQ
	PARAM = DPARAM
	TWOE = E+E
	TWOG = G+G
	PHI = TWOE*((1.0-ESQ*0.125) * SIN(G) + E * 0.625 * SIN(TWOG)
     &         + ESQ * 0.5416667 * SIN(G+TWOG))
	F = G+PHI
	SINF = SIN(F)
	COSF = COS(F)
	DPSI = DPARAM/(DC1+E*COSF)
	PHID = TWOE*CCSGD*((1.0+ESQ*1.5)*COSF+E*(1.25-SINF*SINF*0.5))
	PSID = CCSGD*E*SINF/SQRT(PARAM)

*     PERTURBED HELIOCENTIRC MOTION OF THE EMB.

	D1PDRO = (DC1+PERTR)
	DRD    = D1PDRO*(PSID+DPSI*PERTRD)
	DRLD   = D1PDRO*DPSI*(DCSLD+PHID+PERTLD)
	DTL    = MOD(DML+PHI+PERTL, DC2PI)
	DSINLS = SIN(DTL)
	DCOSLS = COS(DTL)
	DXHD   = DRD*DCOSLS-DRLD*DSINLS
	DYHD   = DRD*DSINLS+DRLD*DCOSLS

*    INFLUENCE OF ECCENTRICITY, EVECTION AND VARIATION ON THE
*	GEOCENTRIC MOTION OF THE MOON

	PERTL  = 0.0
	PERTLD = 0.0
	PERTP  = 0.0
	PERTPD = 0.0
	DO 500, K = 1,3
	  A = MOD(DCARGM(1,K)+DT*DCARGM(2,K), DC2PI)
	  SINA = SIN(A)
	  COSA = COS(A)
	  PERTL = PERTL +CCAMPM(1,K)*SINA
  	  PERTLD = PERTLD+CCAMPM(2,K)*COSA
	  PERTP = PERTP +CCAMPM(3,K)*COSA
	  PERTPD = PERTPD-CCAMPM(4,K)*SINA
500	CONTINUE

*    HELIOCENTRIC MOTION OF THE EARTH

	TL    = FORBEL(2)+PERTL
	SINLM = SIN(TL)
	COSLM = COS(TL)
	SIGMA = CCKM/(1.0+PERTP)
	A     = SIGMA*(CCMLD+PERTLD)
	B     = SIGMA*PERTPD
	DXHD  = DXHD+A*SINLM+B*COSLM
	DYHD  = DYHD-A*COSLM+B*SINLM
	DZHD  = -SIGMA*CCFDI*COS(FORBEL(3))

*    BARYCENTRIC MOTION OF THE EARTH

	DXBD = DXHD*DC1MME
	DYBD = DYHD*DC1MME
	DZBD = DZHD*DC1MME
	DO 600, K = 1,4
	  PLON = FORBEL(K+3)
	  POMG = SORBEL(K+1)
	  PECC = SORBEL(K+9)
	  TL   = MOD(PLON+2.0*PECC*SIN(PLON-POMG), CC2PI)
	  SINLP(K) = SIN(TL)
	  COSLP(K) = COS(TL)
 	  DXBD = DXBD+CCPAMV(K)*(SINLP(K)+PECC*SIN(POMG))
	  DYBD = DYBD-CCPAMV(K)*(COSLP(K)+PECC*COS(POMG))
	  DZBD = DZBD-CCPAMV(K)*SORBEL(K+13)*COS(PLON-SORBEL(K+5))
600	CONTINUE

*    TRANSITION TO MEAN EQUATOR OF DATE

	DCOSEP = COS(DEPS)
	DSINEP = SIN(DEPS)
	DYAHD = DCOSEP*DYHD-DSINEP*DZHD
	DZAHD = DSINEP*DYHD+DCOSEP*DZHD
	DYABD = DCOSEP*DYBD-DSINEP*DZBD
	DZABD = DSINEP*DYBD+DCOSEP*DZBD

	IF(IDEQ .EQ. 0) THEN
		DVELH(1) = DXHD
		DVELH(2) = DYAHD
		DVELH(3) = DZAHD
		DVELB(1) = DXBD
		DVELB(2) = DYABD
		DVELB(3) = DZABD
	ELSE

*    GENERAL PRECESSION FROM EPOCH DJE TO DEQ

		DEQDAT = (DJE-DCT0-DCBES)/DCTROP + DC1900
		CALL AX_BARPRE (DEQDAT,DEQ,DPREMA)
		DO 710, N = 1,3
	  	    DVELH(N) = DXHD*DPREMA(N,1)+DYAHD*DPREMA(N,2)+
     &			DZAHD*DPREMA(N,3)
	  	    DVELB(N) = DXBD*DPREMA(N,1)+DYABD*DPREMA(N,2)+
     &			DZABD*DPREMA(N,3)
710		CONTINUE
	END IF
	END
