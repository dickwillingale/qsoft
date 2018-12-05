## @package xscat
# X-ray Scattering Code in Qsoft
from __future__ import print_function
import xscatfor
# Initialisation use Fortran routine directly
xscatfor.qrs_init()
# Shareble object version of Wiscombe Mie scattering subroutine
# and Rayleigh-Gans approximation for scattering by dust
class mievs:
    """Mie scattering and Rayleigh-Gans approximation for scattering by dust
    """
    def list(self):
        """List object"""
        print("qext",self.qext)
        print("qsca",self.qsca)
        print("gqsc",self.gqsc)
        print("pmom",self.pmom)
        print("sforw",self.sforw)
        print("sback",self.sback)
        print("s1",self.s1)
        print("s2",self.s2)
        print("tforw",self.tforw)
        print("tback",self.tback)
        print("spike",self.spike)
def miev0(xx,crefin,pfct,mimc,anya,xmu,nmom,ipolzn,momd,prt):
    """`Wiscombe subroutine <https://github.com/cfinch/Mie_scattering/blob/master/Wiscombe/MIEV0.f>`_

    Computes Mie scattering and extinction efficiencies; asymmetry
    factor;  forward- and backscatter amplitude;  scattering
    amplitudes vs. scattering angle for incident polarization parallel
    and perpendicular to the plane of scattering;
    coefficients in the Legendre polynomial expansions of either the
    unpolarized phase function or the polarized phase matrix;
    some quantities needed in polarized radiative transfer;  and
    information about whether or not a resonance has been hit.

    Input and output variables are described in file MIEV.doc.
    Many statements are accompanied by comments referring to
    references in MIEV.doc, notably the NCAR Mie report which is now
    available electronically and which is referred to using the
    shorthand (Rn), meaning Eq. (n) of the report.

    Args:
        xx:     Mie size parameter ( 2 * pi * radius / wavelength )
        crefin: Complex refractive index ( imag part can be + or -,
                but internally a negative imaginary index is assumed ).
                If imag part is - ,  scattering amplitudes as in Van
                de Hulst are returned;  if imag part is + , complex
                conjugates of those scattering amplitudes are returned
                (the latter is the convention in physics).
                ** NOTE ** In the 'PERFECT' case, scattering amplitudes
                in the Van de Hulst (Ref. 6 above) convention will
                automatically be returned unless  Im(CREFIN)  is
                positive;  otherwise, CREFIN plays no role.
        pfct:   TRUE, assume refractive index is infinite and use
                special case formulas for Mie coefficients  'a'
                and  'b'  ( see Kerker, M., The Scattering of
                Light and Other Electromagnetic Radiation, p. 90 ).
                This is sometimes called the 'totally reflecting',
                sometimes the 'perfectly conducting' case.
                ( see CREFIN for additional information )
        mimc:   (positive) value below which imaginary refractive
                index is regarded as zero (computation proceeds
                faster for zero imaginary index)
        anya:   TRUE, any angles whatsoever may be input through
                XMU.  FALSE, the angles are monotone increasing
                and mirror symmetric about 90 degrees (this option
                is advantageous because the scattering amplitudes
                S1,S2 for the angles between 90 and 180 degrees
                are evaluable from symmetry relations, and hence
                are obtained with little added computational cost.)
        numang: No. of angles at which scattering amplitudes
                S1,S2 are to be evaluated  ( set = 0 to skip
                calculation of S1,S2 ).  Make sure NUMANG does
                not exceed the parameter MAXANG in the program.
        xmu(N): Cosines of angles ( N = 1 TO NUMANG ) at which S1,S2
                are to be evaluated.  If ANYANG = FALSE, then
                (a) the angles must be monotone increasing and
                mirror symmetric about 90 degrees (if 90-A is
                an angle, then 90+A must be also)
                (b) if NUMANG is odd, 90 degrees must be among
                the angles
        nmom:   Highest Legendre moment PMOM to calculate,
                numbering from zero ( NMOM = 0 prevents
                calculation of PMOM )
        ipolzn: POSITIVE, Compute Legendre moments PMOM for the
                Mueller matrix elements determined by the
                digits of IPOLZN, with 1 referring to M1,
                2 to M2, 3 to S21, and 4 to D21 (Ref. 3).
                E.g., if IPOLZN = 14 then only moments for
                M1 and D21 will be returned.
                0,  Compute Legendre moments PMOM for the
                npolarized unnormalized phase function.
                NEGATIVE, Compute Legendre moments PMOM for the
                Sekera phase quantities determined by the
                digits of ABS(IPOLZN), with 1 referring to
                R1, 2 to R2, 3 to R3, and 4 to R4 (REF. 4).
                E.g., if IPOLZN = -14 then only moments for
                R1 and R4 will be returned.
                ( NOT USED IF  NMOM = 0 )
        momd: Determines first dimension of PMOM, which is dimensioned
                internally as PMOM( 0:MOMDIM, * ) (second dimension must
                be the larger of unity and the highest digit in
                IPOLZN; if not, serious errors will occur).
                Must be given a value, even if  NMOM = 0.  Minimum: 1.
        prt(L): Print flags (LOGICAL).  L = 1  prints  S1,S2, their
                squared absolute values, and degree of polarization,
                provided NUMANG is non-zero.   L = 2  prints all
                output variables other than  S1,S2.

    Returns:
        list containing the following

    | **qext**:       (REAL) extinction efficiency factor  ( Ref. 2, Eq. 1A )
    | **qsca**:       (REAL) scattering efficiency factor  ( Ref. 2, Eq. 1B )
    | **gqsc**:       (REAL) asymmetry factor times scattering efficiency
    |  ( Ref. 2, Eq. 1C )  ( allows calculation of radiation
    |   pressure efficiency factor  QPR = QEXT - GQSC )
    |  NOTE -- S1, S2, SFORW, SBACK, TFORW, AND TBACK are calculated
    |  internally for negative imaginary refractive index;
    |  for positive imaginary index, their complex conjugates
    |  are taken before they are returned, to correspond to
    |  customary usage in some parts of physics ( in particular,
    |  in papers on CAM approximations to Mie theory ).
    | **pmom(M,NP)**: (REAL) moments  M = 0 to NMOM  of unnormalized NP-th
    |  phase quantity PQ  ( moments with  M .GT. 2*NTRM  are
    |  zero, where  NTRM = no. terms in Mie series =
    |  XX + 4*XX**1/3 + 1 ) 
    | **PQ( MU, NP )**: = sum( M=0 to infinity ) ( (2M+1)
    |  * PMOM( M,NP ) * P-sub-M( MU ) )
    |  WHERE  MU = COS( scattering angle )
    |  P-sub-M = M-th Legendre polynomial 
    |  and the definition of 'PQ' is as follows:
    |  IPOLZN.GT.0:  PQ(MU,1) = CABS( S1(MU) )**2
    |  PQ(MU,2) = CABS( S2(MU) )**2
    |  PQ(MU,3) = RE( S1(MU)*CONJG( S2(MU) ) )
    |  PQ(MU,4) = - IM( S1(MU)*CONJG( S2(MU) ) )
    |  ( called M1, M2, S21, D21 in literature )
    |  IPOLZN=0:  PQ(MU,1) = ( CABS(S1)**2 + CABS(S2)**2 ) / 2
    |  ( the unnormalized phase function )
    |  IPOLZN.LT.0:  PQ(MU,1) = CABS( T1(MU) )**2
    |  PQ(MU,2) = CABS( T2(MU) )**2
    |  PQ(MU,3) = RE( T1(MU)*CONJG( T2(MU) ) )
    |  PQ(MU,4) = - IM( T1(MU)*CONJG( T2(MU) ) )
    |  ( called R1, R2, R3, R4 in literature )
    |  The sign of the 4th phase quantity is a source of
    |  confusion.  It flips if the complex conjugates of
    |  S1,S2  or  T1,T2  are used, as occurs when a
    |  refractive index with positive imaginary part is
    |  used (see discussion below).  The definition above
    |  is consistent with a negative imaginary part.
    |  See Ref. 5 for correct formulae for PMOM ( Eqs. 2-5
    |  of Ref. 3 contain typographical errors ).  Ref. 5 also
    |  contains numerous improvements to the Ref. 3 formulas.
    |  NOTE THAT OUR DEFINITION OF MOMENTS DIFFERS FROM REF. 3
    |  in that we divide out the factor (2M+1) and number the
    |  moments from zero instead of one.
    |  ** WARNING **  Make sure the second dimension of PMOM
    |  in the calling program is at least as large as the
    |  absolute value of IPOLZN.
    |  For small enough values of XX, or large enough values
    |  of M,  PMOM  will tend to underflow.  Thus, it is
    |  unwise to assume the values returned are non-zero and,
    |  for example, to divide some quantity by them.
    | **s1(N),s2(N)**:  (COMPLEX) Mie scattering amplitudes at angles specified
    |   by XMU(N) ( N=1 to NUMANG )  ( Ref. 2, Eqs. 1d-e ).
    | **sforw**:  (COMPLEX) forward-scattering amplitude S1 at
    |   0 degrees.  ( S2(0 deg) = S1(0 deg) )
    | **sback**:  (COMPLEX) backscattering amplitude S1 at
    |   180 degrees.   ( S2(180 deg) = - S1(180 deg) )
    | **tforw(I)**:   (COMPLEX) values of
    |   I=1:  T1 = ( S2 - (MU)*S1 ) / ( 1 - MU**2 )
    |   I=2:  T2 = ( S1 - (MU)*S2 ) / ( 1 - MU**2 )
    |   At angle theta = 0 ( MU = COS(theta) = 1 ), where the
    |   expressions on the right-hand side are indeterminate.
    |   ( these quantities are required for doing polarized
    |   radiative transfer (Ref. 4, Appendix). )
    | **tback(I)**:   (COMPLEX) values of  T1 (for I=1) or  T2 (for I=2) at
    |   angle  theta = 180 degrees ( MU = COS(theta) = - 1 ).
    | **spike**:     (REAL) magnitude of the smallest denominator of
    |   either Mie coefficient (a-sub-n or b-sub-n),
    |   taken over all terms in the Mie series past
    |   N = size parameter XX.  Values of SPIKE below
    |   about 0.3 signify a ripple spike, since these
    |   spikes are produced by abnormally small denominators
    |   in the Mie coefficients (normal denominators are of
    |   order unity or higher).  Defaults to 1.0 when not
    |   on a spike.  Does not identify all resonances
    |   (we are still working on that).
    """
    a=xscatfor.miev0(xx,crefin,pfct,mimc,anya,xmu,nmom,ipolzn,momd,prt)
    b=mievs()
    b.qext=a[0]
    b.qsca=a[1]
    b.gqsc=a[2]
    b.pmom=a[3]
    b.sforw=a[4]
    b.sback=a[5]
    b.s1=a[6]
    b.s2=a[7]
    b.tforw=a[8]
    b.tback=a[9]
    b.spike=a[10]
    return b
#print(xrayfor.miev0.__doc__)
# X-ray scattering from interstellar dust
def dustrings(data,derr,dsou,ekv,srate,dts,td,zd,amin,amax,qa,sig1):
    """Model fitting to dust X-ray scattering halo rings

    Args:
        data:      data array of rings surface brightness cts/s/str
        derr:      array of errors on data
        dsou:      distance to source PC
        ekv:       array of energies keV (equally spaced across observed band)
        srate:     source spectrum cts/s/keV
        dts:       source burst duration
        td:        delay time of observations secs
        zd:        fraction of source distance to rings
        amin,amax: grain size radius range microns
        qa:        grain size distribution index N(a)=A.a^-qa
        sig1:      differential cross-section of 1 grain cm2, 1 keV, 0.1 microns

    Returns:
        the following

    | **angs**:       array of angles
    | **ndust**:      N dust columns cm-2
    | **edust**:      errors on N dust columns cm-2
    | **model**:      model cts/s/str in rings
    | **chisq**:      Chi-Squared between data and model
    | **ndof**:       ndof
    """
    return xscatfor.qra_dustrings(data,derr,dsou,ekv,srate,dts,td,zd,
    amin,amax,qa,sig1)
class duststat:
    """Statistic of fit of model to dust rings"""
    def __init__(self,data,derr,dsou,ekv,srate,dts,td,zd,sig1):
        self.data=data
        self.derr=derr
        self.dsou=dsou
        self.ekv=ekv
        self.srate=srate
        self.dts=dts
        self.td=td
        self.zd=zd
        self.sig1=sig1
        self.ncall=0
    def __call__(self,pars):
        """input parameters as for dustrings()

        Returns:
            chi-squared value
        """
        angs,ndust,edust,model,chisq,ndof=dustrings(self.data,self.derr,
        self.dsou, self.ekv,self.srate,self.dts,self.td,self.zd,
        pars[0]*pars[1],pars[1],pars[2],self.sig1)
        self.ncall=self.ncall+1
        print(self.ncall,pars,chisq)
        return chisq
def dustthetascat(td,ds,zd):
    """ Scattering angle radians

    Args:
        td:    delay time seconds after source flare (assumed delta function)
        ds:    distance to source parsecs (convert to m - 3.086e16/parsec)
        zd:    fractional distance to dust 

    Returns:
        angle in radians
    """
    return xscatfor.qra_dustthetascat(td,ds,zd)
# X-ray optical properties
class xoptdat: pass
def xopt(mspec,rho,ekev,itype):
    """X-ray optical properties of a material

    Args:
        mspec:        specification of composition
        rho:          density gm/cm**3
        ekev:         array of photon energies in keV
        itype:        data source (0=Cromer, 1=Henke)

    Returns:
        list with the following

    | **alpha**:     array of real parts dielectric constant
    | **gamm**:      array of imaginary parts dielectric constant
    | **absl**:      array of absortion lengths cm-1
    | **f1**:        array of real parts scattering factor
    | **f2**:        array of imaginary parts scattering factor
    """
    a=xscatfor.qrt_xopt(len(mspec),mspec,rho,ekev,itype)
    b=xoptdat() 
    b.alpha=a[0]
    b.gamma=a[1]
    b.absl=a[2]
    b.f1=a[3]
    b.f2=a[4]
    return b
def xfresnel(alpha,gamma,angs):
    """Calculate X-ray reflectivity using Fresnel's equations

    Args:
        alpha:        real incremental part dielectric constant
        gamma:        imaginary part of dielectric constant
        angs:         incidence angles (degrees)

    Returns:
        list of following:

    | **rs**:       sigma reflectivity
    | **rp**:       pi reflectivity
    | **runp**:     unpolarized reflectivity

    If angs(I) out of range 0-90 degrees returns zero reflectivity
    """
    a=xscatfor.qrt_xfresnel(alpha,gamma,angs)
    b=xoptdat()
    b.rs=a[0]
    b.rp=a[1]
    b.runp=a[2]
    return b
def mlayer(angs,ekev,nr,ni,d,nper):
    a=xscatfor.qrt_mlayer(angs,ekev,nr,ni,d,nper)
    b=xoptdat()
    b.rsig=a[0]
    b.rpi=a[1]
    b.tsig=a[2]
    b.tpi=a[3]
    return b
