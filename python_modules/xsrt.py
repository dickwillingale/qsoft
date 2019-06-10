## @package xsrt
# X-ray, electron and proton sequential ray tracing
import xsrtfor 
import numpy as np
# Initialisation and reseting
def init():
    """Initialisation of Fortran common blocks"""
    xsrtfor.qrt_init()
def reset():
    """Reseting of Fortran common blocks to initial condition"""
    xsrtfor.qrt_init()
init()
def rseed(iseed):
    """Set random number seed

    Args:
        iseed:      integer seed
    """
    xsrtfor.qr_rseed(iseed)
def kingrand(alpha,nr):
    """Get samples from a modified Cauchy distribution (King profile)

    Args:
        alpha:       index of modified Cauchy distribution
        nr:          number of samples to return

    Returns:
        | array of random samples from distribution
    """
    ran=np.array([nr])
    return xsrtfor.qrt_kingrand(alpha,ran)
# Proton/electron path tracing
class field: pass
def bfield(dm,pdx,pdy,pdz,ddx,ddy,ddz,px,py,pz):
    """Calculation of magnetic field for array of dipoles

    Args:
        dm:    dipole moments (Gauss cm3)
        pdx:   x positions of dipoles (cm)
        pdy:   y positions of dipoles (cm)
        pdz:   z positions of dipoles (cm)
        ddx:   x direction cosines of dipole moments
        ddy:   y direction cosines of dipole moments
        ddz:   z direction cosines of dipole moments
        px:    x positions for calculation
        py:    y positions for calculation
        pz:    z positions for calculation 

    Returns:
        | list of following
        |  **bf**     magnitude of B-field Gauss
        |  **dx**     direction cosines of B-field in x direction
        |  **dy**     direction cosines of B-field in y direction
        |  **dz**     direction cosines of B-field in z direction
        |  **rmin**   minimum distance from dipoles
    """
    a=xsrtfor.qrt_bfield(dm,pdx,pdy,pdz,ddx,ddy,ddz,px,py,pz)
    b=field()
    dx=a[0]
    dy=a[1]
    dz=a[2]
    b.bf=np.sqrt(dx**2+dy**2+dz**2)
    b.dx=dx/b.bf
    b.dy=dy/b.bf
    b.dz=dz/b.bf
    b.rmin=a[3]
    return b
class eptrace: pass
def prtathena(dm,pdx,pdy,pdz,ddx,ddy,ddz,ekv,tsig,xaper,xdiv,
    rrings,trings,drings,rdet,maxst):
    """Trace proton through Athena telescope with magnetic diverter

    Args:
        dm:     dipole moments (Gauss cm3)
        pdx:    x positions of dipoles (cm)
        pdy:    y positions of dipoles (cm)
        pdz:    z positions of dipoles (cm)
        ddx:    x direction cosines of dipole moments
        ddy:    y direction cosines of dipole moments
        ddz:    z direction cosines of dipole moments
        ekv:    proton energy keV
        tsig:   rms width of beam degrees
        xaper:  x position of mirror aperture cm (focal length)
        xdiv:   x position of diverter input aperture cm
        nrings: number of rings in diverter
        rrings: radius of rings in diverter cm
        trings: radial thickness of rings in diverter cm
        drings: axial depth of rings in diverter cm
        rdet:   radius of detector cm (axial position XDET=0.0)
        maxst:  maximum number of steps along path

    Returns:

        | list of following
        |  **npath**  number of points along path
        |  **iqual**  type of path
        |            0 hits active detector
        |            1 too close to dipole
        |            2 hits telescope tube
        |            3 hits mirror aperture
        |            4 hits diverter
        |            5 hits focal plane beyond detector
        |            6 hit maximum number of steps
        |  **xp**     x positions along path cm
        |  **yp**     y positions along path cm
        |  **zp**     z positions along path cm
    """
    a=xsrtfor.qrt_prtathena(dm,pdx,pdy,pdz,ddx,ddy,ddz,ekv,tsig,xaper,xdiv,
    rrings,trings,drings,rdet,maxst)
    b=eptrace()
    b.npath=a[3]
    b.iqual=a[4]
    b.xp=a[0][:b.npath]
    b.yp=a[1][:b.npath]
    b.zp=a[2][:b.npath]
    return b
def eltmxt(dm,pdx,pdy,pdz,ddx,ddy,ddz,ekv,tsig,xaper,xdiv,
    apy,apz,apsiz,wdiv,ddiv,sdet,maxst):
    """Trace electrons through SVOM MXT telescope with magnetic diverter

    Args:
        dm:     dipole moments (Gauss cm3)
        pdx:    x positions of dipoles (cm)
        pdy:    y positions of dipoles (cm)
        pdz:    z positions of dipoles (cm)
        ddx:    x direction cosines of dipole moments
        ddy:    y direction cosines of dipole moments
        ddz:    z direction cosines of dipole moments
        ekv:    proton energy keV
        tsig:   rms width of beam degrees
        xaper:  x position of mirror aperture cm (focal length)
        xdiv:   x position of diverter input aperture cm
        apy:    y centre of apertures cm
        apz:    z centre of apertures cm
        apsiz:  aperture size (3.8) cm
        wdiv:   extra width of diverter apertures (0.2) cm
        ddiv:   axial depth of diverter cm
        sdet:   size of detector cm (axial position XDET=0.0)
        maxst:  maximum number of steps along path
    
    Returns:
        | list of the following
        |  **npath**  number of points along path
        |  **iqual**  type of path
        |           0 hits active detector
        |           1 too close to dipole
        |           2 hits telescope tube
        |           3 hits mirror aperture
        |           4 hits diverter
        |           5 hits focal plane beyond detector
        |           6 hit maximum number of steps
        |  **xp**     x positions along path cm
        |  **yp**     y positions along path cm
        |  **zp**     z positions along path cm
    """
    a=xsrtfor.qrt_eltmxt(dm,pdx,pdy,pdz,ddx,ddy,ddz,ekv,tsig,xaper,xdiv,
    apy,apz,apsiz,wdiv,ddiv,sdet,maxst)
    b=eptrace()
    b.npath=a[3]
    b.iqual=a[4]
    b.xp=a[0][:b.npath]
    b.yp=a[1][:b.npath]
    b.zp=a[2][:b.npath]
    return b
# Sequential Ray Tracing
class srtdat: pass
def srtreset():
    """Reset Fortran common blocks to intial condition"""
    xsrtfor.qrt_init()
def srtlist():
    """List all current xsrt parameters"""
    xsrtfor.qrt_list()
def shift(iss,pl):
    """Shift position of surface element

    Args:
        iss:      surface element index
        pl:       3 vector shift (dx,dx,dz)
    """
    xsrtfor.qrt_shift(iss,pl)
def rotate(iss,pl,ax,angle):
    """Rotate surface element

    Args:
        iss:      surface element index
        pl:       position of rotation centre (x,y,z)
        ax:       rotation axis (ax,ay,az)
        angle:    rotation angle degrees
    """
    xsrtfor.qrt_rotate(iss,pl,ax,angle)
def kbs(pcen,pnor,raxi,ipack,rmin,rmax,flen,csize,pitch,wall,
    plmin,plmax,idf,iq):
    """Set up a Silicon Kirkpatrick-Baez stack array

    Args:
        pcen:      centre of telescope aperture
        pnor:      normal to aperture
        razi:      reference axis at centre of aperture
        ipack:     packing 0 single module, 1 sunflower, 2 cartesian, 3 wide
                   field cartesian
        rmin:      minimum radius for aperture of constellation
        rmax:      maximum radius for aperture of constellation
        flen:      focal length (-ve for 2nd stack)
        csize:     size of each module in constellation
        pitch:     pitch of K-B slots
        wall:      wall thickness of K-B slots
        plmin:     minimum axial length of K-B slots
        plmax:     maximum axial length of K-B slots
        idf:       deformation index
        iq:        surface quality index

    Returns:
        | **rc**        radius of each module
        | **pc**        azimuth of each module
        | **tc**        rotation of each module
        | **ac**        axial length of each module
        | **nset**      number of module coordinates returned

    Deformation:
        | matrix is (5,nmod) specifying 5 deformations per module.
        |   1 displacement of module vertex in **razi** mm
        |   2 displacement of module vertex in **pnor X razi** mm
        |   3 displaement of module vertex in **pnor** mm
        |   4 rms Gaussian in-plane figure errors radians
        |   5 rms Gaussian out-of-plane figure errors radians
    """
    nmax=2000
    a=xsrtfor.qrt_kbs(pcen,pnor,raxi,ipack,rmin,rmax,flen,csize,
    pitch,wall,plmin,plmax,idf,iq,nmax)
    return a
def sle(pcen,pnor,raxi,flen,cwidth,cheight,pitch,wall,pl,nmir,idf,iq):
    """Defines Schmidt configuration lobster eye stack

    Args:
        pcen:    centre of telescope aperture
        pnor:    normal to aperture
        raxi:    reference axis at centre of aperture
        flen:    focal length (-ve for second stack) measured to front aperture
        cwidth:  cell width = mirror width (its active part)
        cheight: cell height
        pitch:   pitch including mirrors thickness = pdd+wall
        wall:    mirror thickness
        pl:      mirror axial length (depth)
        idf:     deformation index
        iq:      surface quality index

    Returns:
        | **rc**        radius of each module
        | **pc**        azimuth of each module
        | **tc**        rotation of each module
        | **ac**        axial length of each module
        | **nset**      number of module coordinates returned
    """
    nmax=2
    a=xsrtfor.qrt_sle(pcen,pnor,raxi,flen,cwidth,cheight,pitch,wall,
    pl,nmir,idf,iq,nmax)
    return a
def surface(iss,it,ekev,srgh,fmin,pind,alpha,gamma,
    angs,refs,gpitch,dhub,order):
    """Set surface quality parameters

    Args:
        is:     surface quality index
        it:     surface type
                1 refl. (Fresnel), 2 refl. (look-up), 3 refract, 4 diffract
        elev:   photon energy keV
        srgh:   Specific roughness (A**2 mm)
                if -ve then rms figure gradient error radians
        fmin:   Minimum surface spatial frequency
        pind:   Surface roughness power spectrum index
        alpha:  real part of diel. constant or refractive index ratio N1/N2
        gamma:  imaginary part of dielectric constant
        angs:   incidence angles (degrees increasing) (QTYPE 2 and 4)
        refs:   reflectivity values (QTYPE 2 and 4)
        gpitch: d-spacing for grating mm (QTYPE 4)
        dhub:   distance from ruling hub to surface reference point (QTYPE 4)
                if <1.0 mm then d-spacing gradient across ruling (QTYPE 4)
        order:  diffraction order (QTYPE 4)

    The X-ray optical constants **alpha** and **gamma** can be calculated for
    a material of specified composition using the function xscat.xopt().
    It **it=1** the X-ray reflectivity is calculated using the same
    code as in the function xscat.
    """
    nref=len(angs)
    xsrtfor.qrt_surface(iss,it,ekev,srgh,fmin,pind,alpha,gamma,
    nref,angs,refs,gpitch,dhub,order)
def mirror(idd,idf,iq,anml,arfx,apos,alim,nsurf):
    """Set up a plane mirror

    Args:
        idd:      aperture type
        idf:      deformation index
        iq:       surface quality index
        anml:     surface normal
        arfx:     surface reference axis
        apos:     surface reference position
        alim:     limits array (see aperture() )
        nsurf:    number of subsequent surfaces ID=2

    The aperture limits of the mirror are specified in the same way as
    for a top. See function aperture().
    """
    xsrtfor.qrt_mirror(idd,idf,iq,anml,arfx,apos,alim,nsurf)
def moa(pcen,pno,rax,rcur,xyap,pitch,wall,plen,idf,iq):
    """Set up a cylindrical Micro Optic Array

    Args:
        pcen:       centre of MOA
        pno:        normal at centre of MOA
        rax:        cylindrical axis of MOA
        rcur:       radius of curvature of cylinder
        xyap:       half width of MOA aperture
        pitchnput:  pitch of slots
        wall:       wall thickness between slots
        plen:       depth of slots (thickness MOA)
        idf:        deformation index
        iq:         surface quality index
    """
    xsrtfor.qrt_moa(pcen,pno,rax,rcur,xyap,pitch,wall,plen,idf,iq)
def baffle(xmin,xmax,rad,ax,ar,rp,iq):
    """Set up a cylindrical baffle

    Args:
        xmin:      axial position of base
        xmax:      axial position of top
        rad:       radius of cylinder
        ax:        axis direction
        ar:        reference direction perpendicular to axis
        rp:        position of vertex
        iq:        surface quality index

    The axial positions of the base and top, **xmin** and **xmax**
    are local coordinates wrt **ap** along the **ax** direction.
    """
    xsrtfor.qrt_baffle(xmin,xmax,rad,ax,ar,rp,iq)
def lens(idd,idf,iq,anml,arfx,apos,rap,r1,r2,refind,thickq):
    """Set up a lens

    Args:
        id:      lens type 1 spherical, 2 cylindrical
        idf:     deformation index
        iq:      surface quality index
        anml:    surface normal 
        arfx:    surface reference axis
        apos:    surface reference position
        rap:     radius of aperture
        r1,r2:   radii of curvature of lens surfaces
        refind:  refractive index of lens material (or n2/n1)
        thick:   lens thickness

    Surface quality:
        This function sets up 2 surface qualities **iq** and **iq+1** to
        represent the entrance and exit surfaces.
    """
    xsrtfor.qrt_lens(idd,idf,iq,anml,arfx,apos,rap,r1,r2,refind,thickq)
def prism(idd,idf,iq,anml,arfx,apos,rap,d1,d2,refind,thick):
    """Set up a prism

    Args:
        id:     1 small angle, 2 right-angle
        idf:    deformation index
        iq:     quality index
        anml:   entrance surface normal 
        arfx:   reference axis
        apos:   reference position
        rap:    aperture radius
        d1:     small angle radians on entry side (if ID=1)
        d2:     small angle radians on exit side (if ID=1)
        refind: refractive index of material or n2/n1
        thick:  thickness

    Surface quality:
        This function sets up 2 surface qualities **iq** and **iq+1** to
        represent the entrance and exit surfaces.
    """
    xsrtfor.qrt_prism(idd,idf,iq,anml,arfx,apos,rap,d1,d2,refind,thick)
def opgrat(idd,defi,iq,al,fpos,zpos,gpos,alim):
    """Set up a single off-plane grating

    Args:
       id:      1 radial, 2 nested rad., 4 cart., 6 slats
       defi:    deformation index
       iq:      surface quality index
       al:      azimuthal exit angle of zeroth order on cone degrees
       fpos:    position of primary focus (not reflected)
       zpos:    position of zeroth order focus (reflected)
       gpos:    position of centre of grating
       alim:    limits on surface
    
    Returns:
       | **adir**     grating normal
       | **rdir**     grating ruling
       | **dpos**     grating hub
       | **graz**     grating grazing angle radians
       | **gam**      grating cone angle radians
       | **dhub**     grating hub distance

    The grating parameters **gpitch**, **dhub** and **order** and the
    reflectivity as a function of incidence angle are set by
    the surface quality (index **iq**).
    """
    nlim=len(alim)
    a=qrt_opgrat(idd,defi,iq,al,fpos,zpos,gpos,nlim,alim)
    return a
def c1nest(ax,ar,ff,iq,ib,idf,pl,ph,hl,hh,rpl,rph,rhl,rhh,tin,tj,tout):
    """Set up conical approximation to a Wolter type I nest

    Args:
        ax:          optical axis
        ar:          reference axis
        ff:          position of focus
        pl:          array of low axial position of parabola
        ph:          array of high axial position of parabola
        hl:          array of low axial position of hyperbola
        hh:          array of high axial position of hyperbola
        rp:          array of radii parabola near join
        rp:          array of radii parabola at input aperture
        rhl:         array of radii hyperbola at exit aperture
        rjh:         array of radii hyperbola near join
        tin:         array of thicknesses of shells at input aperture
        tj:          array of thicknesses of shells at join plane (or near join)
        tout:        array of thicknesses of shells at the output aperture

    The innermost shell is a dummy and provides an inner stop for the nest

    Deformation:
        A deformation sub-matrix is required for each shell (not
        including the innermost shell)
        Each matrix covers a grid of points along the axis in mm and azimuth
        in radians. The axial range covers both the 1st and 2nd reflection
        surfaces.
        The values in the matrices are radial displacments mm at the grid
        points.
    """
    wrk=np.array([(len(rpl)-1)*2+12])
    xsrtfor.qrt_c1nest(ax,ar,ff,iq,ib,idf,pl,ph,hl,hh,rpl,rph,rhl,rhh,tin,tj,
    tout,wrk)
def w1nest(xj,rj,ra,pl,ph,hl,hh,tin,tj,tout,ax,ar,ff,defi,iq,ib):
    """Set up a Wolter Type I nest

    Args:
        xj:     axial position of join plane
        rj:     array of radii of shells at join
        ra:     ratio of grazing angles
        pl:     low axial position of parabola
        ph:     high axial position of parabola
        hl:     low axial position of hyperbola
        hh:     high axial position of hyperbola
        ti:     array of thicknesses of shells at input aperture
        tj:     array of thicknesses of shells at join plane
        tout:   array of thicknesses of shells at the output aperture
        ax:     direction of axis
        ar:     reference axis in aperture
        ff:     position of focus
        defi:   deformation index
        iq:     reflecting surface quality index
        ib:     back of shells surface quality index

    The innermost shell is a dummy and provides an inner stop for the nest

    Deformation:
        A deformation sub-matrix is required for each shell (not
        including the innermost shell)
        Each matrix covers a grid of points along the axis in mm and azimuth
        in radians. The axial range covers both the 1st and 2nd reflection
        surfaces.
        The values in the matrices are radial displacments mm at the grid
        points.
    """
    xsrtfor.qrt_w1nest(xj,rj,ra,pl,ph,hl,hh,tin,tj,tout,ax,ar,ff,defi,iq,ib)
def wolter2(rp,gp,rh,gh,rm,fovr,ax,ar,ff,idf,iq):
    """Set up Wolter Type II surfaces

    Args:
        rp:         maximum radius of parabola
        gp:         grazing angle (degrees) at maximum radius on parabola
        rh:         maximum radius of hyperbola
        gh:         grazing angle (degrees) at maximum radius on hyperbola
        rm:         minimum radius of parabola
        fovr:       radius of field of view degrees
        ax:         direction of axis of telescope
        ar:         reference direction perpendicular to axis
        ff:         position of focus of telescope
        id:         deformation index
        iq:         surface quality index

    Deformation:
        A single deformation matrix specifies radial displacement errors mm
        over a grid of points in  axial mm and azimuthal radians.
        The grid covers the axial positions along both the 1st and 2nd
        reflection surfaces.
    """
    xsrtfor.qrt_wolter2(rp,gp,rh,gh,rm,fovr,ax,ar,ff,idf,iq)
def spider(cone,apos,anml,arfx,nsec,cwid,awid):
    """Set up a support spider

    Args:
        cone:       90-half cone angle degrees (0.0 for plane)
        apos:       axial position of vertex (centre)
        anml:       direction of normal to aperture (optic axis)
        afrx:       reference axis in aperture
        nsec:       number of sectors (number of arms)
        cwid:       constant arm width
        awid:       angular arm width degrees

    The surface figure of the spider aperture is a cone.
    The width of the arms is given by **cwid+aw*R**
    where **R** is the radial distance
    from the axis of the cone and **aw** is **awid**
    degrees converted into radians.
    """
    xsrtfor.qrt_spider(cone,apos,anml,arfx,nsec,cwid,awid)
def trace(ideb,riris,iopt):
    """Perform ray tracing

    Args:
        ideb:       debugging level (0 none)
        riris:      radius about centre of detector for analysis
                    if 0.0 then no analysis of detected distribution
        iopt:       controls form of output

    **iopt**:
        | -2 save traced.dat and detected.dat files
        | -1 save detected.dat
        | 0 don't save files or adjust focus
        | 1 adjust focus and save detected.dat
        | 2 adjust focus and save detected.dat and traced
        | Only rays with **iopt** reflections are used in adjustment

    Returns:
        | list of the following
        |  **area**    detected area within RIRIS
        |  **dshft**   axial shift to optimum focus (0.0 if IOPT<=0)
        |  **ybar**    y centroid of detected distribution
        |  **zbar**    z centroid of detected distribution
        |  **rms**     rms radius of detected distribution
    """
    a=xsrtfor.qrt_trace(ideb,riris,iopt)
    b=srtdat()
    b.area=a[0]
    b.dshft=a[1]
    b.ybar=a[2]
    b.zbar=a[3]
    b.rms=a[4]
    return b
def source(it,sd,sp,ap,an,ar,al,apry,nr,idef):
    """Set up source 

    Args:
        it:        source type 
        sd:        source direction cosines
        sp:        source position
        ap:        aperture position
        an:        aperture normal
        ar:        aperture reference axis
        al:        aperture limits
        apry:      area per ray - if 0 then use NR
        nr:        number of rays
        idef:      deformation index
    **it** values:
      | 1 point source at infinity,        radial limits
      | 2 point source at infinity,        cartesian limits
      | 3 point source at finite distance, radial limits
      | 4 point source at finite distance, cartesian limits

    Deformation:
      The deformation index is used to specify a pixel array which spreads
      out a point source into an angular distribution. The deformation data are
      set using two functions xsrt.deform() and xsrt.defmat().
      If the source is at infinity the x and y sample arrays must be in radians
      measured from the direction **sd** along the reference axis **ar** and the
      other axis.  (**an** cross **ar**).
      If the source is at a finite distance then x and y are displacements
      in mm (or whatever distance unit is used) of the position **sp** along
      the reference axis **ar** and the other axis (**an** cross **ar**).
    """
    xsrtfor.qrt_source(it,sd,sp,ap,an,ar,al,apry,nr,idef)
def detector(idd,dpos,dnml,drfx,dlim,radet):
    """Set up detector

    Args:
        id:       detector type
        dpos:     detector position
        dnml:     detector normal
        drfx:     detector reference axis
        dlim:     detector limits
        radet:    radius of curvature of spherical detector

    **id** values:
      | 1 planar detector, radial limits
      | 2 planar detector, cartesian limits
      | 3 spherical detector, radial limits
      | 4 spherical detector, cartesian limits
    """
    xsrtfor.qrt_detector(idd,dpos,dnml,drfx,dlim,radet)
def fresnel(nreal,kimag,angs):
    """Calculate reflectivity using Fresnel's equations

    Args:
        nreal:      real part of refractive index
        kimag:      imaginary part of refractive index
        angs:       array of incidence angles (degrees range 0-90)
    
    Returns:
        | list of following
        |  **rs**          sigma reflectivity
        |  **rp**          pi reflectivity
        |  **runp**        unpolarized reflectivity

    | Reference "Handbook of Optical Constants of Solids" Ed. Edward D.Palik
    | Academic Press 1985, page 70
    | If angs(I) out of range 0-90 degrees returns zero reflectivities
    """
    a=xsrtfor.qrt_fresnel(nreal,kimag,angs)
    b=srtdat()
    b.rs=a[0]
    b.rp=a[1]
    b.runp=a[2]
    return b
def deform(idd,it,nm,nx,ny):
    """Set up surface deformation dimensions

    Args:
        id:      deformation index
        it:      deformation type (1 matrix, 2 sin() in x and/or y)
        nm:      number of sub-matrices
        nx:      number of x samples
        ny:      number of y samples

    Note when **it=2** the sin() function requires 3 parameters, amplitude (**amp**),
    wavelength (**lam**) and phase (**phi**) all in mm (or the same units)
    **defx=amp.sin(2.pi.(x-phi)/lam)**

    If **it** -ve then uses **abs(it)** but expects a 2nd deformation to be
    defined which will be added to the current deformation.
    """
    xsrtfor.qrt_defs(idd,it,nm,nx,ny)
def defmat(idd,im,xsam,ysam,zdef):
    """Load deformation matrix

    Args:
        id:      deformation index
        im:      sub-matrix index (runs from 1-N)
        xsam:    x values
        ysam:    y values
        zdef:    deformation matrix

    Note for deformation of a surface of revolution X is axial direction and
    Y is azimuth.

    Used when **it=1** in deform().
    """
    nw=len(xsam)*len(ysam)
    w1=np.empty(nw,dtype=float)
    w2=np.empty(nw,dtype=float)
    xsrtfor.qrt_mats(idd,im,xsam,ysam,zdef,w1,w2)
def defparxy(idd,im,xpar,ypar):
    """Sets parameters of deformation defined by parameters in x and y axes.

    Args:
        idd:  deformation index
        im:   sub-surface index
        xpar: parameters of deformations in x axis
        ypar: parameters of deformations in y axis

    Empty field [] can be used as xpar or ypar if there are no deformations
    in the corresponding axis.

    Use this routine to set parameteric deformations, e.g.
    when **it=2** in deform() the sin() function requires 3 parameters,
    amplitude (**amp**), wavelength (**lam**) and phase (**phi**) all in mm
    (or the same units) **defx=amp.sin(2.pi.(x-phi)/lam)**
    """
    xsrtfor.qrt_defparxy(idd,im,xpar,ypar)
def sqpore(pcen,pnorm,raxis,rcur,ipack,rap,pitch,wall,plen,idf,
    iq,plmin,plmax,fibre,ar):
    """Set up slumped square pore MPO

    Args:
        pcen:     position of centre of plate
        rnorm:    normal at centre of plate
        raxis:    reference axi at centre of plate
        rcur:     radius of curvature
        ipack:    pore packing
                    1 cart, 2 rad, 3 waff, 4 octag, 5 rand, 6 MIXS, 7 NFL
        rap:      half width of plate aperture
        pitch:    pitch of pores on a cartesian grid
        wall:     pore wall thickness
        plen:     length of pores
        idf:      deformation index
        iq:       surface quality index
        plmin:    minimum pore length
        plmax:    maximum pore length
        fibre:    size of fibre bundle in packing
        ar:       array of additional parameters specifying plate apertures

    Additional parameter array:
      | used for ipack 6 or 7
      | RC,PC  radius and azimuth of each plate
      | TC     rotation angle of plate
      | WC     width of each plate x
      | HC     height of each plate y 
      | AC     axial length of each plate
      | CC     spare
      | GC     spare
    """
    xsrtfor.qrt_sqpore(pcen,pnorm,raxis,rcur,ipack,rap,pitch,wall,plen,
    idf,iq,plmin,plmax,fibre,ar)
def sqmpoarr(pcen,pnorm,raxis,rcur,hwid,idf,ar):
    """Set up an array of square pore MPOs

    Args:
        pcen:     position of centre of array
        rnorm:    normal at centre of array
        raxis:    reference axi at centre of array
        rcur:     radius of curvature of array
        hwid:     half width of array aperture
        idf:      deformation index for array
        ar:       array of additional parameters for each MPO in array

    Additional parameter array:
      | there are 15 parameters for each MPO
      | XP     x position of each MPO
      | YP     y position of each MPO
      | TC     rotation angle of each MPO
      | WC     width of each MPO delx 
      | HC     height of each MPO dely 
      | AC     axial length of each MPO (thickness)
      | CU     radius of curvature of each MPO
      | MF     size of multifibres in MPO
      | PP     pitch of pores (pore size + wall thickness)
      | WA     wall thickness between pores
      | SQ     reflecting surface quality index for pores
      | BU     bias angle x radians
      | BV     bias angle y radians
      | BZ     MPO efficiency
      | SP     MPO tilt error index

    Deformation:
      | matrix (5,nmpo) specifying 5 deformations for each MPO
      |    1 scaling of intrinsic slumping distortion - 1.0 for model
      |    2 +ve amplitude thermoelastic pore axial pointing errors radians
      |      -ve amplitude fixed pattern tilt errors radians
      |    3 Gaussian rms pore rotation errors (about pore axis) radians
      |    4 pore shear error amplitude mm within multifibre structure 
      |    5 +ve Gaussian rms pore tilt/figure error radians (x2.37=FWHM)
      |      -ve Modified Cauchy pore tilt/figure errors radians (x2=FWHM)
      |      both applied independently in 2 tilt axes
      |      When tilt error index 1 then Lorentzian, otherwise King profile
    """
    xsrtfor.qrt_sqmpoarr(pcen,pnorm,raxis,rcur,hwid,idf,ar)
def sipore(pcen,pnorm,raxis,flen,rpitch,apitch,wall,
    rm,pm,tm,wm,hm,am,cm,gm,wfr,a2j,idf,iq):
    """Set up Silicon Pore Optics

    Args:
        pcen:      position of centre of aperture (above join plane)
        rnorm:     normal at centre of aperture
        raxis:     reference axis at centre of aperture
        flen:      focal length
        rpitch:    pore radial pitch
        apitch:    pore azimuthal pitch
        wall:      wall thickness
        rm:        array of module radii
        pm:        array of module azimuths (radians)
        tm:        array of module rotations (radians normally 0)
        wm:        array of module widths (radial mm)
        hm:        array of module heights (azimuthal mm)
        am:        array of module lengths (axial pore length mm)
        cm:        array of module curvature signatures
                   0 conical, 1 Wolter, 2 curve-plane, 3 constant...
        gm:        array of module grazing angle ratios
        wrf:       module frame width (surrounding module aperture)
        a2j:       aperture to join plane axial distance
        idf:       deformation index
        iq:        reflecting surface quality
    """
    xsrtfor.qrt_sipore(pcen,pnorm,raxis,flen,rpitch,apitch,wall,
    rm,pm,tm,wm,hm,am,cm,gm,wfr,a2j,idf,iq)
def spoarr(pcen,pnorm,raxis,flen,a2j,
    rm,pm,tm,wm,hm,am,cm,gm,rpitch,wall,apitch,wfr,siq,biq,idf):
    """Set up Silicon Pore Optics Wolter I array

    Args:
        pcen:      position of centre of aperture (above join plane)
        rnorm:     normal at centre of aperture
        raxis:     reference axis at centre of aperture
        flen:      focal length
        a2j:       aperture to join plane axial distance
        rm:        array of module radii
        pm:        array of module azimuths (radians)
        tm:        array of module rotations (radians normally 0)
        wm:        array of module widths (radial mm)
        hm:        array of module heights (azimuthal mm)
        am:        array of module lengths (axial pore length mm)
        cm:        array of module curvature signatures
                   0 conical, 1 Wolter, 2 curve-plane, 3 constant...
        gm:        array of module grazing angle ratios
        rpitch:    array of module pore radial pitch
        wall:      array of module wall thickness
        apitch:    array of module pore azimuthal pitch
        rwi:       array of module pore rib thickness
        wrf:       array of module frame widths (surrounding module aperture)
        siq:       array of surface quality indices
        biq:       array of non-reflecting surface quality indices
        idf:       deformation index

    Deformation:
        | matrix is (6,nmod) specifying 6 deformations per module.
        |   1 displacement of module in direction **raxis** mm
        |   2 displacement of module in direction of **rnorm X raxis** mm
        |   3 error in focal length of module
        |   4 Gaussian rms in-plane figure error radians
        |   5 Gaussian rms out-of-plane figure error radians
        |   6 width mm along axial edges of module in which figure degrades
    """
    xsrtfor.qrt_spoarr(pcen,pnorm,raxis,flen,a2j,
    rm,pm,tm,wm,hm,am,cm,gm,rpitch,wall,apitch,rwi,wfr,siq,biq,idf)
def aperture(idd,idf,ap,an,ar,alim,nsurf):
    """Set up an aperture stop

    Args:
        id:    aperture type
        idf:   deformation index
        ap:    position of aperture
        an:    normal to aperture plane
        ar:    reference axis in aperture plane
        alim:  limit values (depend on id see above)
        nsurf: number of subsequent surfaces per aperture (id=2)

    **id** values:
      |   1 Single annulus, radial limits (aref,rmin,rmax,0,0,0)
      |   2 Nested annuli, radial limits  (aref,rmin1,rmax1,rmin2,0,0)
      |     aref is an axial reference position used for deformation
      |   3 Hole Cartesian limits (xmin,ymin,xmax,ymax,0,0)
      |   4 Block Cartesian limits (xmin,ymin,xmax,ymax,0,0)
      |   5 Cartesian grid limits (pitchx pitchy ribx riby,0,0)
      |   6 Radial/azimuthal sector limits (rmin,rmax,amin,amax,0,0)
      |     amin and amax in radians range 0-2pi
      |   7 Parallelogram limits (xmin,ymin,xmax,ymax,dx,0)
      |     where dx is the shear in X over the distance ymax-ymin
      |   8 Aperture for MCO test station limits (hsize dols,0,0,0,0)

    The parameter **nsurf** is used for radially nested apertures so that
    the code knows how many surface elements to skip if a ray penitrates
    a particular annulus.

    If the limits are radial the deformation is defined as a vector (1-d
    sampless in azimuth) and is applied in the radial direction.
    
    If the limits are cartesian, radial/azimuthal or parallelogram
    the deformation is defined over a matrix (2-d) and is applied in the
    direction of the normal.

    No deformation is used for **id=8**.
    """
    xsrtfor.qrt_aperture(idd,idf,ap,an,ar,alim,nsurf)
def aparray(ap,an,ar,xhw,yhw):
    """Set up an array of rectangular apertures

    Args:
        ap:    positions of aperture centres
        an:    normals to aperture planes
        ar:    reference axes in aperture planes
        xhw:   half widths in reference axes
        yhw:   half widths in other axes

    The apertures are tested in sequence until a ray is found to go through.
    Control then jumps to the next surface after the apertures.
    """
    apt=np.transpose(ap)
    ant=np.transpose(an)
    art=np.transpose(ar)
    xsrtfor.qrt_aparray(apt,ant,art,xhw,yhw)
def bafflearray(ap,an,ar,rmin,rmax,amin,amax):
    """Set up an array of baffle elements

    Args:
        ap:    positions of aperture centres
        an:    normals to aperture planes
        ar:    reference axes in aperture planes
        rmin:  minimum radius of sectors
        rmax:  maximum radius of sectors
        amin:  minimum azimuth of sectors (radians)
        amax:  maximum azimuth of sectors (radians)

    The sector elements are tested in sequence until the ray is found to hit
    an elements or miss all elements
    Control then jumps to the next surface after the apertures.
    """
    apt=np.transpose(ap)
    ant=np.transpose(an)
    art=np.transpose(ar)
    xsrtfor.qrt_bafflearray(apt,ant,art,rmin,rmax,amin,amax)
def elips(org,axs,cen,xmin,xmax,amin,amax,smb,rab,ide,iq):
    """Set up elliptical grazing indidence mirror

    Args:
        org:       local origin on surface of ellipse
        axs:       reference axis in aperture
        cen:       focus of ellipse
        xmin:      minimum axial position
        xmax:      maximum axial position
        amin:      minimum azimuth radians
        amax:      maximum azimuth radians
        smb:       conic coefficient
        rab:       conic coefficient
                   conic equation of form r**2=rab**2.x**2+smb**2
        ide:       deformation index
        iq:        reflecting surface quality
    """
    xsrtfor.qrt_elips(org,axs,cen,xmin,xmax,amin,amax,smb,rab,ide,iq)
