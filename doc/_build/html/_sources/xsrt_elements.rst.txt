Optical Elements and Coordinates
********************************

Each of the elements, including the source and detector, are specified
by:

* 3 surface reference vectors - origin position, surface normal at origin and reference tangent at origin
* The surface figure - planar, spherical or conic section - parameters to define the curvature etc. of the figure
* The surface boundary - circles or rectangles in local surface coordinates
* A surface quality - source of rays, detection, reflection, diffraction, scattering, refraction, absorption
* The surface deformation - a grid of displacements defined in local surface coordinates

A full list of all the currently defined elements is produced by the 
function srtlist().

Individual elements referenced by the surface element index
can be shifted and rotated using shift() and rotate().

The data base of elements can be cleared to the initial condition (no
elements defined) using the function reset(). If elements are
repeatedly defined within a procedure (for instance within a loop)
the safe and prefered option is to reset()
everything and redefine all elements each time they are required.

**Coordinates**

There is no fixed coordinate system and elements can be set/defined at
any orientation. However it is conventional to use the X-axis as the
optical axis (which defines the direction of paraxial rays) and the
Y-axis as the nominal tangent reference axis. Of course given elements
may not be aligned exactly with the X-axis and Y-axis.
In most cases the local coordinates in the detector plane are nominally
aligned with the Y-axis and Z-axis. Rays are usually traced from right
to left travelling in the -X direction but this is not necessary and
it is possible for rays to bounce back and forth as in, for example,
a cassegrain system.

The source is always the first element in the sequence. All other
elements are placed in sequence as they are defined. If the source()
function is used repeatedly the source specification will be overwritten
each time.
If the detector command is used repeatedly a new detector will be
added to the sequence each time and all detectors defined will be active.

Local surface coordinates
are specified using the tangent plane
to the surface at the point defined as the surface origin. For a sphere
points on this tangent plane are projected
onto the surface along the normal to the surface at the surface
origin (Lambert's projection). The local x-axis is specified by a tangent vector
at the surface origin. The local y-axis is the cross product of the
normal and tangent vector at the surface origin.

The coordinates used for the limits of apertures and stops are given
in the docstrings of the xsrt.aperture() function.

The local coordinates used for surfaces of revolution generated from conic
sections (hyperbola, parabola, ellipse) depend on whether the surface
is designated as being "normal" or "grazing" incidence. For normal
incidence they are defined in a similar way to the planar or spherical
surfaces as given above. For grazing incidence a cylindrical coordinate
system is used where the axis is the normal to the surface at
the surface origin and the azimuth is the rotation about this axis
with zero at the surface reference axis at the surface origin. Local
coordinates are given as axial position and azimuthal position (radians).
The limits of such surfaces are specified by axial and/or radial limits
corresponding to the bottom and top edges of the surface of revolution.
