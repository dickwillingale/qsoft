Source and Detector
*******************

**Source of Rays**

The source of rays consists of an annular or rectangular aperture on a planar
surface. The origin of each ray is a random point within
the aperture. The direction of the rays is specified either by
a source at infinity, a source at a finite distance or diffuse.
For a source at infinity all rays are parallel with the direction
set by direction cosines. A source at a finite distance is specified
by a position vector somewhere behind the aperture. Diffuse rays
are generated so as to give a uniform random distribution over a hemisphere.
The total number of rays generated is either set explicitly or by
using an aperture area per ray.

Only one source can be specified. If the source command is used in a
loop then the source will change on each pass through the loop.

The deformation index is used to specify a pixel array which spreads
out a point source into an angular distribution. The deformation data are
set using two functions xsrt.deform() and xsrt.defmat().
If the source is at infinity the x and y sample arrays must be in radians
measured from the direction **sd** along the reference axis **ar** and the
other axis.  (**an** cross **ar**).
If the source is at a finite distance then x and y are displacements
in mm (or whatever distance unit is used) of the position **sp** along
the reference axis **ar** and the other axis (**an** cross **ar**).

**Detector**

The detector consists of an annular or rectangular aperture on a planar
or spherical surface. More than one detector can be specified for an instrument.
Each detector defined will
occupy a given position within the sequence of optical elements specified.
