Deformations
************

Deformations of surfaces can be specified using matricies which either span
a grid of points in the local coordinate system of the surface or
are indexed using integer labels for sectors or areas on a surface.
A set of deformations pertaining to a single surface or group of
related surfaces are given a deformation index (integer 1,2,3...).
The positions of the deformation grid points in local coordinates are
specified by two 1-dimensional arrays.

Alternatively deformations can be specified using functions with parameters set
separately for local x and y coordinates.

A surface deformation is applied along the normal to the surface. The
deformation value is interpolated from the 2-d grid of points.

Radial deformations for annula apertures are specified by a vector, sampling
in 1-d in azimuth, and the deformation is applied as a perturbation in the radial
direction.

The function xsrt.deform() is used to set up the type and dimensionality of a
particular deformation and must be the first call. The component matricies or
function parameters are then set using calls to xsrt.defmat() or
xsrt.defparxy().

A deformation applied to the source() spreads the point source into a pixel
array. See **Source of Rays**.
