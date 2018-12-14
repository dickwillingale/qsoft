Surface Quality, Reflectivity and Scattering
********************************************

Several surface qualities can be set up for the simulation of a given
instrument. Each is referenced using a surface quality index
(integer 1,2,3...). The type of surface can be reflecting (with reflectivity
specified using Fresnel's equations or a lookup table), refracting
or diffracting. The roughness of the surface can also be specified
using a power law distribution.

The X-ray optical constants **alpha** and **gamma** can be calculated for
a material of specified composition using the function xscat.xopt().
Within the ray tracing the reflectivity is calculated using these
constants using the same code as in function xscat.xfresnel().

The reflectivity as a function of incidence angle in other energy bands can be
calculated from the real and imaginary part of the refractive index using
the function fresnel().

Stops which are intended to block radiation have a surface quality index
set to 0. When rays hit such surfaces they are terminated (absorbed).
Detectors have surface quality index -1. If a ray hits such a surface
it is terminated (detected).
The source aperture surface has quality index -2.
The quality indices of the source, stops and detectors are set automatically.
As ray tracing proceeds rays are stored for further analysis. Each position
along a ray where an intersection with a surface element occured is
labelled with the quality index of the surface.

For a grating the surface type is it=4.
In this case the ruling direction is specified by
the surface element axis and dhub controls the geometry. dhub < 1 in-plane
in which the dhub specifies the d-spacing gradient across the ruling
and dhub > 1 off-plane where the d-spacing gradient along the ruling is
determined from the distance to the hub.
