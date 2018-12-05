xsrt.functions
**************

**Utility functions**

* rotate() Rotate surface element
* shift() Shift position of surface element
* rseed() Set random number seed
* srtlist() List all current xsrt parameters
* setreset() Reset Fortran common blocks to inital condition

**Apertures, baffles and support structure**

* aperture() Set up an aperture stop
* baffle() Set up a cylindrical baffle
* spider() Set up a support spider

**Wolter systems**

* w1nest() Set up a Wolter Type I nest
* c1nest() Set up conical approximation to a Wolter type I nest
* sipore() Set up Silicon Pore Optics
* spoarr() Set up Silicon Pore Optics array
* wolter2() Set up Wolter Type II surfaces

**Square pore and Kirkpatrick-Baez systems**

* sqpore() Set up slumped square pore MPOs
* sqmpoarr() Set up an array of square pore MPOs
* kbs() Set up a Silicon Kirkpatrick-Baez stack array
* sle() Set up a Schmidt lobster eye stack

**Common optical elements**

* lens() Set up a lens
* prism() Set up a prism
* mirror() Set up a plane mirror
* elips() Set up elliptical grazing indidence mirror
* moa() Set up a cylindrical Micro Optic Array
* opgrat() Set up a single off-plane grating

**Surface quality and deformations**

* surface() Set surface quality parameters
* fresnel() Calculate reflectivity using Fresnel’s equations
* deform() Set up surface deformation dimensions
* defmat() Load deformation matrix
* defparxy() Set up deformation defined by parameters in x and y axes

**Source of rays, tracing rays, detecting rays**

* source() Set up source of rays
* trace() Perform ray tracing
* detector() Set up detector

**Tracing charged particles through magnetic fields**

* bfield() Calculation of magnetic field for array of dipoles
* eltmxt() Trace electrons through SVOM MXT telescope with magnetic diverter
* prtathena() Trace proton through Athena telescope with magnetic diverter

.. automodule:: xsrt
   :members:
