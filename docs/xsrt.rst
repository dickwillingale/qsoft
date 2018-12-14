xsrt   - Sequential Ray Tracing
*******************************

The code was written for modelling grazing incidence X-ray telescopes but it
works at normal incidence and includes many common optical elements.

The code is *sequential* in the sense
that rays encounter the optical elements in the order that they are specified.
However, by setting flags associated
with each element the sequence can be controlled dynamically to handle
multiple interactions, between different optical elements, in any sequence.

.. toctree::

   xsrt_example_scripts
   xsrt_elements
   xsrt_source_detector
   xsrt_random
   xsrt_deformations
   xsrt_surfaces
   xsrt_telescopes
   xsrt_tracing
   xsrt_functions
