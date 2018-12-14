Monte Carlo and Random Numbers
******************************

The starting positions of rays, X-ray scattering angles from surface roughness and
some surface figure errors/deformations are generated using random numbers.
The sequence of random numbers used will be different each run unless the
random number seed is set using function rseed(). If the same seed value is set
before calling the function trace() exactly the same random sequence will be
generated for the ray tracing and the results will be identical.
