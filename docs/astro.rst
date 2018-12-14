astro  - Astronomy Applications
*******************************

X-ray spectra. X-ray spectral data are usually presented as histograms of 
X-ray counts vs. energy (keV). The functions btoc() and ctob() provide conversion between energies at n bin centres and n+1 bin boundaries.

* brems() Bremsstrahlung spectrum
* habs() X-ray absorption by a Hydrogen column density
* setabnd() Set abundances for XSPEC routines used in absorption and optical depth calculations
* btoc() Convert n+1 bin boundaries to n bin centres
* ctob() Convert n bin centres to n+1 bin boundaries

X-ray optical depth. The local ISM and cosmic IGM gas is modelled as a hydrogen column including heavier elements at specified abundances and ionisation state.

* igmtau() X-ray optical depth of IGM gas
* iigmtau() X-ray optical depth of ionized IGM gas
* iigmtauvz() X-ray optical depth of ionized IGM gas vs. redshift
* ismtau() X-ray optical depth of cold ISM gas
* iismtau() X-ray optical depth of ionized ISM gas
* lyftau() X-ray optical depth of Lyman Forest
* lyftauvz() X-ray optical depth of Lyman Forest vs. redshift

Cosmology and redshift.

* cosmo() Calculation of cosmological quantities (luminosity distance etc.)
* kcorrb() K-correction using numerical integration of the Band function

.. toctree::

   astro_functions
