images  - Image Processing
**************************

Positions in images. The coordinates of an image field are set up
using the function setfield(). The *current position* can be set using the
function setpos(). This position is then used by functions like beam().
The position is
given in pixel coordinates (a real scale running 0-NCOLS in X and 0-NROWS in Y)
and local X,Y (usually taken as mm). If, in addition to setfield(),
sky coordinates are set using the function setsky() the *current position* is
also specified in local 
azimuth and elevation (degrees), Celestial RA,DEC (degrees J2000), Ecliptic
EA,EL (degrees) and Galactic LII, BII (degrees). The sky coordinates can be
set using different projections, Plate Carre, Aitoff (Hammer) or Lambert
(equatorial aspect of Azimuthal equal-area).

* setfield() Set local coordinates for image field
* setsky() Set up sky coordinates for image field
* setpos() Set current position in image field
* getpos() Get current position in image field
* toxy() Convert position to local xy coordinates
* plt_show_locator() Get local coordinate positions using the cursor

Analysis of a beam containing a source or PSF. The beam is centred at the
*current positon* (see above).

* beam() Analysis of source above background within a circular beam
* sqbeam() Analysis of source above background within a square beam
* lecbeam() Analysis of source above background in a lobster eye cross-beam

Creation of images from event lists or 2-d functions.
In Python these function return a 2-d array. In R they return
an image object which contains the array and ancillary information.

* binxy() x-y event binning to form an image array
* lebin() Create an image from an event list binning using lobster eye psf
* lecimage() Create an image array of the lobster eye cross-beam
* lepsf() Create an image of the lobster eye PSF

Drawing over images. The function hamgrid() and lamgrid() will only work if
sky coordinates have been set up using setsky().

* rectangles() Draw rectangles
* hamgrid() Plot a Hammer projection grid on figure
* lamgrid() Plot a Lambert projection grid on figure

1-d profiles

* gaussian() Gaussian profile
* king_profile() King function (modified Lorentzian) profile
* lorentzian() Lorentzian profile

Model function fitting using a statistic. Can be used for fitting of PSF profiles to image data or more generally for fitting data with a model function.

* srchmin() Search for minimum statistic and return best fit parameters and confidence limits of the parameters
* peakchisq() Chi-squared for image peak fitting
* quaderr() Quadratic estimator for confidence limit

.. toctree::

   images_functions
