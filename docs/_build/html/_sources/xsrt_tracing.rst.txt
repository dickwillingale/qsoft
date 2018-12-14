Ray Tracing and Saving Rays
***************************

Once the source, detector and other elements have been defined rays can be
traced through the instrument using the function trace(). The form of the
output is controlled by the parameter **iopt**.

* -2 save traced.dat and detected.dat files
* -1 save detected.dat
* 0 don't save files or adjust focus
* 1 adjust focus and save detected.dat
* 2 adjust focus and save detected.dat and traced.dat
* Only rays with **iopt** reflections are used in adjustment

The files traced.dat and detected.dat are ASCII tabulations. 

When **iopt** is +ve then the detector position which gives the best focus is
determined. Only rays which have **iopt** reflections and impact the detector
within a radius **riris** of the centre of the detector are included in the
analysis. The detector is shifted along the normal direction to find
the axial position of minimum rms radial spread. The results of this
analysis are returned as:

* **area**    detected area within RIRIS
* **dshft**   axial shift to optimum focus (0.0 if IOPT<=0)
* **ybar**    y centroid of detected distribution
* **zbar**    z centroid of detected distribution
* **rms**     rms radius of detected distribution

The file traced.dat contains the paths of all the rays. It can be very large
so should not be saved unless required for detailed analysis.

* **RXP,RYP,RZP**  positions of points along each ray
* **AREA** aperture area associated with ray
* **IQU**  quality index -2 at source, 1 reflected, 0 absorbed, -1 detected

Note: in the tabulation the beginning of each ray is identified using **IQU=-2**
and the end using **IQU=0** absorbed or **IQU=-1** detected. Using these data
you can plot the paths of all the rays.

The file detected.dat contains information about the detected rays.

* **XD,YD,ZD**  the detected position for each ray
* **XC,YC,ZC**  the direction cosines for each ray
* **XR,YR,ZR**  the position of the last interaction before detection
* **YDET,ZDET** the local detected position on detector
* **AREA**  the aperture area associated with the ray
* **IREF** the number of reflections suffered by the ray

The position **XR,YR,ZR** is used to indicate where the ray came from.

The following snippets of code show how an image of the detected rays can be
generated in Python or R.

.. code-block:: python

    import numpy as np
    import images
    import xsrt
    ...
    ...
    # half width of image mm
    hwid=5.0
    # trace all the rays
    results=xsrt.trace(0,rdet,-2)
    # Create an image of the detected area
    XD,YD,ZD,XC,YC,ZC,XR,YR,ZR,YDET,ZDET,AREA,IREF=np.loadtxt("detected.dat",
        skiprows=1,unpack=True)
    arr=images.binxy(YDET,ZDET,0,AREA,-hwid,hwid,-hwid,hwid,nx,ny) 

.. code-block:: R

    # half width of image mm
    hwid<- 5.0
    # trace all the rays
    results<- qrt_trace(0,rdet,-2)
    # Create an image of the detected area
    detpos<-read.table("detected.dat",header=TRUE)
    aim<-qri_binxy(detpos$YDET,detpos$ZDET,0,detpos$AREA,-hwid,hwid,nx,-hwid,hwid,ny)

In Python **arr** is an image array. In R **aim** is an image object which
contains the image array **aim$data_array**.
In both cases the function images.binxy() is used to bin up the aperture area
associated with each ray into an image (2-d histogram). The effective area
is found by summing up areas of the image.
