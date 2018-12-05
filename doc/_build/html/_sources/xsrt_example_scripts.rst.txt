Example Scripts
***************

The following Python and R scripts illustrate how the ray tracing is done.

The sequence of elements is:

        source --> mirrors/stops/lenses/gratings/etc. --> detector

.. code-block:: python

    #!/usr/bin/env python
    # Use Swift XRT geometry as an example
    from __future__ import print_function
    import sys
    import numpy as np
    import images
    import xsrt
    # Useful vectors
    sn=np.array([1,0,0])
    nn=np.array([0,0,0])
    rx=np.array([0,1,0])
    # Set look-up table reflectivity to 1.0
    angs=np.array([0,90])
    refs=np.array([1,1])
    # Support spiders
    spi=np.array([3838.8,0,0])
    tp=np.array([3800,0,0])
    conea=10.05
    nsec=12
    cwid= 0.0
    awid= 3.0
    edf2=np.array([3200,0,0])
    edf1=np.array([3161.2,0,0])
    # Wolter I shell parameters
    fl= 3500
    ph= 3800
    hl= 3200
    ra= 1.0
    ns= 13
    rj=np.array([146.880,140.980,135.320,129.890,124.670,119.660,114.850,
    110.240,105.810,101.560,97.490,93.560,90.833])
    tt=np.array([1.25,1.20,1.15,1.10,1.05,1.00,0.95,0.90,0.85,0.80,0.75,0.70,0.70])
    # Source
    di=np.array([-1,0,0])
    rlim=np.array([92.0,151.0,0,0,0,0])
    nray= 10000
    # Detector
    rdet= 30
    dlim=np.array([0,rdet,0,0])
    dpos=np.array([0,0,0])
    # image paramters
    nx= 100
    ny= 100
    hwid= 5.0
    # Ray tracing calls
    xsrt.reset()
    xsrt.source(1,di,nn,spi,sn,rx,rlim,0.0,nray,0)
    xsrt.surface(1,2,0.0,0.0,0.0,0.0,0.0,0.0,angs,refs,0,0,0)
    xsrt.spider(-conea,spi,sn,rx,nsec,cwid,awid)
    xsrt.spider(0.0,tp,sn,rx,nsec,cwid,awid)
    xsrt.w1nest(fl,rj,ra,fl,ph,hl,fl,tt,tt,tt,sn,rx,nn,0,1,0)
    xsrt.spider(0.0,edf2,sn,rx,nsec,cwid,awid)
    xsrt.spider(conea,edf1,sn,rx,nsec,cwid,awid)
    xsrt.detector(1,dpos,sn,rx,dlim,0.0)
    results=xsrt.trace(0,rdet,-2)
    # Create an image of the detected area
    XD,YD,ZD,XC,YC,ZC,XR,YR,ZR,YDET,ZDET,AREA,IREF=np.loadtxt("detected.dat",skiprows=1,unpack=True)
    arr=images.binxy(YDET,ZDET,0,AREA,-hwid,hwid,-hwid,hwid,nx,ny)
    # Analyse beam to get total collecting area
    images.setfield(nx,-hwid,hwid,ny,-hwid,hwid)
    images.setpos(2,[0,0])
    bb=images.beam(a,hwid,0,0)
    area=bb.flux/100.
    print("area cm^2",area)

.. code-block:: R

    #!/usr/bin/env Rscript
    # Use Swift XRT geometry as an example
    # Useful vectors
    sn<-c(1,0,0)
    nn<-c(0,0,0)
    rx<-c(0,1,0)
    # Set look-up table reflectivity to 1.0
    angs<- c(0,90)
    refs<- c(1,1)
    # Support spiders
    spi<- c(3838.8,0,0)
    tp<- c(3800,0,0)
    conea<- 10.05
    nsec<- 12
    cwid<- 0.0
    awid<- 3.0
    edf2<- c(3200,0,0)
    edf1<- c(3161.2,0,0)
    # Wolter I shell parameters
    fl<- 3500
    ph<- 3800
    hl<- 3200
    ra<- 1.0
    ns<- 13
    rj<-   c(146.880,140.980,135.320,129.890,124.670,119.660,114.850,110.240,
    105.810,101.560,97.490,93.560,90.833)
    tt<- c(1.25,1.20,1.15,1.10,1.05,1.00,0.95,0.90,0.85,0.80,0.75,0.70,0.70)
    # Source
    di<- c(-1,0,0)
    rlim<- c(92.0,151.0)
    nray<- 10000
    # Detector
    rdet<- 30
    dlim=c(0,rdet,0,0)
    dpos<- c(0,0,0)
    # image paramters
    nx<- 100
    ny<- 100
    hwid<- 5.0
    # Ray tracing calls
    qrt_reset()
    qrt_source(1,di,nn,spi,sn,rx,rlim,0.0,nray,0)
    qrt_surface(1,2,0.0,0.0,0.0,0.0,0.0,0.0,angs,refs,0,0,0)
    qrt_spider(-conea,spi,sn,rx,nsec,cwid,awid)
    qrt_spider(0.0,tp,sn,rx,nsec,cwid,awid)
    qrt_w1nest(fl,rj,ra,fl,ph,hl,fl,tt,tt,tt,sn,rx,nn,0,1,0)
    qrt_spider(0.0,edf2,sn,rx,nsec,cwid,awid)
    qrt_spider(conea,edf1,sn,rx,nsec,cwid,awid)
    qrt_detector(1,dpos,sn,rx,dlim,0.0)
    results<- qrt_trace(0,rdet,-2)
    # Create an image of the detected area
    detpos<-read.table("detected.dat",header=TRUE)
    a<-qri_binxy(detpos$YDET,detpos$ZDET,0,detpos$AREA,-hwid,hwid,nx,-hwid,hwid,ny)
    # Analyse beam to get total collecting area
    bb<-qri_beam(a$data_array,hwid,0,0)
    area<-bb$flux/100.
    cat("area cm^2",area,"\n")
