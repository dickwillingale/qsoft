#!/usr/bin/env python
# Plot area vs keV for SVOM MXT
from __future__ import print_function
import numpy as np
import matplotlib.pylab as plt
# read tabulation file
a=np.loadtxt("svom_mxt_area_vs_kev_v4.dat",skiprows=1)
ekev=a[:,0]
farea=a[:,1]
tarea=a[:,2]
fwhm=a[:,3]
# Set ranges for plot
emin=0.1
emax=10.0
amin=0
amax=80
# 
plt.figure()
plt.axis([emin,emax,amin,amax])
plt.title('SVOM MXT effective area')
plt.xlabel('E keV')
plt.ylabel(r'area cm$^2$')
plt.plot(ekev,tarea)
plt.plot(ekev,farea,color='red')
plt.show()
