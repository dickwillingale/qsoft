#!/usr/bin/env python
# Plot the reflectivity curves for Ir with B4C overcoat
from __future__ import print_function
import numpy as np
import xscat
import matplotlib.pylab as plt
from ir_b4cover import *
# Set up array of energies
emin= 0.1
emax= 15.0
esam= 0.05
ekev= np.arange(emin,emax,esam)
# Get reflectivities
iangle,refs=ir_b4c(ekev)
#
ne= ekev.size
na=iangle.size
# Plot
plt.figure()
for i in range(na):
    plt.plot(ekev,refs[i,:],"k-")
plt.xlabel("keV")
plt.ylabel("reflectivity")
plt.title("Iridium + B4C overcoat")
plt.show()
