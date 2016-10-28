#!/usr/bin/env python
# Plot the Athena large adaptor SPO module layout
from __future__ import print_function
import numpy as np
import matplotlib.pylab as plt
# Set up SPO modules
from athena_ladapt_modules import *
# Function to return the aperture of each module
def rectangle(x,y,w,h,t):
    cth=np.cos(t)
    sth=np.sin(t)
    xa=w*cth
    xb=w*sth
    ya=h*cth
    yb=h*sth
    x=np.array([x+xa+yb,x+xa-yb,x-xa-yb,x-xa+yb,x+xa+yb])
    y=np.array([y+xb-ya,y+xb+ya,y-xb+ya,y-xb-ya,y+xb-ya])
    return x,y
# Loop for modules
plt.figure()
plt.axes().set_aspect('equal','datalim')
for i in range(nm):
    x,y=rectangle(xm[i],ym[i],wm[i]/2,hm[i]/2,tm[i])
    plt.plot(x,y,"k-")
plt.show()
