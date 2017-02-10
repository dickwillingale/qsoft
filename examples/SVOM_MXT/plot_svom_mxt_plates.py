#!/usr/bin/env python
# Plot layout of SVOM MXT plates
from __future__ import print_function
import numpy as np
import matplotlib.pylab as plt
# Set up plates
from svom_mxt_plates_v4 import *
# Function to return the rectangular aperture of each module
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
# Loop to plot plates
plt.figure()
plt.axes().set_aspect('equal','datalim')
for i in range(nd):
    x,y=rectangle(xm[i],ym[i],wm[i]/2,hm[i]/2,tm[i])
    plt.plot(x,y,"k-")
plt.show()
