## @package images
# Functions for image display and manipulation
import numpy as np
import matplotlib.pylab as plt
# Python QSOFT locator
def plt_show_locator(fg,npos):
    """Get local coordinate positions using the cursor

    Args:
        fg:     figure id returned by plt.figure()
        npos:   number of positions to be retured (clicked)
    
    Returns:
        xx,yy 2 arrays containing npos local coordinates
    
    This routine provides a basic level of interaction with a displayed figure.
    It is used in place of a simple plt.show() call so that npos local coordinate
    positions can be selected interactively using the cursor from a displayed image
    and returned in arrays to the user.
    """
    xx=np.empty(npos)
    yy=np.empty(npos)
    ipos=0
    def onclick(event):
        nonlocal xx
        nonlocal yy
        nonlocal ipos
        nonlocal npos
        xx[ipos]=event.xdata
        yy[ipos]=event.ydata
        ipos=ipos+1
        if ipos==npos:
            plt.close()
    fg.canvas.mpl_connect("button_press_event",onclick)
    plt.show()
    return xx,yy
