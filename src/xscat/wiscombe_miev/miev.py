# Mie Scattering Code for Python
# Shareble object version of Wiscombe MIEV subroutine
from __future__ import print_function
import mievfor
import numpy as np
#
print(mievfor.miev0.__doc__)
#
class mievs:
    def list(self):
        print("qext",self.qext)
	print("qsca",self.qsca)
	print("gqsc",self.gqsc)
	print("pmom",self.pmom)
	print("sforw",self.sforw)
	print("sback",self.sback)
	print("s1",self.s1)
	print("s2",self.s2)
	print("tforw",self.tforw)
	print("tback",self.tback)
	print("spike",self.spike)
def miev0(xx,crefin,perfct,mimcut,anyang,xmu,nmom,ipolzn,momdim,prnt):
    a=mievfor.miev0(xx,crefin,perfct,mimcut,anyang,xmu,nmom,ipolzn,momdim,prnt)
    b=mievs()
    b.qext=a[0]
    b.qsca=a[1]
    b.gqsc=a[2]
    b.pmom=a[3]
    b.sforw=a[4]
    b.sback=a[5]
    b.s1=a[6]
    b.s2=a[7]
    b.tforw=a[8]
    b.tback=a[9]
    b.spike=a[10]
    return b

