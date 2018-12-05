Writing FITS Files
******************

FITS files can be created using a Python or R script as
illustrated below. Python uses methods in class fitsfile and
class fitshdu. R uses functions qr_fits*().

.. code-block:: python

    #!/usr/bin/env python
    import numpy as np
    import qfits
    #
    a=qfits.fitsfile("new")
    # Primary array
    hdu=qfits.fitshdu(0,0)
    arr=np.arange(100)+101.1
    arr.shape=[20,5]
    hdu.data_array=arr
    a.hdu.append(hdu)
    # Extension Arrays
    hdu=qfits.fitshdu(1,0)
    iarr=np.arange(100)+1
    iarr=iarr[::-1]
    iarr.shape=[20,5]
    hdu.data_array=iarr
    # Add some comments
    hdu.cr.append("This is a test comment")
    hdu.cr.append("This is some history")
    # keywords
    hdu.kw["TESTD"]=55.5
    hdu.kc["TESTD"]="a test double value"
    hdu.kw["TESTJ"]=28
    hdu.kc["TESTJ"]="a test integer value"
    hdu.kw["TESTL"]=True
    hdu.kc["TESTL"]="a true logical value"
    hdu.kw["TESTF"]=False
    hdu.kc["TESTF"]="a false logical value"
    hdu.kw["TESTS"]="a string"
    hdu.kc["TESTS"]="a test string value"
    a.hdu.append(hdu)
    # Tables
    hdu=qfits.fitshdu(2,2)
    ix=np.arange(4)+1
    x=np.array(ix*10+ix*.1,float)
    ix=np.array(ix*10,int)
    y=["aaaaaaaaaaa","bbbbbbbbbbb","ccccccccccc","ddddddddddd"]
    z=np.array([True,True,False,True])
    q=np.array([9,8,7,6],np.int8)
    s=np.array([55+22j,66+23j,77+24j,88+25j])
    xx=np.array(np.arange(12)+10,float)
    xx.shape=[4,3]
    b_bits=["1000000111","0111001101","0010000010","0001101100"]
    rnames=["r1","r2","r3","r4"]
    short=np.array([1.1,2.2,3.3])
    hdu.table["ix"]=ix
    hdu.units["ix"]="splogs"
    hdu.table["x"]=x
    hdu.table["xx"]=xx
    hdu.table["y"]=y
    hdu.table["z"]=z
    hdu.table["q"]=q
    hdu.table["s"]=s
    hdu.table["b_bits"]=b_bits
    hdu.table["short"]=short
    hdu.table["rnames"]=rnames
    a.hdu.append(hdu)
    #
    a.nhdu=3
    # Print a summary and save
    a.display()
    a.save("test_fitsnew.fits")

.. code-block:: R

    #!/usr/bin/env Rscript
    qr_fitsnew("fitswrite_test.fits")
    # Primary array
    arr<-seq(from=10.1,to=20,length.out=100)
    dim(arr)<-c(5,20)
    qr_fitsparrd(arr)
    # Extension Arrays
    iarr<-100:1
    dim(iarr)<-c(5,20)
    qr_fitsparrj(iarr)
    # Comments and history cards
    qr_fitspcom("This is a test comment")
    qr_fitsphis("This is some history")
    qr_fitsempty()
    # keywords
    qr_fitspkeyd("TESTD",55.5,"a test double value")
    qr_fitspkeyj("TESTJ",28,"a test integer value")
    qr_fitspkeyl("TESTL",T,"a true logical value")
    qr_fitspkeyl("TESTF",F,"a false logical value")
    qr_fitspkeys("TESTS","a string","a test string value")
    # Tables
    x<- c(10.10,20.20,30.30,40.40)
    ix<- as.integer(c(10,20,30,40))
    y<- c("aaaaaaaaaaa","bbbbbbbbbbb","ccccccccccc","ddddddddddd")
    z<- c(T,T,F,T)
    q<-as.raw(c(9,8,7,6))
    s<- complex(real=c(55,66,77,88),imaginary=c(22,23,24,25))
    rnames<- c("r1","r2","r3","r4")
    tt<- data.frame(x,ix,y,z,q,s,stringsAsFactors=F,row.names=rnames)
    qr_fitspobj(tt,"dframe")
    qr_fitsclose()
    # Print a summary
    fitsin<- qr_fitsread("fitswrite_test.fits")
    qr_fitsprint(fitsin)
