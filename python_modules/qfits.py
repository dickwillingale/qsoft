# qfits - qpy fitsio
from __future__ import print_function
import qfitsfor 
import numpy as np
from datetime import datetime
# Initialisation
def init():
    qfitsfor.qrf_init()
init()
# Get column from table on FITS file
def fitsgcol(ic,typ,nrows,rp,vr):
# ic        column number
# typ        column data type
# nrows        number of column rows
# rp        repeat count
# vr        If rp zero then vr[nrows] is variable count for each row
# The following parameters are used:
#      nc - the number of calls required to get complete column
#      ne -length returned per call
#      qt 1 complete column returned in 1 call
#      qt 2 1 call per row fixed width returned as a list
#      qt 3 variable width column 1 call per row returned as a list
    if rp==0:
# Variable width column and request for each row
        nc=nrows
        ne=np.array(vr,int)
        vc=[]
        qt=3
    else:
# Fixed width column
        if typ==4 or typ==8:
# character - get 1 string per row because passing character
#             arrays to/from R is not recommended.
# bit - get 1 character string per row - charact ers set '0 or '1'
            nc=nrows
            ne=np.ones(nrows,int)*int(rp)
            vc=[]
            qt=2
        else:
            nc=1
            ne=np.array([nrows*rp],int)
            qt=1
# Loop to get columnn values
    for ii in range(nc):
        if qt==1:
            vc=fitsgcolv(typ,ic,ii+1,ne[ii])
            if rp>1:
                vc.shape=[nrows,int(rp)]
        else:
            vc.append(fitsgcolv(typ,ic,ii+1,ne[ii]))
    return vc
# Get column values from FITS file
def fitsgcolv(typ,ic,ii,ne):
    if ne==0:
        return None
    if typ==1:
        return qfitsfor.qr_fitsgcvj(0,ic,ii,ne)
    elif typ==2:
        return qfitsfor.qr_fitsgcvd(np.nan,ic,ii,ne)
    elif typ==3:
        a=qfitsfor.qr_fitsgcvl(ic,ii,ne)
        return a.astype(np.bool)
    elif typ==4:
        a=qfitsfor.qr_fitsgcs(ic,ii,ne)
        if a[2]==1:
            s=None
        else:
            s=a[0]
        return s
    elif typ==5:
        return qfitsfor.qr_fitsgcvc(np.nan,ic,ii,ne)
    elif typ==6:
        return qfitsfor.qr_fitsgcvm(np.nan,ic,ii,ne)
    elif typ==7:
        a=qfitsfor.qr_fitsgcvj(0,ic,ii,ne)
        return a.astype(np.int8)
    elif typ==8:
        return qfitsfor.qr_fitsgcx(ic,ii,ne)
# Set current fits HDU and return info
class hduinfo: pass
def fitshduinfo(ihdu):
    a=qfitsfor.qr_fitshdu(ihdu,10)
    b=hduinfo()
    b.hdutype=a[0]
    b.naxis=a[1]
    b.naxes=a[2][:a[1]]
    if b.naxis>0:
        b.nrows=b.naxes[0]
        b.nels=np.prod(b.naxes)
    else:
        b.nrows=0
        b.nels=0
    if b.naxis>1:
        b.ncols=b.naxes[1]
    else:
        b.ncols=0
    if b.hdutype==0: b.ncols=1
# numpy uses the C array index order, FITS uses the Fortran order
    b.naxes=b.naxes[::-1]
    b.nkeys=a[3]
    return b
# Get key data
def fitsgetkey(ikey):
    a=qfitsfor.qr_fitsgetkey(ikey)
    ki=a[1]
    key=a[0][:ki].decode()
    si=a[3]
    sval=a[2][:si].decode()
    jval=a[4]
    dval=a[5]
    lval=a[6]
    ktype=a[7]
    ci=a[9]
    cval=a[8][:ci].decode()
    return key,ki,sval,si,jval,dval,lval,ktype,cval,ci
# Get data types
def fitstypes(hdutype,ncols):
    a=qfitsfor.qr_fitstypes(hdutype,ncols)
    ctype=a[0]
    rp=a[1]
    return ctype,rp
# Open fits file for read/write
def fitsupdate(filename):
    return qfitsfor.qr_fitsopen(len(filename),filename,1)
# Get column info
class fitstab: pass
def fitscolnam(ic,rr,nrows):
    a=qfitsfor.qr_fitscolnam(ic,rr,nrows,500)
    b=fitstab()
    b.colnam=a[0].decode()
    b.iname=a[1]
    b.vrep=a[2]
    return b
# Put array onto fits file as primary or extension data array
def fitsparr(arr):
    nels=np.array(arr.shape,int)
    nel=arr.size
# Reverse the array indices for FITS which is Fortran like
    nelr=nels[::-1]
    arr.shape=[nel]
    if arr.dtype==float:
        qfitsfor.qr_fitsparrd(nelr,arr)
    else:
        qfitsfor.qr_fitsparrj(nelr,arr)
    arr.shape=nels
# Define fitshdu object
class fitshdu:
    def __init__(self,ihdu,hdutype):
        self.ihdu=ihdu
        self.htype=hdutype
# Comment records are held in a list
        self.cr=[]
# Keywords are held in a dictionary
        self.kw={}
# Keyword comments are held in a dictionary
        self.kc={}
# data_array or table
        if(hdutype==0):
            self.data_array=None
        else:
            self.table={}
            self.units={}
# List fits HDU object
    def display(self):
        print("************* HDU",self.ihdu,"***************")
        print("hdu.htype      ",self.htype)
        for ic in range(len(self.cr)):
            print(" ",self.cr[ic])
        for key in self.kw.keys():
            if self.kc[key]!=None:
                print(key," = ",self.kw[key]," : ",self.kc[key])
            else:
                print(key," = ",self.kw[key])
        if self.htype==0:
            if self.data_array is None:
                print("hdu.data_array")
            else:
                print("hdu.data_array",self.data_array.dtype)
            print(self.data_array)
        else:
            print("hdu.table")
            for key in self.table:
                if type(self.table[key])==list:
                    print(key,"list")
                else:
                    print(key,self.table[key].dtype)
                print("  ",self.table[key])
# Define fitsfile object
class fitsfile:
    def display(self):
        print("filename       ",self.filename)
        print("date           ",self.date)
        print("nhdu           ",self.nhdu)
        for i in range(self.nhdu):
            self.hdu[i].display()
    def __init__(self,filename):
        self.filename=filename
        self.date=str(datetime.now())
        self.nhdu=0
        self.hdu=[]
# Creating a new file return
        if filename=="new": return
# Read contents of existing fits file into object
# Open file and get number of header units
        nhdu=qfitsfor.qr_fitsopen(len(filename),filename,0)
        self.nhdu=nhdu
# Loop for all header units
        for i in range(nhdu):
# Get info for this header unit
            a=fitshduinfo(i+1)
            hdutype=a.hdutype
            naxis=a.naxis
            nrows=a.nrows
            ncols=a.ncols
            naxes=a.naxes
            nkeys=a.nkeys
            nels=a.nels
# Append instance of fitshdu() to current list
            hdu=fitshdu(i,hdutype)
            hdu.htype=hdutype
            self.hdu.append(hdu)
# Loop for all keywords/records (including comments etc)
            for k in range(nkeys):
                key,ki,sval,si,jval,dval,lval,ktype,cval,ci=fitsgetkey(k+1)
                if ktype==1:
                    val=jval
                elif ktype==2:
                    val=dval
                elif ktype==3:
                    val=bool(lval)
                elif ktype==4:
                    if si>0:
                        val=sval[:si]
                    else:
                        val=" "
                if key==""or key=="COMMENT"or key=="HISTORY"or key=="CONTINUE":
                    if ci>0:
                        hdu.cr.append(cval[:ci])
                else:
                    hdu.kw[key]=val
                    if ci>0:
                        hdu.kc[key]=cval[:ci]
                    else:
                        hdu.kc[key]=None
# Get data types
            ctype,rp=fitstypes(hdutype,ncols)
            if hdutype==0 and naxis>0:
# Primary data array
                if ctype==1:
                    dd=qfitsfor.qr_fitsgpvj(0,nels)
                    dd.shape=naxes[:naxis]
                    hdu.data_array=dd
                else:
                    dd=qfitsfor.qr_fitsgpvd(np.nan,nels)
                    dd.shape=naxes[:naxis]
                    hdu.data_array=dd
            if hdutype>0 and naxis>0:
# data tables saved as dictionaries
                for ic in range(ncols):
                    a=fitscolnam(ic+1,rp[ic],nrows)
                    colnam=a.colnam[:a.iname]
                    if ctype[ic]==8: colnam=colnam+"_bits"
                    vr=a.vrep
# Look for TUNIT keyword for this column
                    tkw="TUNIT"+str(ic+1)
                    if tkw in hdu.kw:
                        hdu.units[colnam]=hdu.kw[tkw]
                    else:
                        hdu.units[colnam]=""
# Get column data
                    hdu.table[colnam]=fitsgcol(ic+1,ctype[ic],nrows,rp[ic],vr)
# Close fits file
        qfitsfor.qr_fitsclose()
# Save to  file
    def save(self,filename):
# Open new fits file
        qfitsfor.qr_fitsnew(len(filename),filename)
        self.filename=filename
        self.date=str(datetime.now())
# Loop for all HDU
        for i in range(self.nhdu):
            hdu=self.hdu[i]
            if hdu.htype==0:
                if hdu.data_array is None:
                    qfitsfor.qr_fitsempty()
                else:
                    fitsparr(hdu.data_array)
            else:
                if len(hdu.table)==0:
                    qfitsfor.qr_fitsempty()
                else:
                    fitsptab(hdu)
            for ic in range(len(hdu.cr)):
                    qfitsfor.qr_fitspcom(len(hdu.cr[ic]),hdu.cr[ic])
            for key in hdu.kw.keys():
                va=hdu.kw[key]
                co=hdu.kc[key]
                ty=type(va)
                if ty==str:
                    qfitsfor.qr_fitspkeys(key,va,co)
                elif ty==float:
                    qfitsfor.qr_fitspkeyd(key,va,co)
                elif ty==int:
                    qfitsfor.qr_fitspkeyj(key,va,co)
                elif ty==bool:
                    qfitsfor.qr_fitspkeyl(key,va,co)
# Close fits file
        qfitsfor.qr_fitsclose()
# Create binary table HDU on fits file
def fitsptab(hdu):
    tab=hdu.table
    units=hdu.units
#    names=tab.keys()
    names=list(tab)
    ncols=len(names)
    print(names)
    nrows=len(tab[names[0]])
    ity=np.empty(ncols,int)
    ire=np.empty(ncols,int)
    iwi=np.empty(ncols,int)
    iss=np.empty(ncols,int)
    for i in range(ncols):
        nam=names[i]
        if type(tab[nam])==list:
            ty=type(tab[nam][0])
            irr=len(tab[nam][0])
            iss[i]=0
        else:
            ty=tab[nam].dtype
            if tab[nam].ndim>1:
                irr=tab[nam].shape[1]
                iss[i]=1
            else:
                irr=1
                iss[i]=0
        if ty==np.int8:
                ity[i]=7
                ire[i]=irr
                iwi[i]=1
        elif ty==np.int32:
            ity[i]=1
            ire[i]=irr
            iwi[i]=4
        elif ty==np.int64:
            ity[i]=1
            ire[i]=irr
            iwi[i]=8
        elif ty==float:
            ity[i]=2
            ire[i]=irr
            iwi[i]=8
        elif ty==bool:
            ity[i]=3
            ire[i]=irr
            iwi[i]=4
        elif ty==str:
            if names[i].endswith("_bits"):
                ity[i]=8
                ire[i]=irr
                iwi[i]=0
            else:
                ity[i]=4
                ire[i]=1
                iwi[i]=0
            for k in range(nrows):
                ll=len(tab[nam][k])
                if ll>iwi[i]:
                    iwi[i]=ll
        elif ty==np.complex64:
            ity[i]=5
            ire[i]=irr
            iwi[i]=8
        elif ty==np.complex128:
            ity[i]=6
            ire[i]=irr
            iwi[i]=16
# Set binary table header in FITS file
    qfitsfor.qr_fitspbtab(nrows,ity,ire,iwi,ncols)
# Put column data into table
    for i in range(ncols):
        nam=names[i]
        if nam in units:
             cunit=units[nam]
        else:
             cunit=""
        if iss[i]==1:
# Vectorize 2-D column values ready for put call
             ss=tab[nam].shape
             nn=tab[nam].size
             tab[nam].shape=[nn]
        if ity[i]==1:
# Put integer column data on to fits file
            qfitsfor.qr_fitspcolj(i+1,tab[nam],len(nam),nam,len(cunit),cunit)
        elif ity[i]==2:
# Put double precision column data on to fits file
            qfitsfor.qr_fitspcold(i+1,tab[nam],len(nam),nam,len(cunit),cunit)
        elif ity[i]==3:
# Put logical column data on to fits file
            lco=tab[nam].astype(int)
            qfitsfor.qr_fitspcoll(i+1,lco,len(nam),nam,len(cunit),cunit)
        elif ity[i]==4:
# Put string column data (list) on to fits file 
            for k in range(nrows):
                s=tab[nam][k]
                ls=len(s)
                qfitsfor.qr_fitspcols(i+1,k+1,ls,s,len(nam),nam,len(cunit),cunit)
        elif ity[i]==5:
# Put complex column data on to fits file
            qfitsfor.qr_fitspcolc(i+1,tab[nam],len(nam),nam,len(cunit),cunit)
        elif ity[i]==6:
# Put double complex column data on to fits file
            qfitsfor.qr_fitspcolm(i+1,tab[nam],len(nam),nam,len(cunit),cunit)
        elif ity[i]==7:
# Put raw byte column data on to fits file
            b=tab[nam].astype(int)
            qfitsfor.qr_fitspcolb(i+1,b,len(nam),nam,len(cunit),cunit)
        elif ity[i]==8:
# Put bit column data on to fits file
# Each row of bits is a string of "1" or "0" characters
# Strip "_bits" from nam
            ln=len(nam)-5
            fna=nam[:ln]
            for k in range(nrows):
                s=tab[nam][k]
                ls=len(s)
                qfitsfor.qr_fitspcolx(i+1,k+1,ls,s,len(fna),fna,len(cunit),cunit)
        if iss[i]==1:
# Reset shape for 2-D column values
             tab[nam].shape=ss
