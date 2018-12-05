New Surface Elements
********************
In order to introduce a  new type of ray tracing optical element you
should check the list of available surface types. For example a
spherical mirror at normal incidence using radial limits
would be implemented using TYPE=15. If the required surface type
exists then you only need to
write a QRT\_new routine which provides the QSOFT xsrt interface.
This is easily done
by copying and editing an existing routine, e.g. QRT\_SPIDER.

If you want to add a new element which contains multiple surface elements then
you should start with a routine like QRT_SPOARR. In this case
you will also require
new surface element routines e.g. SRT_SU30 and SRT_SPOARR are
required by QRT_SPOARR.

The steps required to compose, compile and link are as follows:

1. Go to the source directory $QSOFT/src/xsrt

2. Write the new QRT\_new routine as file qrt\_new.f.

        Note down the parameters names required by this routine.

3. Edit the Makefile to include qrt\_new.f in the source file list.

4. Use make to compile the new routine and link the shareble library

        $ make

5. Edit the xsrt.R and xsrt.py script files to include the new function.
   
        If you are going to use IDL move into the qIDL directory and create
        a qrt_new.pro file to define the function for IDL.

        You can use the definitions already present in these files
        to see how the interface works.

If the new optical element is not supported by any existing surface
type then a new type must be invented. The programmer must write
a new SRT\_SUnn routine and modify and existing
or produce a new SRT\_type routine.
A call to the new SRT\_SUnn must also be included in the inner loop
of the srt\_trc.f file. The new SRT\_SUnn and SRT\_type routines
must be edited into the makefile. Otherwise the process is
the same as indicated above.

It is important that the parameters gathered by QRT\_new are packed into
common in the right order so that the relevant surface routine
(SRT\_PLNE etc.) access the parameters correctly. The programmer
must check this by reading the comment lines at the start of
the relevant surface routine.

The routine SRT\_SETF is used to push the parameters into common.
This has the following interface:

.. code-block:: fortran

    *+SRT_SETF      Set surface form and limits parameters
        SUBROUTINE SRT_SETF(NS,IT,NP,P,IDEF,IQ,IH,IM,ISTAT)
        IMPLICIT NONE
        INTEGER NS,IT,NP,IDEF(2),IQ,IH,IM,ISTAT
        DOUBLE PRECISION P(NP)
    *NS     input   surface number (0 for new entry)
    *IT     input   surface type
    *NP     input   number of parameters
    *P      input   array of parameters
    *IDEF   input   deformation
    *IQ     input   surface quality
    *IH     input   hit index (-ve for next in sequence)
    *IM     input   miss index (-ve for next in sequence)
    *ISTAT  in/out  returned status
    *-Author Dick Willingale 1996-Dec-6

NS=0 if you want the surface to be allocated the next free
index in the sequence. IT is the surface index and determines
which SRT\_SUnn routine is going to be called in the ray tracing loop.
Note that the parameters are held in a double precision array.
IDEF and IQ are deformation and surface quality indices that
have already be set by DEFORM and SURFACE commands. If IDEF=0
no deformation will be used. If IQ=0 then the surface will act as a stop.
IH and IM are used
to steer the sequence in the ray tracing. They specify which
surface in the sequence should be next depending on whether or
not the present surface is hit or missed. In most cases IH=-1 and
IM=-1. Examples of cases where a more complicated behaviour is
required are SRT\_PORE and QRT\_SQPORE.

6. Use make to install the new libraries and scripts.
   
        $ make install
