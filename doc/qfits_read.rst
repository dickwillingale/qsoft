Reading FITS Files
******************

The entire contents of a FITS file can be read into a Python or R object
using a single call as illustrated below.

.. code-block:: python

    #!/usr/bin/env python
    import qfits
    filename="test_fitsnew.fits"
    fitsin=qfits.fitsfile(filename)
    # Print summary
    fitsin.display()
    # Access the number of HDUs in FITS file
    print("number of hdu",fitsin.nhdu)
    # The primary array
    print(fitsin.hdu[0].data_array)
    # All the keywords on extension 1
    print("HDU 1 keywords",fitsin.hdu[1].kw)
    # A Particular keyword on extension 1
    print("HDU 1 keyword NAXIS1",fitsin.hdu[1].kw["NAXIS1"])
    # The complete table on extension 2
    print("HDU 2 table",fitsin.hdu[2].table)
    # A particular column from table on extension 2
    print("HDU 2 table column s",fitsin.hdu[2].table["s"])

.. code-block:: R

    #!/usr/bin/env Rscript
    filename<-"test_fitsnew.fits"
    fitsin<- qr_fitsread(filename)
    # Print a summary
    qr_fitsprint(fitsin)
    # Get the number of HDU (primary + extensions)
    cat("number of HDU",fitsin$NHDU,"\n")
    # Print the primary data arra
    print(fitsin$primary$DATA_ARRAY)
    # Access a particular keyword in extension 2
    cat("HDU extension 2 keyword TESTD",fitsin$extension[[2]]$TESTD,"\n")
    # All the keywords on extension 2
    print(fitsin$extension[[2]])
    # The complete table on extension 2
    print(fitsin$extension[[2]]$table)
    # A particular column from table on extension 2
    print(fitsin$extension[[2]]$table[,3])

Because the indexing of lists and arrays starts at zero in Python and 1 in R
the internal structure of the object returned is different. 

In Python the primary HDU is hdu[0] and extensions are hdu[1], hdu[2] etc..

In R the primary is called primary and the extensions are a list
accessed as extension[[1]], extension[[2]] etc..

In Python the keywords are stored in a dictionary and
particular keywords are accessed using a name index.

In R the
keywords are stored as named variables within the extension list.

The summary listings produced by fitsfile.display() in Python and
qr_fitsprint() in R can be used to reveal the way in which the contents of
the FITS file are stored in memory.

Details on how to access all the elements of FITS files can be found in the
source for the Python class fitsfile and R function qr_fitsread() which are
defined in $QSOFT/src/qfits/qfits.py and $QSOFT/src/qfits/qfits.R.
