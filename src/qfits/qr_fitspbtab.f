*+QR_FITSPBTAB        put binary table extension
        SUBROUTINE QR_FITSPBTAB(NROWS,NCOLS,ITY,IRE,IWI)
        IMPLICIT NONE
        INTEGER NROWS,NCOLS
        INTEGER ITY(NCOLS),IRE(NCOLS),IWI(NCOLS)
*NROWS        input        number of rows in table
Cf2py  intent(in) nrows
*NCOLS        input        number of columns in table
Cf2py  intent(in) ncols
*ITY        input        column data type
*                1 integer 32 bit (4 bytes)
*                2 double 64 bit (8 bytes)
*                3 logical 32 bit (4 bytes)
*                4 character 8 bit (1 bytes) per character
*                5 complex 2*32 bit (8 bytes)
*                6 complex 2*64 bit (16 bytes)
*                6 raw 8 bit (1 byte)
*		 7 bit array (1 byte per 8 bit)
Cf2py  intent(in) ity
*IRE        input        column repeat count
Cf2py  intent(in) ire
*IWI        input        width of elements in bytes
Cf2py  intent(in) iwi
*-Author: Dick Willingale 2014-Dec-26
        INCLUDE 'QR_COM'
        INTEGER IHTYPE,NMAX,I,II,IC
        PARAMETER (NMAX=10000)
        CHARACTER*6 DUM1(NMAX)
        CHARACTER*11 DUM2(NMAX)
        CHARACTER*1 DUM3(NMAX)
        CHARACTER*1 CTYP
C
        IF(ISTAT.NE.0) RETURN
C
        IF(IFITS.EQ.0) THEN
                WRITE(*,*) 'QR_FITSPBTAB - fits file not open'
                ISTAT=1
                RETURN
        ENDIF
C Check maximum dimension
        IF(NCOLS.GT.NMAX) THEN
            WRITE(*,*) 'QR_FITSPBTAB - NCOLS',NCOLS,'too large >',NMAX
            ISTAT=1
            RETURN
        ENDIF
C Create new binary table extension
        DO I=1,NCOLS
                WRITE(DUM1(I),'(1A,I5.5)') 'C',I
                II=IRE(I)
                IF(ITY(I).EQ.1) THEN
                        CTYP='J'
                ELSEIF(ITY(I).EQ.2) THEN
                        CTYP='D'
                ELSEIF(ITY(I).EQ.3) THEN
                        CTYP='L'
                ELSEIF(ITY(I).EQ.4) THEN
                        CTYP='A'
                        II=II*IWI(I)
                ELSEIF(ITY(I).EQ.5) THEN
                        CTYP='C'
                ELSEIF(ITY(I).EQ.6) THEN
                        CTYP='M'
                ELSEIF(ITY(I).EQ.7) THEN
                        CTYP='B'
                ELSEIF(ITY(I).EQ.8) THEN
                        CTYP='X'
                ENDIF
                IF(II.LT.10) THEN
                        WRITE(DUM2(I),'(I1,1A)') II,CTYP
                        IC=2
                ELSEIF(II.LT.100) THEN
                        WRITE(DUM2(I),'(I2,1A)') II,CTYP
                        IC=3
                ELSEIF(II.LT.1000) THEN
                        WRITE(DUM2(I),'(I3,1A)') II,CTYP
                        IC=4
                ELSEIF(II.LT.10000) THEN
                        WRITE(DUM2(I),'(I4,1A)') II,CTYP
                        IC=5
                ELSE
                        WRITE(DUM2(I),'(I5,1A)') II,CTYP
                        IC=6
                ENDIF        
                IF(ITY(I).EQ.4) THEN
                        IF(IWI(I).LT.10) THEN
                                WRITE(DUM2(I)(IC+1:),'(I1)') IWI(I)
                        ELSEIF(IWI(I).LT.100) THEN
                                WRITE(DUM2(I)(IC+1:),'(I2)') IWI(I)
                        ELSEIF(IWI(I).LT.1000) THEN
                                WRITE(DUM2(I)(IC+1:),'(I3)') IWI(I)
                        ELSEIF(IWI(I).LT.10000) THEN
                                WRITE(DUM2(I)(IC+1:),'(I4)') IWI(I)
                        ELSE
                                WRITE(DUM2(I)(IC+1:),'(I5)') IWI(I)
                        ENDIF
                ENDIF
                DUM3(I)=' '
        ENDDO
        CALL FTIBIN(IFITS,NROWS,NCOLS,DUM1,DUM2,DUM3,'',
     +        0,ISTAT)
        IEXT=IEXT+1
C Move to header
        CALL FTMAHD(IFITS,IEXT,IHTYPE,ISTAT)
        END
