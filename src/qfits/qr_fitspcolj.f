*+QR_FITSPCOLJ        put integer as table column
        SUBROUTINE QR_FITSPCOLJ(NCOL,NEL,COL,NN,CNAME,NU,CUNIT)
        IMPLICIT NONE
        INTEGER NCOL,NEL,NN,NU
        INTEGER COL(NEL)
        CHARACTER*(NN) CNAME
        CHARACTER*(NU) CUNIT
*NCOL    input        column number
Cf2py  intent(in) ncol
*NEL     input        number of column elements
Cf2py  intent(in) nel
*COL     input        column array
Cf2py  intent(in) col
*NN      input        number of characters in column name
Cf2py  intent(in) nn
*CNAME   input        column name
Cf2py  intent(in) cname
*NU      input        number of characters in column unit
Cf2py  intent(in) nu
*CUNIT   input        column unit
Cf2py  intent(in) cunit
*-Author: Dick Willingale 2014-Dec-26
        INCLUDE 'QR_COM'
        INTEGER IHTYPE,NR,NC,IW,IR,IC
        CHARACTER*(10) KWORD
C
        IF(ISTAT.NE.0) RETURN
C
        IF(IFITS.EQ.0) THEN
                WRITE(*,*) 'QR_FITSPCOLJ - fits file not open'
                ISTAT=1
                RETURN
        ENDIF
C Get type of current extension
        CALL FTGHDT(IFITS,IHTYPE,ISTAT)
C Check extension type is binary table
        IF(IHTYPE.NE.2) THEN
                WRITE(*,*)
     +         'QR_FITSPCOLJ - current extension not binary table'
                ISTAT=1
                RETURN
        ENDIF
C Check dimensions
        CALL FTGNRW(IFITS,NR,ISTAT)
        CALL FTGTCL(IFITS,NCOL,IC,IR,IW,ISTAT)
        IF(NEL.GT.NR*IR) THEN
           WRITE(*,*) 'QR_FITSPCOLD - nrow',NEL,'too large >',NR,'*',IR
           ISTAT=1
           RETURN
        ENDIF
        CALL FTGNCL(IFITS,NC,ISTAT)
        IF(NCOL.GT.NC) THEN
            WRITE(*,*) 'QR_FITSPCOLJ - ncol',NCOL,'too large >',NC
            ISTAT=1
            RETURN
        ENDIF
C Put in keywords for column
        IF(NCOL.LT.10) THEN
                IW=6
                WRITE(KWORD(6:IW),'(I1)') NCOL
        ELSEIF(NCOL.LT.100) THEN
                IW=7
                WRITE(KWORD(6:IW),'(I2)') NCOL
        ELSEIF(NCOL.LT.1000) THEN
                IW=8
                WRITE(KWORD(6:IW),'(I3)') NCOL
        ELSEIF(NCOL.LT.10000) THEN
                IW=9
                WRITE(KWORD(6:IW),'(I4)') NCOL
        ELSE
                IW=10
                WRITE(KWORD(6:IW),'(I5)') NCOL
        ENDIF
        KWORD(1:5)='TTYPE'
        CALL FTUKYS(IFITS,KWORD(1:IW),CNAME,'',ISTAT)
        IF(NU.GT.0) THEN
            KWORD(1:5)='TUNIT'
            CALL FTUKYS(IFITS,KWORD(1:IW),CUNIT,'',ISTAT)
        ENDIF
C Put array into column
        CALL FTPCLJ(IFITS,NCOL,1,1,NEL,COL,ISTAT)
        END
