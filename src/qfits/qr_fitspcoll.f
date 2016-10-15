*+QR_FITSPCOLL        put logical as table column
        SUBROUTINE QR_FITSPCOLL(NCOL,NEL,COL,NN,CNAME,NU,CUNIT)
        IMPLICIT NONE
        INTEGER NCOL,NEL,NN,NU
        INTEGER COL(NEL)
        CHARACTER*(NN) CNAME
        CHARACTER*(NU) CUNIT
*NCOL   input        column number
Cf2py  intent(in) ncol 
*NEL    input        number of column elements
Cf2py  intent(in) nel
*COL    input        column array (logical passed as integer from R)
Cf2py  intent(in) col
*NN     input        number of characters in column name
Cf2py  intent(in) nn
*CNAME  input        column name
Cf2py  intent(in) cname
*NU     input        number of characters in column unit
Cf2py  intent(in) nu
*CUNIT  input        column unit
Cf2py  intent(in) cunit
*-Author: Dick Willingale 2014-Dec-28
        INCLUDE 'QR_COM'
        INTEGER IHTYPE,NR,NC,IW,I,J,IR,IC
        INTEGER NDAT,NREP,NWID
        CHARACTER*(10) KWORD
        INTEGER NMAX
        PARAMETER(NMAX=1000)
        LOGICAL IL(NMAX)
C
        IF(ISTAT.NE.0) RETURN
C
        IF(IFITS.EQ.0) THEN
                WRITE(*,*) 'QR_FITSPCOLL - fits file not open'
                ISTAT=1
                RETURN
        ENDIF
C Get type of current extension
        CALL FTGHDT(IFITS,IHTYPE,ISTAT)
C Check extension type is binary table
        IF(IHTYPE.NE.2) THEN
                WRITE(*,*)
     +         'QR_FITSPCOLL - current extension not binary table'
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
            WRITE(*,*) 'QR_FITSPCOLL - ncol',NCOL,'too large >',NC
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
C Get column info
        CALL FTGTCL(IFITS,NCOL,NDAT,NREP,NWID,ISTAT)
        IF(NREP.GT.NMAX) THEN
                WRITE(*,*) 'QR_FITSPCOLL - repeat',NREP,
     +                'too large >',NMAX
                ISTAT=1
                RETURN
        ENDIF
C Put array into column
        DO I=1,NEL,NREP
                DO J=1,NREP
                        IF(COL(I+J-1).EQ.1) THEN
                                IL(J)=.TRUE.
                        ELSE
                                IL(J)=.FALSE.
                        ENDIF
                ENDDO
                CALL FTPCLL(IFITS,NCOL,I,1,NREP,IL,ISTAT)
        ENDDO
        END
