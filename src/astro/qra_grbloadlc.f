*+QRA_GRBLOADLC Load GRB light curve data into common
        SUBROUTINE QRA_GRBLOADLC(NCBAT,NTBAT,TBAT,RBAT,ERBAT,
     +        NCXRT,NTXRT,TXRT,RXRT,ERXRT,MXRT)
        IMPLICIT NONE
        INTEGER NCBAT,NTBAT,NCXRT,NTXRT
        INTEGER MXRT(NTXRT)
        DOUBLE PRECISION TBAT(NTBAT),RBAT(NCBAT,NTBAT)
        DOUBLE PRECISION ERBAT(NCBAT,NTBAT)
        DOUBLE PRECISION TXRT(NTXRT),RXRT(NCXRT,NTXRT)
        DOUBLE PRECISION ERXRT(NCXRT,NTXRT)
*NCBAT        input        number of BAT energy channels
*NTBAT        input        number of BAT time samples
*TBAT        input        array BAT times
*RBAT        input        array BAT rates
*ERBAT        input        array BAT rate errors
*NCXRT        input        number of XRT energy channels
*NTXRT        input        number of XRT time samples
*TXRT        input        array XRT times
*RXRT        input        array XRT rates
*ERXRT        input        array XRT rate errors
*MXRT        input        mode values for XRT times (1 WT, 2 PC, 0 ignore)
*-Author Dick Willingale 2013-Jan-24
        include 'QR_COM'
        include 'GRBS_COM'
        double precision Tin
        integer I,J
C
        IF(ISTAT.NE.0) RETURN
C
        IF(NCBAT.GT.nc1) THEN
                write(*,*) 'QRA_GRBLOADLC too many BAT channels'
                ISTAT=1
                RETURN
        ENDIF
        IF(NCXRT.GT.nc2) THEN
                write(*,*) 'QRA_GRBLOADLC too many XRT channels'
                ISTAT=1
                RETURN
        ENDIF
        IF(NTBAT.GT.n1) THEN
                write(*,*) 'QRA_GRBLOADLC too many BAT times'
                ISTAT=1
                RETURN
        ENDIF
        IF(NTXRT.GT.n2) THEN
                write(*,*) 'QRA_GRBLOADLC too many XRT times'
                ISTAT=1
                RETURN
        ENDIF
        if(npars.gt.0) then
C Set injection time as start of 1st pulse
                Tin=pars(9)-pars(7)
        else
C Set injection time as trigger
                Tin=0.0
        endif
C Set linear fitting
        il=0
C Note if fitting in log flux then only accept +ve flux values and
C times after Tin
        mc1=NCBAT
        m1=NTBAT
        DO J=1,NTBAT
                if(il.eq.0.or.(il.eq.1.and.RBAT(1,J).gt.0.0
     +                .and.TBAT(J).gt.Tin)) then
                           x1(J)=TBAT(J)
                        DO I=1,mc1
                                   y1(I,J)=RBAT(I,J)
                                   ye1(I,J)=ERBAT(I,J)
                        ENDDO
                        ix1(J)=1
                else
                        ix1(J)=0
                   endif
        ENDDO
        mc2=NCXRT
        m2=NTXRT
        DO J=1,NTXRT
                IF(MXRT(J).NE.0) THEN
                    if(il.eq.0.or.(il.eq.1.and.RXRT(1,J).gt.0.0
     +                        .and.TXRT(J).gt.Tin)) then
                           x2(J)=TXRT(J)
                        mx2(J)=MXRT(J)
                        ix2(J)=1
                        DO I=1,mc2
                                   y2(I,J)=RXRT(I,J)
                                   ye2(I,J)=ERXRT(I,J)
                        ENDDO
                    endif
                ELSE
                    mx2(J)=0
                    ix2(J)=0
                ENDIF
        ENDDO
        END
