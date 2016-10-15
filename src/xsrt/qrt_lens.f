*+QRT_LENS Set up lens
        SUBROUTINE QRT_LENS(ID,IDF,IQ,ANML,ARFX,APOS,RAP,R1,R2,
     +  REFIND,THICK)
        IMPLICIT NONE
        INTEGER ID,IDF,IQ
        DOUBLE PRECISION ANML(3),ARFX(3),APOS(3),RAP,R1,R2,REFIND,THICK
*ID      input        lens type 1 spherical, 2 cylindrical
*IDF     input        deformation index
*IQ      input        surface quality index
*ANML    input        surface normal 
*ARFX    input        surface reference axis
*APOS    input        surface reference position
*RAP     input        radius of aperture
*R1,R2   input        radii of curvature of lens surfaces
*REFIND  input        refractive index of lens material (or n2/n1)
*THICK   input        lens thickness
Cf2py  intent(in) ID,IDF,IQ,ANML,ARFX,APOS,RAP,R1,R2,REFIND,THICK
*-Author Dick Willingale 2012-Jun-28
        INTEGER IT,IDEF(2),NP,MAXRS,IQQ
        PARAMETER (MAXRS=20)
        DOUBLE PRECISION PL(MAXRS),TMIN1,TMIN2
        INCLUDE 'QR_COM'
C
        IF(ISTAT.NE.0) RETURN
C
        IDEF(1)=IDF
        THICK=THICK*0.5
C set surface quality for refracting surface
        CALL SRT_SETT(IQ,3,1.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
     +        REFIND,0.0D0,0,0.0D0,0.0D0,ISTAT)
        CALL SRT_SETT(IQ+1,3,1.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
     +        1.D0/REFIND,0.0D0,0,0.0D0,0.0D0,ISTAT)
        IF(ID.EQ.1) THEN
C Spherical lens, radial limits
                IDEF(2)=0
C Check thickness
                IF(ABS(R1).GT.RAP) THEN
                        TMIN1=ABS(R1)-SQRT(R1**2-RAP**2)
                ELSE
                        TMIN1=ABS(R1)
                ENDIF
                IF(ABS(R2).GT.RAP) THEN
                        TMIN2=ABS(R2)-SQRT(R2**2-RAP**2)
                ELSE
                        TMIN2=ABS(R2)
                ENDIF
                IF(THICK.LT.TMIN1.OR.THICK.LT.TMIN2) THEN
                     WRITE(*,*) 'QRT_LENS error -  too thin for radii'
                     ISTAT=1
                     RETURN
                ENDIF
C specify first surface
                IF(R1.GE.0.0) THEN
                        PL(1)=-ANML(1)
                        PL(2)=-ANML(2)
                        PL(3)=-ANML(3)
                        PL(4)=-ARFX(1)
                        PL(5)=-ARFX(2)
                        PL(6)=-ARFX(3)
                        PL(7)=APOS(1)+ANML(1)*(R1-THICK)
                        PL(8)=APOS(2)+ANML(2)*(R1-THICK)
                        PL(9)=APOS(3)+ANML(3)*(R1-THICK)
                        IQQ=IQ
                ELSE
                        PL(1)=ANML(1)
                        PL(2)=ANML(2)
                        PL(3)=ANML(3)
                        PL(4)=ARFX(1)
                        PL(5)=ARFX(2)
                        PL(6)=ARFX(3)
                        PL(7)=APOS(1)+ANML(1)*(R1-THICK+TMIN1)
                        PL(8)=APOS(2)+ANML(2)*(R1-THICK+TMIN1)
                        PL(9)=APOS(3)+ANML(3)*(R1-THICK+TMIN1)
                        IQQ=IQ+1
                ENDIF
                IF(R1.NE.0.0) THEN
                        PL(10)=ABS(R1)
                        PL(11)=0.0
                        PL(12)=RAP
                        IT=15
                        NP=12
                ELSE
                        PL(10)=0.0
                        PL(11)=RAP
                        IT=6
                        NP=11
                ENDIF
                CALL SRT_SETF(0,IT,NP,PL,IDEF,IQQ,-1,-1,ISTAT)
C specify second surface
                IF(R2.GT.0.0) THEN
                        PL(1)=-ANML(1)
                        PL(2)=-ANML(2)
                        PL(3)=-ANML(3)
                        PL(4)=-ARFX(1)
                        PL(5)=-ARFX(2)
                        PL(6)=-ARFX(3)
                        PL(7)=APOS(1)+ANML(1)*(R2+THICK-TMIN2)
                        PL(8)=APOS(2)+ANML(2)*(R2+THICK-TMIN2)
                        PL(9)=APOS(3)+ANML(3)*(R2+THICK-TMIN2)
                        IQQ=IQ+1
                ELSE
                        PL(1)=ANML(1)
                        PL(2)=ANML(2)
                        PL(3)=ANML(3)
                        PL(4)=ARFX(1)
                        PL(5)=ARFX(2)
                        PL(6)=ARFX(3)
                        PL(7)=APOS(1)+ANML(1)*(R2+THICK)
                        PL(8)=APOS(2)+ANML(2)*(R2+THICK)
                        PL(9)=APOS(3)+ANML(3)*(R2+THICK)
                        IQQ=IQ
                ENDIF
                IF(R2.NE.0.0) THEN
                        PL(10)=ABS(R2)
                        PL(11)=0.0
                        PL(12)=RAP
                        IT=15
                        NP=12
                ELSE
                        PL(10)=0.0
                        PL(11)=RAP
                        IT=6
                        NP=11
                ENDIF
                CALL SRT_SETF(0,IT,NP,PL,IDEF,IQQ,-1,-1,ISTAT)
            ELSE
C Cylindrical lens, cartesian limits
                IDEF(2)=0
C Check thickness
                IF(R1.GT.RAP) THEN
                        TMIN1=R1-SQRT(R1**2-RAP**2)
                ELSE
                        TMIN1=R1
                ENDIF
                IF(R2.GT.RAP) THEN
                        TMIN2=R2-SQRT(R2**2-RAP**2)
                ELSE
                        TMIN2=R2
                ENDIF
                IF(THICK.LT.TMIN1.OR.THICK.LT.TMIN2) THEN
                    WRITE(*,*) 'QRT_LENS error -  too thin for radii'
                    ISTAT=1
                    RETURN
                ENDIF
C Specify first surface
                IF(R1.GT.0.0) THEN
                        PL(1)=-ARFX(1)
                        PL(2)=-ARFX(2)
                        PL(3)=-ARFX(3)
                        PL(4)=-ANML(1)
                        PL(5)=-ANML(2)
                        PL(6)=-ANML(3)
                        PL(7)=APOS(1)+ANML(1)*(R1-THICK)
                        PL(8)=APOS(2)+ANML(2)*(R1-THICK)
                        PL(9)=APOS(3)+ANML(3)*(R1-THICK)
                        IQQ=IQ
                ELSEIF(R1.LT.0.0) THEN
                        PL(1)=ARFX(1)
                        PL(2)=ARFX(2)
                        PL(3)=ARFX(3)
                        PL(4)=ANML(1)
                        PL(5)=ANML(2)
                        PL(6)=ANML(3)
                        PL(7)=APOS(1)+ANML(1)*(R1-THICK+TMIN1)
                        PL(8)=APOS(2)+ANML(2)*(R1-THICK+TMIN1)
                        PL(9)=APOS(3)+ANML(3)*(R1-THICK+TMIN1)
                        IQQ=IQ+1
                ELSE
                        PL(1)=-ANML(1)
                        PL(2)=-ANML(2)
                        PL(3)=-ANML(3)
                        PL(4)=-ARFX(1)
                        PL(5)=-ARFX(2)
                        PL(6)=-ARFX(3)
                        PL(7)=APOS(1)+ANML(1)*(-THICK)
                        PL(8)=APOS(2)+ANML(2)*(-THICK)
                        PL(9)=APOS(3)+ANML(3)*(-THICK)
                        IQQ=IQ
                ENDIF
                IF(R1.NE.0.0) THEN
                        PL(10)=0.0
                        PL(11)=0.0
                        PL(12)=R1**2
                        PL(13)=-RAP
                        PL(14)=-ASIN(RAP/ABS(R1))
                        PL(15)=RAP
                        PL(16)=ASIN(RAP/ABS(R1))
                        IT=19
                        NP=16
                ELSE
                        PL(10)=-RAP
                        PL(11)=-RAP
                        PL(12)=RAP
                        PL(13)=RAP
                        IT=5
                        NP=13
                ENDIF
                CALL SRT_SETF(0,IT,NP,PL,IDEF,IQQ,-1,-1,ISTAT)
C specify second surface
                IF(R2.GT.0.0) THEN
                        PL(1)=-ARFX(1)
                        PL(2)=-ARFX(2)
                        PL(3)=-ARFX(3)
                        PL(4)=-ANML(1)
                        PL(5)=-ANML(2)
                        PL(6)=-ANML(3)
                        PL(7)=APOS(1)+ANML(1)*(R2+THICK-TMIN2)
                        PL(8)=APOS(2)+ANML(2)*(R2+THICK-TMIN2)
                        PL(9)=APOS(3)+ANML(3)*(R2+THICK-TMIN2)
                        IQQ=IQ+1
                ELSEIF(R2.LT.0.0) THEN
                        PL(1)=ARFX(1)
                        PL(2)=ARFX(2)
                        PL(3)=ARFX(3)
                        PL(4)=ANML(1)
                        PL(5)=ANML(2)
                        PL(6)=ANML(3)
                        PL(7)=APOS(1)+ANML(1)*(R2+THICK)
                        PL(8)=APOS(2)+ANML(2)*(R2+THICK)
                        PL(9)=APOS(3)+ANML(3)*(R2+THICK)
                        IQQ=IQ
                ELSE
                        PL(1)=ANML(1)
                        PL(2)=ANML(2)
                        PL(3)=ANML(3)
                        PL(4)=ARFX(1)
                        PL(5)=ARFX(2)
                        PL(6)=ARFX(3)
                        PL(7)=APOS(1)+ANML(1)*(THICK)
                        PL(8)=APOS(2)+ANML(2)*(THICK)
                        PL(9)=APOS(3)+ANML(3)*(THICK)
                        IQQ=IQ
                ENDIF
                IF(R2.NE.0.0) THEN
                        PL(10)=0.0
                        PL(11)=0.0
                        PL(12)=R2**2
                        PL(13)=-RAP
                        PL(14)=-ASIN(RAP/ABS(R2))
                        PL(15)=RAP
                        PL(16)=ASIN(RAP/ABS(R2))
                        IT=19
                        NP=16
                ELSE
                        PL(10)=-RAP
                        PL(11)=-RAP
                        PL(12)=RAP
                        PL(13)=RAP
                        IT=5
                        NP=13
                ENDIF
                CALL SRT_SETF(0,IT,NP,PL,IDEF,IQQ,-1,-1,ISTAT)
        ENDIF
        END
