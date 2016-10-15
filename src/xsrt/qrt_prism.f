*+QRT_PRISM Set up prism
        SUBROUTINE QRT_PRISM(ID,IDF,IQ,ANML,ARFX,APOS,RAP,D1,D2,
     +  REFIND,THICK)
*ID        input        1 small angle, 2 right-angle
*IDF       input        deformation index
*IQ        input        quality index
*ANML      input        entrance surface normal 
*ARFX      input        reference axis
*APOS      input        reference position
*RAP       input        aperture radius
*D1        input        small angle radians on entry side (if ID=1)
*D2        input        small angle radians on exit side (if ID=1)
*REFIND    input        refractive index of material or n2/n1
*THICK     input        thickness
Cf2py  intent(in) ID,IDF,IQ,ANML,ARFX,APOS,RAP,D1,D2,REFIND,THICK
        IMPLICIT NONE
        INTEGER ID,IDF,IQ
        DOUBLE PRECISION ANML(3),ARFX(3),APOS(3),RAP,D1,D2,REFIND,THICK
*-Author Dick Willingale 20l2-Jun-28
        INTEGER IT,IDEF(2),NP,MAXRS,IQQ
        PARAMETER (MAXRS=20)
        DOUBLE PRECISION PL(MAXRS),TMIN1,TMIN2,S1,S2,C1,C2
        INCLUDE 'QR_COM'
C
        IF(ISTAT.NE.0) RETURN
C
        IDEF(1)=IDF
C set surface quality for internal and external refracting surface
        CALL SRT_SETT(IQ,3,1.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
     +        REFIND,0.0D0,0,0.0D0,0.0D0,ISTAT)
        CALL SRT_SETT(IQ+1,3,1.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
     +        1.D0/REFIND,0.0D0,0,0.0D0,0.0D0,ISTAT)
        IF(ID.EQ.2) THEN
C right-angled prism
                IDEF(2)=0
C specify first surface
                PL(1)=ANML(1)
                PL(2)=ANML(2)
                PL(3)=ANML(3)
                PL(4)=ARFX(1)
                PL(5)=ARFX(2)
                PL(6)=ARFX(3)
                PL(7)=APOS(1)+ANML(1)*(-THICK)
                PL(8)=APOS(2)+ANML(2)*(-THICK)
                PL(9)=APOS(3)+ANML(3)*(-THICK)
                IQQ=IQ+1
                PL(10)=-RAP
                PL(11)=-RAP
                PL(12)=RAP
                PL(13)=RAP
                IT=5
                NP=13
                CALL SRT_SETF(0,IT,NP,PL,IDEF,IQQ,-1,-1,ISTAT)
C specify second surface
                PL(1)=(ANML(1)+ARFX(1))*0.5
                PL(2)=(ANML(2)+ARFX(2))*0.5
                PL(3)=(ANML(3)+ARFX(3))*0.5
                PL(4)=(ARFX(1)-ANML(1))*0.5
                PL(5)=(ARFX(2)-ANML(2))*0.5
                PL(6)=(ARFX(3)-ANML(3))*0.5
                PL(7)=APOS(1)
                PL(8)=APOS(2)
                PL(9)=APOS(3)
                IQQ=IQ
                PL(10)=-RAP
                PL(11)=-RAP
                PL(12)=RAP
                PL(13)=RAP
                IT=5
                NP=13
                CALL SRT_SETF(0,IT,NP,PL,IDEF,IQQ,-1,-1,ISTAT)
C Specify 3rd surface
                PL(1)=-ARFX(1)
                PL(2)=-ARFX(2)
                PL(3)=-ARFX(3)
                PL(4)=ANML(1)
                PL(5)=ANML(2)
                PL(6)=ANML(3)
                PL(7)=APOS(1)+ARFX(1)*(-THICK)
                PL(8)=APOS(2)+ARFX(2)*(-THICK)
                PL(9)=APOS(3)+ARFX(3)*(-THICK)
                IQQ=IQ+1
                PL(10)=-RAP
                PL(11)=-RAP
                PL(12)=RAP
                PL(13)=RAP
                IT=5
                NP=13
                CALL SRT_SETF(0,IT,NP,PL,IDEF,IQQ,-1,-1,ISTAT)
            ELSE
C small angle prism
                IDEF(2)=0
C Check thickness
                S1=SIN(D1)
                S2=SIN(D2)
                TMIN1=THICK*S1/(S1+S2)
                TMIN2=THICK*S2/(S1+S2)
                IF(TMIN1.LT.RAP*S1.OR.TMIN2.LT.RAP*S2) THEN
                    WRITE(*,*) 'QRT_PRISM error - too thin for angles'
                    ISTAT=1
                    RETURN
                ENDIF
C Specify first surface
                C1=COS(D1)
                C2=COS(D2)
                PL(1)=-ANML(1)*C1+ARFX(1)*S1
                PL(2)=-ANML(2)*C1+ARFX(2)*S1
                PL(3)=-ANML(3)*C1+ARFX(3)*S1
                PL(4)=-ANML(1)*S1-ARFX(1)*C1
                PL(5)=-ANML(2)*S1-ARFX(2)*C1
                PL(6)=-ANML(3)*S1-ARFX(3)*C1
                PL(7)=APOS(1)+ANML(1)*(-TMIN1)
                PL(8)=APOS(2)+ANML(2)*(-TMIN1)
                PL(9)=APOS(3)+ANML(3)*(-TMIN1)
                IQQ=IQ
                PL(10)=-RAP
                PL(11)=-RAP
                PL(12)=RAP
                PL(13)=RAP
                IT=5
                NP=13
                CALL SRT_SETF(0,IT,NP,PL,IDEF,IQQ,-1,-1,ISTAT)
C specify second surface
                PL(1)=ANML(1)*C2+ARFX(1)*S2
                PL(2)=ANML(2)*C2+ARFX(2)*S2
                PL(3)=ANML(3)*C2+ARFX(3)*S2
                PL(4)=-ANML(1)*S2+ARFX(1)*C2
                PL(5)=-ANML(2)*S2+ARFX(2)*C2
                PL(6)=-ANML(3)*S2+ARFX(3)*C2
                PL(7)=APOS(1)+ANML(1)*(TMIN2)
                PL(8)=APOS(2)+ANML(2)*(TMIN2)
                PL(9)=APOS(3)+ANML(3)*(TMIN2)
                IQQ=IQ
                PL(10)=-RAP
                PL(11)=-RAP
                PL(12)=RAP
                PL(13)=RAP
                IT=5
                NP=13
                CALL SRT_SETF(0,IT,NP,PL,IDEF,IQQ,-1,-1,ISTAT)
        ENDIF
        END
