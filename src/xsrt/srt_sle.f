*+SRT_SLE    Dynamically set up a slot for Schmidt Lobster, derived from  SRT_KBS for Dynamically set up a rectangular slot for K-B stack
        SUBROUTINE SRT_SLE(DIR,ORG,HPOS,CAX,CAR,CVR,PC,IDEF,OHIT,ISTAT)
        IMPLICIT NONE
        DOUBLE PRECISION DIR(3),ORG,HPOS(3),CAX(3),CAR(3),CVR(3),PC(50)
        INTEGER IDEF(2),ISTAT
*DIR        input        direction of ray
*ORG        input   origin of ray
*HPOS        input        position of hit on aperture plane reference surface
*CAX    input   axis of optic
*CAR    input   reference direction in aperture
*CVR    input   centre of curvature of front aperture 
*PC     input   remaining parameters
*IDEF   input   deformation indices
*OHIT        output  TRUE if the ray is lost
*ISTAT  in/out  returned status
* Parameters for slots are:
* PC(1) radius of spherical aperture surface if 1st stack 2*FLEN+PLMAX
*                                             if 2nd stack 2*FLEN
* PC(2) outer radius of aperture
* PC(3) focal length FLEN (-ve for 2nd stack)
* PC(4) pitch of slots
* PC(5) wall thickness of slots
* PC(6) minimum axial slot length (PLMIN)
* PC(7) maximum axial slot length (PLMAX)
* PC(8) local x of ray intersection
* PC(9) local y of ray intersection
* PC(10) this surface index
* PC(11) surface quality
* PC(12) number of mirrors
* PC(13) angle distance between adjacent mirrors
* PC(14) angle position of first mirror
*
* Note: -ve focal length indicates 2nd stack
*
*-Author of original code for KB: Dick Willingale 2008-Nov-20
* Modification for SLE: Vladimir Tichy
*
C        INCLUDE 'SRT_COM'
*
        DOUBLE PRECISION PI
        PARAMETER (PI=3.14159265359)
        DOUBLE PRECISION VY(3),X,Y,XP(3),YP(3),ZP(3),POS(3)
        DOUBLE PRECISION XC,YC,MCAX(3),DELT,RXY
        DOUBLE PRECISION PP(20),HX,HY,HL,THETA,PHI,RCYL
        DOUBLE PRECISION DX,DY,DZ,DB,DC,DDX,DDY,RAN(2)
        DOUBLE PRECISION XX,FL,YY,PL,ELH,ELW,ZS,SHWID,SHLEN
        INTEGER I,KSUR,IQ,NSLOT,ISTACK,ITRY
        DOUBLE PRECISION RM,PM,TM,WM,HM,CM,GM
        DOUBLE PRECISION VR(3), VN(3)
        INTEGER MODHIT,ICURV
        LOGICAL RAYLOST
        LOGICAL OHIT
        DOUBLE PRECISION ALPHA,ZALPHA,EPSILON
        INTEGER NMIR
        DOUBLE PRECISION DNSLOT
        DOUBLE PRECISION WXX, WYY
        DOUBLE PRECISION PITCH,LPITCH,WALL
        DOUBLE PRECISION CSDSTSQ
        INTEGER WNSLOT
C
C Modified deformation index for (R)eflecting plane
        INTEGER IDEFR(2)
C Modified deformation index for non-reflecting = (B)lack plane
        INTEGER IDEFB(2)
c        LOGICAL THRU
C
        OHIT=.TRUE.
c        THRU=.FALSE.
C
C#define VDEBUG
C#define XDEBUG
c#define WDEBUG
C
C Take parameters from global field. It is congregated here to
C easy check what is used and what is not.
        RCYL=PC(1)
        FL=PC(3)
        PITCH=PC(4)
        WALL=PC(5)
        NMIR=PC(12)
        EPSILON=PC(13)
        ZALPHA=PC(14)
C
        IDEFR(1)=IDEF(1)
        IDEFB(1)=0
        IDEFB(2)=0
C
c Set ISTACK and set FL to |FL|
        if(FL.GT.0.0) THEN
                ISTACK=1
        else
                ISTACK=2
                FL=-FL
        endif
c
C Find module index
        CALL SRT_SFINDMOD(RCYL,PC(8),PC(9),MODHIT,RM,PM,TM,WM,HM,PL,
     +        CM,GM,XX,YY)
C
        LPITCH=PITCH*(RCYL+PL/2)/RCYL
        CSDSTSQ=RCYL*RCYL-LPITCH*LPITCH/4.0
c
C for debugging:
Cifdef VDEBUG
C        WRITE(*,*) 'ISTACK:',ISTACK,' RCYL:',RCYL
C        WRITE(*,*) 'CAX   :',CAX(1),',',CAX(2),' ',CAX(3)
C        WRITE(*,*) 'CAR   :',CAR(1),',',CAR(2),' ',CAR(3)
C        WRITE(*,*) 'CVR   :',CVR(1),',',CVR(2),' ',CVR(3)
C        DO I=1,14
C                WRITE(*,*) 'PC(',I,')=',PC(I)
C        END DO
C#endif
C
C#ifdef WDEBUG
C        WRITE(*,*) 'MODHIT',MODHIT
C        WRITE(*,*) 'RM',RM
c        WRITE(*,*) 'PM',PM
c        WRITE(*,*) 'TM',TM
c        WRITE(*,*) 'WM',WM
c        WRITE(*,*) 'HM',HM
c        WRITE(*,*) 'PL',PL
c        WRITE(*,*) 'CM',CM
c        WRITE(*,*) 'GM',GM
c        WRITE(*,*) 'XX',XX
c        WRITE(*,*) 'YY',YY
C        WRITE(*,*) 'MODHIT',MODHIT
C        WRITE(*,*) 'RM',RM,rlist(POKUSIND)
C        WRITE(*,*) 'PM',PM,philist(POKUSIND)
C        WRITE(*,*) 'TM',TM,thlist(POKUSIND)
C        WRITE(*,*) 'WM',WM,wlist(POKUSIND)
C        WRITE(*,*) 'HM',HM,hlist(POKUSIND)
C        WRITE(*,*) 'PL',PL,llist(POKUSIND)
C        WRITE(*,*) 'CM',CM,clist(POKUSIND)
C        WRITE(*,*) 'GM',GM,glist(POKUSIND)
C        WRITE(*,*) 'XX',XX
C        WRITE(*,*) 'YY',YY
c#endif
C
        IF (MODHIT.LE.0) RETURN
        CALL SRT_SLEFSLOT(DIR,ORG,CVR,PC,WXX,WYY,ALPHA,NSLOT,RAYLOST,ISTAT)
        IF(RAYLOST) RETURN
c                
        THETA=PM+TM
        ELH=HM*0.5
        ELW=WM*0.5
C Centre of slot  is in local module aperture coordinates WXX,WYY
C        NSLOT=INT((ALPHA-ZALPHA)/EPSILON)
        IF (NSLOT.LT.0) NSLOT=-1
c Warning: mirrors indices are 1,2,3,...,2*NMIR
        IDEFR(2)=2*NSLOT+2
        IF (NSLOT.GE.(NMIR-1)) THEN
                NSLOT=-2
                IDEFR(2)=2*(NMIR-1)+1
        ENDIF
C        IF (NSLOT.GT.NMIR) RETURN
C
c        WRITE(*,*) 'SRT_SLE DEBUG NSLOT=',NSLOT
c
        DNSLOT=NSLOT
C
C Now ALPHA is angle position of centre of slot
C        ALPHA=(DNSLOT+0.5)*EPSILON+ZALPHA
c        WRITE(*,*) '[',NSLOT,']'
        IF(ISTACK.EQ.1) THEN
C orientation =
C                WYY=SIN(ALPHA)*RCYL
C                WXX=0.0
                SHLEN=ELW
                SHWID=ELH
        ELSE
C orientation ||
C                WXX=SIN(ALPHA)*RCYL
C                WYY=0.0
                SHLEN=ELH
                SHWID=ELW
        ENDIF
c
C    
        IF (NSLOT.LT.0) THEN
C        the ray goes between a mirror and housing
                XX=0.0
                YY=0.0
        ELSE
                XX=WXX
                YY=WYY
C Set slot width
C                SHWID=(PC(4)-PC(5))*0.5
                SHWID=LPITCH/2.
        ENDIF
C
C Calculate position of module centre in full aperture coordinates
        XC=RM*COS(PM)
        YC=RM*SIN(PM)
C Transform back to full aperture coordinates
        X=XX*COS(THETA)-YY*SIN(THETA)+XC
        Y=XX*SIN(THETA)+YY*COS(THETA)+YC
C Set slot rotation for 2nd stack
        IF(ISTACK.NE.1) THEN
                THETA=THETA-PI*0.5
        ENDIF
c
C X,Y are centre of slot in local full aperture coordinates
C Find Y reference axis at local origin on aperture
        CALL SRT_VCRS(CAX,CAR,VY)
        CALL SRT_VNRM(VY,ISTAT)
C Calculate position of centre of slot aperture
        DO I=1,3
C                POS(I)=CVR(I)+CAX(I)*SQRT(RCYL**2-X**2-Y**2)
                POS(I)=CVR(I)+CAX(I)*SQRT(CSDSTSQ-X*X-Y*Y)
                POS(I)=POS(I)+CAR(I)*X+VY(I)*Y
        ENDDO
C 
c#ifdef VDEBUG
c        WRITE(*,*) 'POS   :',POS(1),',',POS(2),' ',POS(3)
c#endif
C
c This is to place, where deformations can be applied
c
        DX=0.0
        DY=0.0
        DB=0.0
        DC=0.0
c
C Find direction of slot axis
        DO I=1,3
                ZP(I)=POS(I)-CVR(I)-DX*CAR(I)-DY*VY(I)-DZ*CAX(I)
        ENDDO
        CALL SRT_VNRM(ZP,ISTAT)
C Now set up axes tangential to spherical surface at slot centre
        CALL SRT_VCRS(ZP,CAR,YP)
        CALL SRT_VNRM(YP,ISTAT)
        CALL SRT_VCRS(YP,ZP,XP)
C Rotate reference axes to align with module
C Include out-of-plane figure errors by including slot rotation error
        THETA=THETA+DC
        DO I=1,3
                XP(I)=COS(THETA)*XP(I)+SIN(THETA)*YP(I)
        ENDDO
        CALL SRT_VCRS(ZP,XP,YP)
C Include in-plan8e figure errors by shifting axis in YP direction
        DO I=1,3
                ZP(I)=ZP(I)*SQRT(1.0-DB**2)+YP(I)*DB
        ENDDO
        CALL SRT_VNRM(ZP,ISTAT)
        CALL SRT_VCRS(ZP,XP,YP)
C Set dimensions of slot aperture
        HX=SHLEN
        HY=SHWID
C Set vector parameters for slot aperture
        PP(1)=ZP(1)
        PP(2)=ZP(2)
        PP(3)=ZP(3)
        PP(4)=XP(1)
        PP(5)=XP(2)
        PP(6)=XP(3)
        PP(7)=POS(1)
        PP(8)=POS(2)
        PP(9)=POS(3)
C
        PP(10)=-HX
C        PP(11)=-HY+WALL/2
        PP(11)=-HY
        PP(12)=HX
C        PP(13)=HY-WALL/2
        PP(13)=HY
C Find surface index of slot aperture and set that surface element
        KSUR=PC(10)+1
        CALL SRT_SETF(KSUR,1,13,PP,0,0,0,KSUR+1,ISTAT)
C set half axial length of slot
        HL=PL*0.5
c
c Zahada: kdyz se prida tato radka, zmeni se chovani
c        WRITE(*,*) 'HX HY HL PL:',HX,HY,HL,PL
c
C Get surface quality of slot walls
        IQ=INT(PC(11))
C local x is along slot axis
C
c#ifdef WDEBUG
c        WRITE(*,*) 'XP:',XP(1),',',XP(2),' ',XP(3)
c        WRITE(*,*) 'YP:',YP(1),',',YP(2),' ',YP(3)
c        WRITE(*,*) 'ZP:',ZP(1),',',ZP(2),' ',ZP(3)
c        WRITE(*,*)
c#endif
C
C <1> absorbing wall of slot 
        PP(1)=XP(1)
        PP(2)=XP(2)
        PP(3)=XP(3)
        PP(4)=ZP(1)
        PP(5)=ZP(2)
        PP(6)=ZP(3)
        PP(7)=POS(1)-ZP(1)*HL-XP(1)*HX
        PP(8)=POS(2)-ZP(2)*HL-XP(2)*HX
        PP(9)=POS(3)-ZP(3)*HL-XP(3)*HX
        PP(10)=-HL
        PP(11)=-HY
        PP(12)=HL
        PP(13)=HY
        CALL SRT_SETF(KSUR+1,5,13,PP,IDEFB,0,KSUR+1,-1,ISTAT)
C
c#ifdef XDEBUG
c        WRITE (*,*),'wALL 1'
c        DO I=1,13
c                WRITE(*,*) 'PP(',I,')=',PP(I)
c        END DO
c#endif
C
cC <3> absorbing wall of slot 
        PP(1)=-XP(1)
        PP(2)=-XP(2)
        PP(3)=-XP(3)
        PP(4)=ZP(1)
        PP(5)=ZP(2)
        PP(6)=ZP(3)
        PP(7)=POS(1)-ZP(1)*HL+XP(1)*HX
        PP(8)=POS(2)-ZP(2)*HL+XP(2)*HX
        PP(9)=POS(3)-ZP(3)*HL+XP(3)*HX
        PP(10)=-HL
        PP(11)=-HY
        PP(12)=HL
        PP(13)=HY
        CALL SRT_SETF(KSUR+3,5,13,PP,IDEFB,0,KSUR+1,-1,ISTAT)
C
c#ifdef XDEBUG
c        WRITE (*,*),'wALL 3'
c        DO I=1,13
c                WRITE(*,*) 'PP(',I,')=',PP(I)
c        END DO
c#endif
C
        IF (NSLOT.GE.0) THEN
C        normal slot
C <2> reflecting wall of slot, unchanged code
C        reference position on surface:
        PP(7)=POS(1)-ZP(1)*HL-YP(1)*HY
        PP(8)=POS(2)-ZP(2)*HL-YP(2)*HY
        PP(9)=POS(3)-ZP(3)*HL-YP(3)*HY
C        limits:
        PP(10)=-HL
        PP(11)=-HX
        PP(12)=HL
        PP(13)=HX
C <2> reflecting wall of slot, new code
C add slopeness - reference axis; direction from centre of curvature
C to the reference position
        VR(1)=PP(7)-CVR(1)
        VR(2)=PP(8)-CVR(2)
        VR(3)=PP(9)-CVR(3)
        CALL SRT_VNRM(VR,ISTAT)
        PP(4)=VR(1)
        PP(5)=VR(2)
        PP(6)=VR(3)
C add slopeness - normal
        CALL SRT_VCRS(VR,XP,VN)
        CALL SRT_VNRM(VN,ISTAT)
        PP(1)=-VN(1)
        PP(2)=-VN(2)
        PP(3)=-VN(3)
C Add wall thickness (not exaxt) - necessary to do as the last step
        PP(7)=PP(7)+YP(1)*WALL/2
        PP(8)=PP(8)+YP(2)*WALL/2
        PP(9)=PP(9)+YP(3)*WALL/2
C reflecting:
        CALL SRT_SETF(KSUR+2,5,13,PP,IDEFR,IQ,KSUR+1,-1,ISTAT)
C blacked, for experiments:
C        CALL SRT_SETF(KSUR+2,5,13,PP,0,0,KSUR+1,-1,ISTAT)
c
C <4> reflecting wall of slot 
C
        IDEFR(2)=IDEFR(2)+1
C
C        reference position on surface:
        PP(7)=POS(1)-ZP(1)*HL+YP(1)*HY
        PP(8)=POS(2)-ZP(2)*HL+YP(2)*HY
        PP(9)=POS(3)-ZP(3)*HL+YP(3)*HY
C        limits:
        PP(10)=-HL
        PP(11)=-HX
        PP(12)=HL
        PP(13)=HX
C add slopeness - reference axis; direction from centre of curvature
C to the reference position
        VR(1)=PP(7)-CVR(1)
        VR(2)=PP(8)-CVR(2)
        VR(3)=PP(9)-CVR(3)
        CALL SRT_VNRM(VR,ISTAT)
        PP(4)=VR(1)
        PP(5)=VR(2)
        PP(6)=VR(3)
C add slopeness - normal
        CALL SRT_VCRS(VR,XP,VN)
        CALL SRT_VNRM(VN,ISTAT)
        PP(1)=VN(1)
        PP(2)=VN(2)
        PP(3)=VN(3)
C Add wall thickness (not exaxt) - necessary to do as the last step
        PP(7)=PP(7)-YP(1)*WALL/2
        PP(8)=PP(8)-YP(2)*WALL/2
        PP(9)=PP(9)-YP(3)*WALL/2
C reflecting:
        CALL SRT_SETF(KSUR+4,5,13,PP,IDEFR,IQ,KSUR+1,-1,ISTAT)
C blacked, for experiments:
C        CALL SRT_SETF(KSUR+4,5,13,PP,0,0,KSUR+1,-1,ISTAT)
C
c
        ELSE
c
c#ifdef XDEBUG
c        WRITE(*,*) 'HX HY HL PL:',HX,HY,HL,PL
c        WRITE(*,*) 'POS:',POS(1),',',POS(2),' ',POS(3)
c        WRITE(*,*) 'XP:',XP(1),',',XP(2),' ',XP(3)
c        WRITE(*,*) 'YP:',YP(1),',',YP(2),' ',YP(3)
c        WRITE(*,*) 'ZP:',ZP(1),',',ZP(2),' ',ZP(3)
c        WRITE(*,*)
c#endif
C        the ray goes between mirror and housing
c        RETURN
C <2> 
        PP(1)=YP(1)
        PP(2)=YP(2)
        PP(3)=YP(3)
        PP(4)=ZP(1)
        PP(5)=ZP(2)
        PP(6)=ZP(3)
C <2> 
C        reference position on surface:
        PP(7)=POS(1)-ZP(1)*HL-YP(1)*HY
        PP(8)=POS(2)-ZP(2)*HL-YP(2)*HY
        PP(9)=POS(3)-ZP(3)*HL-YP(3)*HY
C        limits:
        PP(10)=-HL
        PP(11)=-HX
        PP(12)=HL
        PP(13)=HX
C
c#ifdef XDEBUG
c        WRITE (*,*),'wALL 2'
c        DO I=1,13
c                WRITE(*,*) 'PP(',I,')=',PP(I)
c        END DO
c        CALL SRT_SETF(KSUR+2,32,13,PP,IDEFB,0,KSUR+1,-1,ISTAT)
c
c        WRITE(*,*) 'HX HY HL PL:',HX,HY,HL,PL
c        WRITE(*,*) 'POS:',POS(1),',',POS(2),' ',POS(3)
c        WRITE(*,*) 'XP:',XP(1),',',XP(2),' ',XP(3)
c        WRITE(*,*) 'YP:',YP(1),',',YP(2),' ',YP(3)
c        WRITE(*,*) 'ZP:',ZP(1),',',ZP(2),' ',ZP(3)
c        WRITE(*,*)
c#endif
C        the ray goes between mirror and housing
C <4> 
        PP(1)=-YP(1)
        PP(2)=-YP(2)
        PP(3)=-YP(3)
        PP(4)=ZP(1)
        PP(5)=ZP(2)
        PP(6)=ZP(3)
C        reference position on surface:
        PP(7)=POS(1)-ZP(1)*HL+YP(1)*HY
        PP(8)=POS(2)-ZP(2)*HL+YP(2)*HY
        PP(9)=POS(3)-ZP(3)*HL+YP(3)*HY
C        limits:
        PP(10)=-HL
        PP(11)=-HX
        PP(12)=HL
        PP(13)=HX
C
c#ifdef XDEBUG
c        WRITE (*,*),'wALL 4'
c        DO I=1,13
c                WRITE(*,*) 'PP(',I,')=',PP(I)
c        END DO
c
c        WRITE(*,*) 'HX HY HL PL:',HX,HY,HL,PL
c        WRITE(*,*) 'POS:',POS(1),',',POS(2),' ',POS(3)
c        WRITE(*,*) 'XP:',XP(1),',',XP(2),' ',XP(3)
c        WRITE(*,*) 'YP:',YP(1),',',YP(2),' ',YP(3)
c        WRITE(*,*) 'ZP:',ZP(1),',',ZP(2),' ',ZP(3)
c        WRITE(*,*)
c#endif
C        the ray goes between mirror and housing
c
        CALL SRT_SETF(KSUR+4,5,13,PP,IDEFB,0,KSUR+1,-1,ISTAT)
C
        ENDIF
c
C
        OHIT=.FALSE.
        END
