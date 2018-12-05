*+SRT_SLEFSLOT      Finds  Schmidt lobster eye slot
C        SUBROUTINE SRT_SLEFSLOT(DIR,XA,CAX,CAR,CVR,PC,HIT,PH,CNL,
C     +        ISTAT)
        SUBROUTINE SRT_SLEFSLOT(DIR,XA,CVR,PC,XX,YY,
     +        ALPHA,NSLOT,RAYLOST,ISTAT)
        IMPLICIT NONE
        DOUBLE PRECISION DIR(3),XA(3),CVR(3),PC(50)
        DOUBLE PRECISION XX,YY
        DOUBLE PRECISION ALPHA
        LOGICAL RAYLOST
        INTEGER ISTAT
*DIR      input   direction of ray +
*XA       input   origin of ray +
*CVR      input   centre of cylinder +
*PC       input   common parameters +
*RAYLOST  output  .TRUE. if ray does not hit aperture, i.e. if it is lost
*ALPHA    output  angular position of centre of hit slot
*XX       output  local X of centre of hit slot
*YY       output  local Y of centre of hit slot
*ISTAT    in/out  returned status
*
* Parameters used from common field PC:
* PC(1) is the radius of the sphere -> RCYL
* PC(3) is the focal length -> FL
* PC(12) -> NMIR     number of mirrors
* PC(13) -> EPSILON  angle pitch of mirrors
* PC(14) -> ZALPHA   angle position of mirror of index 0
*
* original code for spherical aperture Author Dick Willingale 1996-Nov-28
* For cylindrical aperture modified by Vladimir Tichy
*
* Fasting cure performed.
*
        DOUBLE PRECISION A,B,C,D,VP(3),VR(3),VY(3)
        DOUBLE PRECISION X,Y,R,DN
        DOUBLE PRECISION RCYL,FL
        INTEGER CI,CJ
        INTEGER K
        DOUBLE PRECISION PH(3),CNL(3)
        DOUBLE PRECISION DNSLOT, EPSILON, ZALPHA
        INTEGER NSLOT,NMIR
C PH  = intersection point
C CNL = normal to surface at intersection point
C
C#define VDEBUG
C        
        IF(ISTAT.NE.0) RETURN
C
C
C If calculation is terminated before it's finished, it is automatically
C assigned that the ray is lost for some reason.
C
        RAYLOST=.TRUE.
C
C parameters used from common field:
        RCYL=PC(1)
        FL=PC(3)
        NMIR=PC(12)
        EPSILON=PC(13)
        ZALPHA=PC(14)
C
c#ifdef VDEBUG
c        WRITE (*,*) 'Ray origin      ',XA(1),XA(2),XA(3)
c        WRITE (*,*) 'Ray direction   ',DIR(1),DIR(2),DIR(3)
c        WRITE (*,*) 'Aperture centre ',CVR(1),CVR(2),CVR(3)
c        WRITE (*,*) 'Aperture radius ',RCYL
c#endif
C
c relative position of ray origin to centre of aperture
        DO K=1,3
                VP(K)= XA(K)-CVR(K)
        ENDDO
C
C find cylinder orientation
        IF (FL.GT.0) THEN
C mirror orientation = , cylinder x,z
                CI=1
                CJ=3
        ELSE
C mirror orientation || , cylinder x,y
                CI=1
                CJ=2
        ENDIF
C
C coefficients of quadratic equation
        A=DIR(CI)*DIR(CI)+DIR(CJ)*DIR(CJ)
        B=2*(VP(CI)*DIR(CI)+VP(CJ)*DIR(CJ))
        C=VP(CI)*VP(CI)+VP(CJ)*VP(CJ)-RCYL*RCYL
        D=B*B-(4*A*C)
C
c#ifdef VDEBUG
c        WRITE (*,*) 'D=',D
c#endif
C
        IF (D.LT.0.)         THEN
                WRITE (*,*) 'Error: Discriminant in SRT_CLNDR < 0'
                RETURN
        ENDIF
C if the setup is right, this is the desired root (the smaller one):
        DN=(-B-SQRT(D))/(2*A)
C
C        DN=-B/(2*A)
C        WRITE (*,*) 'DNw=',DN
C        IF (D.GT.0.) THEN
C                D=SQRT(D)
C                WRITE (*,*) 'SQRT D=',D
C                DN=DN-D/(2*A)
C        ENDIF
C
C DN represents the length along the ray between the source and the intersection.
c#ifdef VDEBUG
c        WRITE (*,*) 'DN=',DN
c#endif
        IF (DN.LT.0.) RETURN
c#ifdef VDEBUG
c        WRITE (*,*) '*HIT*'
c#endif
c 
C Find intersection point PH
        DO K=1,3
                PH(K)=XA(K)+DIR(K)*DN
        ENDDO
C
C angle position of hit in cylindrical coordinates
        ALPHA=ASIN((PH(CJ)-CVR(CJ))/RCYL)
c
C find slot index
        IF (ALPHA.GE.ZALPHA) THEN
                NSLOT=INT((ALPHA-ZALPHA)/EPSILON)
        ELSE
                NSLOT=-1
        ENDIF
C convert to double
        DNSLOT=NSLOT
c
C Now ALPHA is angle position of centre of slot
        ALPHA=(DNSLOT+0.5)*EPSILON+ZALPHA
c find local X,Y position of SLOT CENTRE
        IF (FL.GT.0) THEN
C mirror orientation = , cylinder x,z
                XX=0.0
                YY=RCYL*SIN(ALPHA)
        ELSE
C mirror orientation || , cylinder x,y
                XX=RCYL*SIN(ALPHA)
                YY=0.0
        ENDIF
C
c#ifdef VDEBUG
c        WRITE (*,*) 'Alpha',ALPHA
c        WRITE (*,*) '[',PH(1),PH(2),PH(3),PH(2)-XX,PH(3)-YY,']'
c#endif
C  set RAYLOST to indicate that the ray is to be traced 
        RAYLOST=.FALSE.
c
        END
