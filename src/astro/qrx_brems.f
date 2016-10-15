*+QRX_BREMS        Bremsstrahlung with Gaunt factor
        SUBROUTINE QRX_BREMS(NE,E,T,PH)
        IMPLICIT NONE
        INTEGER NE
        DOUBLE PRECISION E(NE),T,PH(NE)
*NE              number of energies
*E               Energies keV
*T               Temperature keV
*PH              Bremsstrahlung continuum at energy E (photons/keV)
* Normalised to 1 ph/kev at 1 keV
*-Author Dick Willingale 2014-July-22
        INTEGER I
        DOUBLE PRECISION QRX_GAUNT,NORM
        IF(T.GT.1.E-15) THEN
                NORM=1.0/EXP(-1.0/T)/QRX_GAUNT(1.0D0,T)
                DO I=1,NE
                        PH(I)=NORM*EXP(-E(I)/T)*QRX_GAUNT(E(I),T)/E(I)
                ENDDO
        ELSE
                DO I=1,NE
                        PH(I)=0.
                ENDDO
        ENDIF
        END
*+QRX_GAUNT        Born approximation to Gaunt factor
        FUNCTION QRX_GAUNT(E,T)
        IMPLICIT NONE
        DOUBLE PRECISION QRX_GAUNT,E,T
*E                Energy keV
*T                Temperature keV
*SPF_GAUNT        Gaunt factor at energy E
*-Author Dick Willingale 2014-July-22
C From Gordon Stewart's code May 1984
        DOUBLE PRECISION X,AIO,Y,T2,Z,Y2,BO,FAC
        X=(E/T)*0.5
        IF(X.LT.2.) THEN
                T2=(X/3.75)**2
                AIO=(((((.0045813*T2+.0360768)*T2+
     +                .2659732)*T2+1.2067492)*T2+
     +                3.0899424)*T2+3.5156229)*T2+1.0
                Y=X*0.5
                Z=LOG(Y)
                Y2=Y**2
                BO=-AIO*Z+(((((.0000074*Y2+.0001075)*Y2+.00262698)*Y2+
     +                .0348859)*Y2+.23069756)*Y2+.4227842)*Y2-.57721566
                QRX_GAUNT=.551329*EXP(X)*BO
        ELSE
                Y=2./X
                BO=(((((.00053208*Y-.0025154)*Y+
     +                .00587872)*Y-.01062446)*Y+
     +                .02189568)*Y-.07832358)*Y+1.25331414
                FAC=BO/SQRT(X)
                QRX_GAUNT=.551329*FAC
        ENDIF
        END
