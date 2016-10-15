*+SRT_QURT      Real roots of a quadratic
        SUBROUTINE SRT_QURT(A,B,C,IROOT,R1,R2)
        DOUBLE PRECISION A,B,C,R1,R2
        INTEGER IROOT
*A,B,C  input   quadratic coefficients
*IROOT  output  number of real roots
*R1,R2  output  real roots
*-Author Dick Willingale 1996-Nov-19
        DOUBLE PRECISION SR
C
        IF(A.EQ.0.0) THEN
                IF(B.NE.0.0) THEN
                        IROOT=1
                        R1=-C/B
                ELSE
                        IROOT=0
                ENDIF
        ELSE
                SR=B**2-4.0*A*C
                IF(SR.LT.0) THEN
                        IROOT=0
                ELSEIF(SR.EQ.0.0) THEN
                        IROOT=1
                        R1=-B/(2.0*A)
                ELSE
                        IROOT=2        
                        SR=SQRT(SR)
                        R1=(-B+SR)/(2.0*A)
                        R2=(-B-SR)/(2.0*A)
                ENDIF
        ENDIF
        END

