*+SX_FACTOR  Computes factorial for a positive integer (up to about 33)
      REAL FUNCTION SX_FACTOR(N)
      INTEGER N
*-Author  Clive Page  1993-OCT-27
      INTEGER I
*
      SX_FACTOR = 1.0
      DO I = 2,N
         SX_FACTOR = SX_FACTOR * REAL(I)
      END DO
      END

