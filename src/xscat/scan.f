*+SCAN scans string for match with any of specified set of characters
      INTEGER FUNCTION SCAN(STRING, SET, BACK)
      CHARACTER STRING*(*), SET*(*)	
      LOGICAL BACK
*STRING  input    Character string to scan
*SET     input    Set of characters to find
*BACK    input    If .TRUE. scans from right to left, else left to right.
*SCAN    function returns position in STRING where first match is found, else
* returns zero if STRING contains none of the characters listed in SET.
*-Author	Clive Page	1991-July-9
*Can be replaced by FORTRAN-90 intrinsic function SCAN when available
      INTEGER K
*
      IF(BACK) THEN
         DO K = LEN(STRING),1,-1
            IF(INDEX(SET,STRING(K:K)) .NE. 0) THEN
               SCAN = K
               GO TO 999
            END IF
         END DO
      ELSE
         DO K = 1,LEN(STRING)
            IF(INDEX(SET,STRING(K:K)) .NE. 0) THEN
               SCAN = K
               GO TO 999
            END IF
         END DO
      END IF
      SCAN = 0
999   CONTINUE
      END
