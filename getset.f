      SUBROUTINE GETSET(CWHAT,CVAR,INT)
C     .. Parameters ..
      INTEGER ISET
      PARAMETER (ISET=3)
C     ..
C     .. Scalar Arguments ..
      INTEGER INT
      CHARACTER CVAR* (*), CWHAT* (*)
C     ..
C     .. Local Scalars ..
      INTEGER I, IOERR
C     ..
C     .. Local Arrays ..
      INTEGER IARY(ISET)
      CHARACTER CARY(ISET)*10
C     ..
C     .. External Subroutines ..
      EXTERNAL ERRMSG
C     ..
C     .. Data statements ..
      DATA CARY/'MXNAT', 'IOOUT', 'IOERR'/
c      DATA IARY/0, 6, 6/
      DATA IARY/0, 0, 0/
      DATA IOERR/0/
C     ..
C PURPOSE: SETS AND GETS INTEGER VARIABLES ("BASKET PROCEDURE")
C FOR USE IN ERROR CHECKING OF PARAMETER STATEMENTS, AND/OR
C SETTING INPUT/OUTPUT INTEGERS. THUS THESE VARIABLES DON'T HAVE
C TO BE TRANSFERED AS FORMAL PARAMETERS OF THE PROCEDURE.
C EACH NEW CVAR HAS TO BE MANUALY SET IN THIS ROUTINE!
C CWHAT ---  'SET' OR  'GET'
C CVAR  --- CURRENTLY SET
C           'MXNAT'
C           'IOOUT'
C
C Copyright (C) 1993, B.T. Draine and P.J. Flatau
C This code is covered by the GNU General Public License.
C***********************************************************************
      IF (CWHAT .EQ. 'SET') THEN
          DO 10 I = 1, ISET
              IF (CVAR .EQ. CARY(I)) THEN
                  IARY(I) = INT
                  RETURN
              END IF
   10     CONTINUE
          CALL ERRMSG('FATAL','GETSET',' WRONG CVAR ')
          RETURN
C
      ELSE IF (CWHAT .EQ. 'GET') THEN
          DO 20 I = 1, ISET
              IF (CVAR .EQ. CARY(I)) THEN
                  INT = IARY(I)
                  RETURN
              END IF
   20     CONTINUE
          CALL ERRMSG('FATAL','GETSET',' WRONG CVAR ')
          RETURN
      ELSE
          CALL ERRMSG('FATAL','GETSET',' WRONG CWHAT ')
      END IF
      RETURN
      END
