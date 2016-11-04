      SUBROUTINE ERRMSG(CSTATS,CSUBRT,CMSGNM)
C Arguments:
      CHARACTER CMSGNM*(*),CSTATS*(*),CSUBRT*(*)
C
C Local variables:
      INTEGER IOERR
C
C External Subroutines:
      EXTERNAL GETSET
C***********************************************************************
C Given:
C       CSTATS = 'WARNING' or 'FATAL'
C       CSUBRT = name of subroutine
C       CMESGN = message
C Prints a warning message in a standardized way,
C and STOPs if CSTATS='FATAL'
C
C Copyright (C) 1993, B.T. Draine and P.J. Flatau
C This code is covered by the GNU General Public License.
C***********************************************************************
      CALL GETSET('GET','IOERR',IOERR)
      IF(CSTATS.EQ.'FATAL')THEN
         WRITE (IOERR,FMT=9000) CSUBRT
 9000    FORMAT(/,' >>>>> FATAL ERROR IN PROCEDURE: ',A)
         WRITE(IOERR,FMT=9010)CMSGNM
 9010    FORMAT(1X,A)
         WRITE(IOERR,FMT=9020)
 9020    FORMAT(' >>>>> EXECUTION ABORTED ')
         STOP
      ELSEIF(CSTATS.EQ.'WARNING')THEN
         WRITE(IOERR,FMT=9030)CSUBRT
 9030    FORMAT(/,' >>>>> WARNING IN PROCEDURE: ',A)
         WRITE(IOERR,FMT=9010)CMSGNM
 9040    FORMAT(1X,A)
      ENDIF
      RETURN
      END
