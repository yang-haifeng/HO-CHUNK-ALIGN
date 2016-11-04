      SUBROUTINE WRIMSG(CSUBRT,CMSGNM)
C Standard procedure for writing messages
C
C Copyright (C) 1993, B.T. Draine and P.J. Flatau
C This code is covered by the GNU General Public License.
C***********************************************************************
C     .. Scalar Arguments ..
      CHARACTER CMSGNM* (*), CSUBRT* (*)
C     ..
C     .. Local Scalars ..
      INTEGER IOOUT
C     ..
C     .. External Subroutines ..
      EXTERNAL GETSET
C     ..
      CALL GETSET('GET','IOOUT',IOOUT)
      WRITE (IOOUT,FMT=9000) CSUBRT, CMSGNM
 9000 FORMAT (' >',A,' ',A)
      RETURN
      END
