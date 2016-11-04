      SUBROUTINE NAMER(IINC,CDESCR,CFLNAM,MXNAM)
C***********************************************************************
C Purpose: to generate file names for specific (wavelength,
c           inclination angle.
C Given:
C         IINC=integer (0-99) identifying size
c         CDESCR = wavelength description string
c         MXNAM = integer dimension of CFLNAM
C Returns:
C         CFLNAM = 7 file names for output
c file created 1999/11/20
c
c modification history:
c
c
C***********************************************************************
      implicit none

C Arguments:
      INTEGER IINC,MXNAM
      CHARACTER CFLNAM(MXNAM)*(*),CDESCR*(*)
C
C Local variables:
      INTEGER idx,IINC0
      CHARACTER N(10)*1,V1(2)*1

	data N/'0','1','2','3','4','5','6','7','8','9'/
C***********************************************************************

C     F_00_ww.dat    for images F,I,Q,U,V

c      finites.out  plflux.out
c      fluxap_ww.out  pfluxap_ww.out  where ww is wavelength

C
C    Set IINC0=IINC-1 to run over range 0-99
      IINC0=IINC-1
      V1(1)=N(IINC0/10+1)
      V1(2)=N(IINC0-10*(IINC0/10)+1)
      
c     if there are any spaces, idx is last non-space character 
	idx = INDEX(CDESCR,' ')-1
      if (idx.eq.0) idx=len(CDESCR)

      cflnam(1) = 'F_'//V1(1)//V1(2)//'_'//CDESCR(1:idx)//'.dat'
      cflnam(2) = 'I_'//V1(1)//V1(2)//'_'//CDESCR(1:idx)//'.dat'
      cflnam(3) = 'Q_'//V1(1)//V1(2)//'_'//CDESCR(1:idx)//'.dat'
      cflnam(4) = 'U_'//V1(1)//V1(2)//'_'//CDESCR(1:idx)//'.dat'
      cflnam(5) = 'V_'//V1(1)//V1(2)//'_'//CDESCR(1:idx)//'.dat'
      cflnam(6) = 'finites_'//CDESCR(1:idx)//'.dat'
      cflnam(7) = 'pflux_'//CDESCR(1:idx)//'.dat'
      cflnam(8) = 'star_'//CDESCR(1:idx)//'.dat'
	cflnam(9) = 'finites_err'//CDESCR(1:idx)//'.dat'
	cflnam(10) = 'vger_'//CDESCR(1:idx)//'.dat'
	cflnam(11) = 'imivger_'//CDESCR(1:idx)//'.dat'
	cflnam(12) = 'imqvger_'//CDESCR(1:idx)//'.dat'
	cflnam(13) = 'imuvger_'//CDESCR(1:idx)//'.dat'
	cflnam(14) = 'imvvger_'//CDESCR(1:idx)//'.dat'

      RETURN
C
      END
