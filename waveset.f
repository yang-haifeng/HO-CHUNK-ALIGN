      subroutine waveset(cwav,cfldust,ctherm,nlam,MXLAM)
c
c brilliant method for handling wavelength loop arrays
c
c file created: 1999/11/20
c
c history:
c
c
      implicit none

	integer MXLAM,nlam,i
	character cwav(MXLAM)*(*),cfldust(MXLAM)*(*),ctherm(MXLAM)*(*)

	open(unit=55,file='wave.in',status='old')

	i=1
010   read(55,*,end=100) cwav(i),cfldust(i),ctherm(i)
      i = i + 1
	goto 010
100   nlam = i-1

      close(55)

      return
	end