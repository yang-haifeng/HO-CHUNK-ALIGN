      real*8 function linterp(xtab,ztab,nx,xval,ix)
c
c  linearly interpolate a 1-D table.  assumes you already called a 'locate'
c  subroutine to find position (ix) in array.
c  
c  b.whitney, space science institute
c  apr 7, 2000

c  history
c  00/12/12  (baw) specific to mcvger16, ztab array is real*4.

      implicit none

      integer ix,nx
      real*8 xtab(nx),xx,xval
      real*4 ztab(nx)

      xx=(xval-xtab(ix))/(xtab(ix+1)-xtab(ix))

      linterp=ztab(ix)+xx*(ztab(ix+1)-ztab(ix))

      return
      end
