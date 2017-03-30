c     **************************************************************** 
      
      subroutine initpDust(ii)

c     called by subroutine disk
c     initialize variables for new photon based on the thermal emission 
c     of the dust

      implicit none

      real*8 xran,costi,sinti,cospn,phinew
      real ran2
      integer ii(3)

      include 'tts.txt'
      include 'stokes.txt'
      include 'random.txt'
      include 'grid.txt'
      include 'taunum.txt'

c     incident radiation chosen to be unpolarized. For now.
      sip=1.d0
      sqp=1.d0
      sup=1.d0
      svp=1.d0

c     sample radius of the location
      xran=ran2(i1)
      xran*nrg
      cosb=
      sinb=

      zp=
      xp=
      yp=

      ux=
      uy=
      uz=

      ii(1)=
      ii(2)=
      ii(3)=

      end
