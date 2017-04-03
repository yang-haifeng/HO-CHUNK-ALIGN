c     **************************************************************** 
      
      subroutine initpDust(ii)

c     called by subroutine disk
c     initialize variables for new photon based on the thermal emission 
c     of the dust

      implicit none

      real*8 xran,thetb
      real ran2
      integer ii(3),ir,it,ip

      include 'tts.txt'
      include 'stokes.txt'
      include 'random.txt'
      include 'grid.txt'
      include 'taunum.txt'

c     incident radiation chosen to be unpolarized. For now.
      sip=1.d0
      sqp=0.d0
      sup=0.d0
      svp=0.d0

c     sample cosb, sinb, location of the photon in latitude
      xran=ran2(i1)
      cosb=1.d0-2.d0*xran
      sinb=dsqrt(1.d0-cosb**2)
c     sample location of the photon in longitude
      xran = ran2(i1)
      lp=r2p*xran

c     sample cost, sint, direction of the photon in theta
      xran=ran2(i1)
      cost=1.d0-2.d0*xran
      sint=dsqrt(1.d0-cost**2)
c     sample cosp, sinp, direction of photon in phi
      xran = ran2(i1)
      phi=r2p*xran
      sinp=dsin(phi)
      cosp=dcos(phi)

c     sample radius of the location
      xran = ran2(i1)
      rtot = (xran**(1.d0/3.d0))*rarr(nrg) ! Uniformly pick a location 
c     rtot = xran*rarr(nrg) ! Uniformly pick a location 
      rsq = rtot*rtot

      zp=rtot*cosb
      xp=rtot*sinb*dcos(lp)
      yp=rtot*sinb*dsin(lp)

      ux=sint*cosp
      uy=sint*sinp
      uz=cost

      call locate(rarr,nrg,nrg,rtot,ir)
      ii(1)=ir
      if (ntg.gt.1) then
         thetb=dacos(cosb)
         call locate(thetarr,ntg,ntg,thetb,it)
         ii(2)=it
      else
         ii(2)=1
         it=1
      endif
      if (npg.gt.1) then
         call locate(phiarr,npg,npg,lp,ip)
         ii(3)=ip
      else
         ii(3)=1
         ip=1
      endif

      end
