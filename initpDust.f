c     **************************************************************** 
      
      subroutine initpDust(ii)

c     called by subroutine disk
c     initialize variables for new photon based on the thermal emission 
c     of the dust

      implicit none

      real*8 xran,thetb,currentEps
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

c     sample cost, sint, direction of the photon in theta
      xran=ran2(i1)
      cost=1.d0-2.d0*xran
      sint=dsqrt(1.d0-cost**2)
c     sample cosp, sinp, direction of photon in phi
      xran = ran2(i1)
      phi=r2p*xran
      sinp=dsin(phi)
      cosp=dcos(phi)

10    xran = ran2(i1)
      ir = int(xran*(nrg-1))+1
      xran = ran2(i1)
      it = int(xran*(ntg-1))+1
      xran = ran2(i1)
      ip = int(xran*(npg-1))+1

      currentEps=massarr(ir,it,ip)*temparr(ir,it,ip)
      xran = ran2(i1)
      if (xran*maxEps>currentEps) then
        goto 10
      endif

c     sample cosb, sinb, location of the photon in latitude
      thetb = 0.5*(thetarr(it)+thetarr(it+1))
      cosb=dcos(thetb)
      sinb=dsin(thetb)
c     sample location of the photon in longitude
      lp=0.5*(phiarr(ip)+phiarr(ip))

      rtot = 0.5*(rarr(ir)+rarr(ir+1))
      rsq = rtot*rtot

      zp=rtot*cosb
      xp=rtot*sinb*dcos(lp)
      yp=rtot*sinb*dsin(lp)

      ux=sint*cosp
      uy=sint*sinp
      uz=cost

      end
