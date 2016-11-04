      subroutine spotset(np,wave,tstar)
      implicit none

      include 'spot.txt'
      include 'stokes.txt'

      real*8 nu,wave,bnustar,bnuspot,f,tstar
      integer np

      thspot=thspot*pi/180.d0
      thspotin=thspotin*pi/180.d0
      th0sp1=(90.d0-spotlat)*pi/180.d0
      phi0sp1=spotlon
      th0sp2=pi-th0sp1
      phi0sp2=phi0sp1+pi

      scspot1=dsin(th0sp1)*dcos(phi0sp1)
      ssspot1=dsin(th0sp1)*dsin(phi0sp1)
      cspot1=dcos(th0sp1)
      scspot2=dsin(th0sp2)*dcos(phi0sp2)
      ssspot2=dsin(th0sp2)*dsin(phi0sp2)
      cspot2=dcos(th0sp2)

      csprad=dcos(thspot)
      cspradin=dcos(thspotin)

      spotflag=0

c     determine ratio of spot photons to non-spot based on size and
c     temperatures
      nu=2.9979d14/wave    !in microns
      print*,'nu',nu
      call plancknu(tstar,nu,bnustar)
      call plancknu(tspot,nu,bnuspot)

c     one spot
      f=(1.d0-csprad)/(1.d0+csprad)*bnuspot/bnustar
c     two spots
      if (nspot.eq.2) f=2.d0*f
c     haven't done formula for ring....

      npspot=float(np)*f/(1.d0+f)

      print*,'ratio of Bspot/Bstar ',bnuspot/bnustar
      print*,'ratio of spot photons to star photons ',f
      print*,'spot photons, star photons, total ',npspot, np-npspot, np
      print*, ' '

      return
      end
