c     ********************************************************

      subroutine initpout

c     called by subroutine disk
c     initialize variables for new photon which scatters in disk
c
c history:
c 00/03/19 (mjw):  set initial stokes vector using Stokes vector
c                      from REAPAR if iveeg.eq.1 (CVEEJ='YES')
c

      implicit none

      real*8 xran,costi,sinti,cospn,phinew
      real ran2

      include 'tts.txt'
      include 'stokes.txt'
      include 'random.txt'
      include 'vger.txt'
      
c     incident radiation chosen to be unpolarized.
c     can put in a desired polarization in the limb darkening
c     subroutine.
      if (iveeg.eq.1) then
         sip=isrc0
         sqp=qsrc0
         sup=usrc0
         svp=vsrc0
      else
         sip=1.d0
         sqp=0.d0
         sup=0.0d0
         svp=0.d0
      endif
      
c     sample uniform angle
      xran=ran2(i1)
      cost=1.d0-2.d0*xran
      sint=dsqrt(1.d0-cost**2)
      xran=ran2(i1)
      phi=r2p*xran
      cosp=dcos(phi)
      sinp=dsin(phi)

      return
      end
      

c     **********************************************************
