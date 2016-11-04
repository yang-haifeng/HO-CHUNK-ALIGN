c     ********************************************************

      subroutine initp

c     called by subroutine disk
c     initialize variables for new photon which scatters in disk

      implicit none

      real*8 xran,costi,sinti,cospn,phinew
      real ran2

      include 'tts.txt'
      include 'stokes.txt'
      include 'random.txt'
      
c     incident radiation chosen to be unpolarized.
c     can put in a desired polarization in the limb darkening
c     subroutine.
      sip=1.d0
      sqp=0.d0
      sup=0.d0
      svp=0.d0

c     sample cosb,sinb, position on star in latitude.
      xran=ran2(i1)
      cosb=1.d0-2.d0*xran
      sinb=dsqrt(1.d0-cosb**2)
      
c     sample angle of exit using appropriate limb darkening law.
      if(limb.eq.0) then
c        isotropic intensity, sample from n=mu
         xran=ran2(i1)
         cost=dsqrt(xran)
      else
         xran=ran2(i1)
         call darkening(xran,cost)
      end if
      sint=dsqrt(1.d0-cost*cost)

c     sample phi
      xran=ran2(i1)
      phi=r2p*xran
      cosp=dcos(phi)
      sinp=dsin(phi)

c     sample azimuthal coord (longitude)
      xran=ran2(i1)
      lp=r2p*xran

c     transform to coordinate system of star
      if(dabs(cosb).lt.0.9999999d0) then
         costi=cost*cosb-sint*sinb*cosp
         if(dabs(costi).gt.1.d0) then
            write(6,*) 'initp: costi gt 1',costi
            if(costi.gt.1.d0) costi=1.d0
            if(cost.lt.-1.d0) costi=-1.d0 
         end if
         sinti=dsqrt(1.d0-costi**2)
         if(sinti.gt.0.0000001d0)then
            cospn=(cost-costi*cosb)/(sinti*sinb)
            if(dabs(cospn).gt.1.d0) then
               if (dabs(cospn).gt.1.01d0) 
     1            write(6,*)'cospn gt 1 in initp',cospn
               if (cospn.lt.0.d0) then
                  cospn=-1.d0
               else
                  cospn=1.d0
               end if
            end if
            if(phi.lt.pi) then
               phinew=dacos(cospn)
               if(phinew.gt.pi) then
                  write(6,*) 'phinew wrong in initp'
                  write(6,*) 'phinew,pi',phinew,pi
                  phinew=pi
               end if
            else
               phinew=r2p-dacos(cospn)
               if(phinew.lt.pi) then
                  write(6,*) 'phinew wrong in initp'
                  write(6,*) 'phinew,pi',phinew,pi
                  phinew=pi
               end if
            end if
         else
            phinew=0.d0
         end if
         phi=phinew
         cost=costi
         sint=sinti
      else
         if(cosb.lt.0.d0) cost=-cost
      end if   
      phi=phi+lp
      if(phi.gt.r2p) phi=phi-r2p
      if(phi.lt.0.d0) write(6,*)'phi lt 0',phi
      if(phi.gt.r2p) write(6,*) 'phi gt 2pi',phi
      sinp=dsin(phi)
      cosp=dcos(phi)

      return
      end
      

c     **********************************************************