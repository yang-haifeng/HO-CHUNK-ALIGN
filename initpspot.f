c     ********************************************************

      subroutine initpspot(iphot)

c     called by subroutine disk
c     initialize variables for new photon which scatters in disk

      implicit none

      real*8 xran,costi,sinti,cospn,phinew
	real*8 sthet,cthet,sphi,cphi,probspot,dotsp1,dotsp2
      real ran2
      integer iphot

      include 'tts.txt'
      include 'stokes.txt'
      include 'random.txt'
	include 'spot.txt'
      
c      print*,'iphot, nspot, npspot ',iphot,nspot,npspot

c     incident radiation chosen to be unpolarized.
c     can put in a desired polarization in the limb darkening
c     subroutine.
      sip=1.d0
      sqp=0.d0
      sup=0.d0
      svp=0.d0

c     sample cosb,sinb, position on star in latitude.
100   continue
      xran=ran2(i1)
      cosb=1.d0-2.d0*xran
      sinb=dsqrt(1.d0-cosb**2)

c     sample azimuthal coord (longitude)
      xran=ran2(i1)
      lp=r2p*xran
	cphi=cos(lp)
	sphi=sin(lp)
	sthet=sinb
	cthet=cosb

c        **************SPOT STUFF***************************
          dotsp1=scspot1*sthet*cphi+
     $            ssspot1*sthet*sphi+
     $            cspot1*cthet
          dotsp2=scspot2*sthet*cphi+
     $            ssspot2*sthet*sphi+
     $            cspot2*cthet
c          if((dotsp1.lt.csprad).and.    ! for two spots
c     $       (dotsp2.lt.csprad))then    ! that are filled in

c          if( ((dotsp1.lt.csprad).or.(dotsp1.gt.cspradin)) .and. 
c     $         ((dotsp2.lt.csprad).or.(dotsp2.gt.cspradin)) )
c     $        then

c          if(dotsp1.lt.csprad)then  ! for one spot filled in
c             probspot=ampl
c          else
c             probspot=1.d0
c          endif
c          xran=ran2(i1)
c          if(xran.gt.probspot) then
c            goto 100
c          endif

c     baw, 6/27/00, doing it differently;
c     use equation 3 from Wood et al 2000.  If iphot < N_spot, put photon
c     source at spot.  otherwise, put it outside spot.

      if (nspot.eq.1) then
c      print*,'dotsp1,csprad ',dotsp1,csprad
c     this is for a single spot.
         if (iphot.lt.npspot) then
c        note , csprad goes from 1 to 0 as theta goes from 0 to 90
            if (dotsp1.lt.csprad) go to 100
         else
            if (dotsp1.gt.csprad) go to 100
         endif
      else
c     this is for 2 spots
         if (iphot.lt.npspot) then
            if (dotsp1.lt.csprad.and.dotsp2.lt.csprad) go to 100
         else
            if (dotsp1.gt.csprad.or.dotsp2.gt.csprad) go to 100
         endif
      endif

c        ***************************************************

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
