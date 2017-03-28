      subroutine disk

c     monte carlo radiative transfer for disk surrounding star.
c     uses cartesian coordinates, arbitrary density distribution,
c     arbitrary disk structure.
c     This program works for extremely (geometrically) thin disks of large
c     radial extent plus tenuous envelope.  In order to do the
c     radiative transfer in the thin disk, the opacity is calculated
c     at each step, and the photon path integration is variable.
c     In the disk, the step size is small for directions perpendicular
c     to the disk, say zmindisk/10.  Outside the
c     disk, z gt. zmaxdisk, the step size is much larger, zmax/100
c     or zmax/200, something on that order.
c     Note that zmindisk is probably about .1 stellar radius.
c     zmax is about 10000 stellar radii, that is, about 100 AU.
c     The envelope extend is even larger for protostars--10^4 AU,
c     so the step size increases with distance from the source.

c     calls these subroutines:
c        stokes
c        initp
c        opacdisk
c        opacinfall
c        dels
c
c history:
c 99/04/12 (mjw): add cputime for photon loop update output
c                 (this breaks the overall CPU clock in newtts)
c 99/11/20 (baw):  big changes...
c
c     *************************************************************

      implicit none

      include 'stokes.txt'
      include 'tts.txt'
      include 'opacin.txt'
      include 'taunum.txt'
      include 'random.txt'
      include 'out.txt'
      include 'spot.txt'
      include 'vger.txt'
      include 'grid.txt'

      real*8 sini(nmu),xran,t,nu,bnu,zmaxi,xmax,angle,dsl,rnri,rnzi
     1   ,fractx,area,rstarcgs,dconst,rin,rout,dr,thetb
     1   ,zout,zin,dz,rad,normstar1,normstar2,tsum,coslp,sinlp,kaps
	1   ,cosnorm,sipnew,fact,rtmp,taunoam,eps,sqpnew,supnew,svpnew

      integer ii(3),ns,nd,ne,icount,ithet,gridflag,iphot,i,nstart,ir
     1   ,nend,iz,outflag,it,ip
      logical ifound
      real ran2

c cpu time variables
      real etime,cpusec,time(2)
      character cmsgnm*70
      external etime
     
      icount=0

c     set up density grid
      call gridset

c     The flux is summed to image(ix,iy,it)
c     polarization to imagei(ixp,iyp,it),imageq...,imageu...

c     inclination arrays.
      do ithet=1,nmu
         sini(ithet)=dsqrt(1.d0-u(ithet)**2)
      end do

      if (ispot.eq.1.and.spotflag.eq.1) call spotset(npsav,wave,tstar)
c     only call spotset once, it will set spotflag=0 after first call

      nscat=0.d0
      tot=0.d0

c     x,y arrays (images)
      fractx=2.d0*rmaxi/float(nx)
c     x,y arrays for peeled-off image
      fractxh=2.d0*rmaxi/float(nxhst)

c     eps for stepping through grid (to step just beyond cell)
      eps=1.d-9

      call flush(6)
      
c     np=npsav
      print*,'np',np
      flux=0.d0
      sflux=0.d0
      aflux=0.d0

      if (itherm.eq.0) then
         ns=np+npout
         nd=0
         ne=0
      else
c     calculate luminosity of star
         nu=2.99792458d14/wave  !in microns
         call plancknu(tstar,nu,bnu)
         print*,'stellar Bnu ',bnu
c     bnu=bnu*4*pi     !convert from energy/s/cm2/Sr to energy/s/cm2
         rstarcgs=rstar*rsol
         area=pi*(rstarcgs)**2
         normstar=area*bnu      !this is an energy/s; normalization 
                                !for output
         print*,'stellar energy at this wavelength at earth',normstar
         print*,'luminosity ',4.d0*area/3.826d33*tstar**4*5.669d-5 
                                !energy/s
         normstar1=normstar*4.d0*pi !flux=pi*bnu area=4*pi*rstar**2
         print*,'luminosity at this wavelength ',normstar1
         ns=np+npout
      endif

c     ********  do loop over stellar photons  *******
      print*,'ns',ns
      print*,'npout',npout

      do iphot=1,ns
         
c     if (iphot.gt.30000.and.iphot.lt.31000) then
c     continue
c     print*,iphot
c     endif
         if(mod((iphot),iwrite).eq.0) then
c     cpu time
            cpusec=etime(time)
c     use explicit format to keep g77 happy
            write(cmsgnm,'(i12,a20,f11.2,a)') iphot,
     $           ' photons completed. ',cpusec,' = CPU time (sec)'
            call WRIMSG('MAIN',cmsgnm)
         end if
         
c     inside illumination
         if (iphot.ge.npout) then
c     initialize stokes parameters, incident angles
            if (ispot.eq.0) then
               call initp
            else
               call initpspot(iphot)
            endif
c     print*,'cost,sint',cost,sint
            coslp=dcos(lp)
            sinlp=dsin(lp)
c     zstar
            zp=cosb*rmin*(1.d0+eps)
c     rstar=xp if all photons start off on x-axis
            rp=sinb*rmin*(1.d0+eps)
            xp=rp*coslp
            yp=rp*sinlp
            rsq=zp**2+rp**2
            rtot=dsqrt(rsq)
            ux=sint*cosp
            uy=sint*sinp
            uz=cost
            iflag=0
            ii(1)=1             !index of radial grid (stellar 
                                !surface is 1)
            if (ntg.gt.1) then
               thetb=acos(cosb)
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
            endif

C            print*,'xp,yp,zp,ux,uy,uz',xp,yp,zp,ux,uy,uz
C            print*,'ir,it,ip',ii(1),ii(2),ii(3)
C            print*,'cost,sint,cosp,sinp',cost,sint,cosp,sinp
C            stop
c     first, peel off direct flux
c     weight photon intensity by angle between normal and 
c     photon direction, since emitted from a surface.
            if (ipeel.eq.1) then
               cosnorm=cosb*coste+(sinb*sinte*(cospe*coslp+
     1              sinpe*sinlp))
               if (limb.eq.0) then
c     intensity constant, energy/Sr proportional to mu
                  sipnew=4.*sip*cosnorm !normalization = 2
                  sqpnew=sqp*sipnew
                  supnew=sup*sipnew
                  svpnew=svp*sipnew
               else
c     intensity goes as (1+mu).  energy/Sr has another factor
c     of mu.
                  sipnew=12.d0/5.d0*sip*(cosnorm+cosnorm*cosnorm)
                  sqpnew=sqp*sipnew
                  supnew=sup*sipnew
                  svpnew=svp*sipnew
                  kaps=1.d0
               endif
c     sipnew=sip	 
               if (cosnorm.gt.0.) then
                   call peeloff(xp,yp,zp,sipnew,sqpnew,supnew,svpnew
     1           ,cost,sint,cosp,sinp,phi
     1           ,pi,r2p,hit,htot,rsq,rtot
     1           ,peak,tsum,ii,idust,iflag,iphot,eps,kaps)   
               endif
               if (cosnorm.gt.0.and.iveeg.eq.1) then
                  call vger(xp,yp,zp,sipnew,sqpnew,supnew,svpnew
     1           ,cost,sint,cosp,sinp,phi
     1           ,pi,r2p,hit,htot,rsq,rtot
     1           ,peak,tsum,ii,idust,iflag,iphot,eps)   
               endif
            endif
         else
c     outside illumination
            call initpout
            ux=sint*cosp
            uy=sint*sinp
            uz=cost
            iflag=0
c     don't peel off direct flux because star is outside image
c     field  . 
c     but send star flux to vger   
            if (iveeg.eq.1) then
               call vger(xse,yse,zse,sip,sqp,sup,svp
     1              ,cost,sint,cosp,sinp,phi
     1              ,pi,r2p,hit,htot,rsq,rtot
     1              ,peak,tsum,ii,idust,iflag,iphot,eps)   
               if (iphot.eq.1) then
                  print*, 'tau from outside star to vger ',
     1                 tsum
                  write(12,*) 'tau from out star to vger ',
     1                 tsum
               endif
            endif
c     now calculate position in envelope of photon of selected
c     direction. note that it may not hit envelope at all.
            call Rdist(ifound,t,ux,uy,cost,xs,ys,zs,rmax)
            if (ifound.eqv..false.) go to 5
            fact=1.0d0+1.d-3
            xp=xs+ux*t*fact
            yp=ys+uy*t*fact
            zp=zs+cost*t*fact
            rsq=xp**2+yp**2+zp**2
            rtot=dsqrt(rsq)
c            rp=dsqrt(xp**2+yp**2)
c            rtmp=dsqrt(xp**2+yp**2+zp**2)	 
         endif
         
c     now calculate scattered flux
c     print*,'before propagate, ir,it,ip',ii(1),ii(2),ii(3)
         call propagate(iphot,sini,rmaxi,fractx,ii,eps)
         
 5    end do
      
c     ***** end of loop over stellar photons  ******
      
      
      write(6,*) 'sini,cosi'
      write(6,*) (sini(i),i=1,nmu),(u(i),i=1,nmu)
      write(6,*) 'fractx,xmax',fractx,rmaxi
      
      write(6,*) 'fraction of phots scattered outside image'
      write(6,*) float(icount)/float(np)
      
      return
      end

c     *********************************************************
