c     *******************************************************
      
      subroutine setup_wave(iwav,cfldust,ctherm,MXWAV)
c     
c     history:
c 00/03/19 (mjw):  set vger position and outside source using
c                      values from parameter file (e.g. REAPAR)
c 00/07/27 (mjw):  fix formatting issue for lf95/gcc (a**-1.5 type)
c 
c     
      
      implicit none
      
      include 'tab.txt'
      include 'tts.txt'
      include 'stokes.txt'
      include 'tabl.txt'
      include 'opacin.txt'
      include 'out.txt'
c     include 'mrncoef.txt'
      include 'random.txt'
      include 'vger.txt'
      include 'taunum.txt'
      
      integer MXWAV
      character cfldust(MXWAV)*(*),ctherm(MXWAV)*(*)
      
      real*8 a1,c,c1
     1     ,z1cgs,rmincgs,rmindcgs,rmaxdcgs,rddustcgs,p
     1     ,dnoam,hnoam,rnoam,taunoam,c0r2h,xnoam,deltax,deltay,deltaz
     $     ,stot,stot2,roota,rootb,rootc,rootd,t,rvg,fact,rtmp
     1     ,rho1,sigz1,rdxcgs,kappa

      
      real*8 check,linterp
      
      integer inoam,iwav,i,iinc
      logical ifound
      character*20 routine
      character filename*80
      
      write(12,*) ' '
      write(12,*) 'wavelength index ',iwav
      write(12,*) 'include thermal emission? ',CTHERM(iwav)
c     call reatab(cfldust(iwav),rlam,kapd,wave)
      filename=cfldust(iwav)
      call read_16scat(wave,filename)
      
      if (ctherm(iwav).eq.'yestherm'.or.ctherm(iwav).eq.'YESTHERM') 
     $     then
         itherm=1
      else if (ctherm(iwav).eq.'notherm'.or.ctherm(iwav).eq.'NOTHERM') 
     1        then
         itherm=0
      else
         print*,'ERROR, ctherm(iwav) not properly defined'
         stop
      endif
      
c     first calculate rmine based on dust sublimation temperature
      if (itherm.eq.1) then
         rmine=(1500.d0/tstar)**(-2.5d0)
         if (rmine.lt.1.) rmine=1.
         print*,'dust destruction radius in rstar ',rmine
         rddust=(1500.d0/tstar)**(-4.d0/3.d0)
         if (rddust.lt.1.d0) rddust=1.d0
         a=0.75
         b=0.
         zmin=.01
         print*,'since itherm=1, setting rmine,rddust,a,b,zmin= '
         print*,rmine,rddust,a,b,zmin
         if(rmaxd.gt.0) then
            zmax=z1*rmaxd**b
         else
            zmax=0.d0
         end if
      else
c     use input values
      endif
      
c     for calculating density in cgs
      z1cgs=rstar*rsol*z1
      rmincgs=rstar*rsol
      rmindcgs=rmind*rmincgs
      rmaxdcgs=rmaxd*rmincgs
      rddustcgs=rddust*rmincgs
c     calculate density from massd=mass of disk. 
c     = 2*pi*(integral(density*r*dr*dz)).	
      a1=1.d0-a  
      p=b-a
      
c     scaling factor for computing scattered light in peeled-off image
      sfactor=1.d0/(4.d0*pi)
c     vgfactor=(rstar*rsol)**2
      vgfactor=rstar**2
      
      if(rmaxd.gt.0.d0) then
         c=r2p**1.5d0*rmincgs**2*z1cgs
         if(p.ne.-2) then
            rho0=massd*msol*(p+2)/c/(rmaxd**(p+2)-rmind**(p+2))
         else
            rho0=massd*msol/c/dlog(rmaxd/rmind)
         end if
         write(6,*)'rho0 of disk',rho0
c     opacity in units of rstar**-1 is kappa*rho*rmin.  c0 is multiplier
c     used in radiative transfer.
         c0=rho0*rmincgs
         c1=rho0*rmincgs**a !for optical depth calcs here
c        find direction-dependent opacity
         call locate(cosbarr,MXTHETB,nthetab,1.d-10,iinc)
         kappa=linterp(cosbarr,kapd,MXTHETB,1.d-10,iinc)
c     taux is region of disk we will finely sample in grid
c     set taux=20 for now,  (not used when massd=0)
         taux=20.d0
         rdx=rmind
         if(a.ne.1) then
            taur=c1*kappa/a1*(rmaxdcgs**a1-rddustcgs**a1)
            if (massd.gt.0.d0) then
               rdxcgs=(taux*a1/(c1*kappa)+rddustcgs**a1)**(1./a1)
               rdx=rdxcgs/rmincgs
            endif
         else
            taur=c1*kappa*dlog(rmaxd/rddust)
            if (massd.gt.0.d0) then
               rdx=rddust*exp(taux/(c1*kappa))
            endif
         end if
         if (rdx.gt.rmaxd) rdx=rmaxd
      else
         c=0.d0
         c0=0.d0
         rho0=0.d0
         taur=0.d0
      end if
      rhod0=rho0
      write(6,*)'c0,c1,taur',c0,c1,taur
      if(rmaxd.gt.0.) then
         call locate(cosbarr,MXTHETB,nthetab,0.999999d0,iinc)
         kappa=linterp(cosbarr,kapd,MXTHETB,0.999999d0,iinc)
         tauzmin=kappa*dsqrt(r2p)*z1cgs*rho0*rddust**p
         tauzmax=kappa*dsqrt(r2p)*z1cgs*rho0*rmaxd**p
         tauz30=kappa*dsqrt(r2p)*z1cgs*rho0*(30*autors)**p
         write(6,*) 'tauzmin,tauzmax,tauz30',tauzmin,tauzmax,tauz30
c     density at rmind,30AU, rmaxd
         rhomin=rho0*(rmind/rmin)**(-a)
         rho30=rho0*(30.d0*autors/rmin)**(-a)
         rhomax=rho0*(rmaxd/rmin)**(-a)
         write(6,*)'density at rmind,30au,rmaxd',rhomin,rho30,rhomax
c     density at 1 AU, sigma at 1 AU
         rho1=rho0*(1.d0*autors/rmin)**(-a)
         sigz1=dsqrt(r2p)*z1cgs*rho0*(1.d0*autors)**p
         print*,'rho_0, Sigma at 1 AU ',rho1,sigz1
      end if
      
c     make a file of disk height and constants vs radius for use with 
c     error function to calculate taus.
      if(rmaxd.gt.0.d0) then
         open(unit=14,file='tau.dat',status='unknown')
         write(14,*)
     1        'radius (r),scale height (h),c0*r**-a*2h,tau at 3h,F(x1)'
         dnoam=rmaxd/10.d0
         do inoam=1,10
            rnoam=(dfloat(inoam)-0.5d0)*dnoam
            hnoam=z1*(rnoam/rmin)**b
            c0r2h=kappa*c0*(rnoam/rmin)**(-a)*2.*hnoam
            if(c0.gt.0.) then
               taunoam=c0r2h*0.0013d0
               xnoam=1.d0-1.d0/c0r2h
            end if
            write(14,*) rnoam,hnoam,c0r2h,taunoam,xnoam
         end do
         close(14)
      end if
      
c     find tau for H=3.5 at rmin
      hnoam=z1*rmind**b
      c0r2h=kappa*c0*(rmind)**(-a)*2.*hnoam
      if(c0.gt.0.) then
         taunoam=c0r2h*0.0002d0
      end if
      write(12,*) 'optical depth to H=3.5 at rmind',taunoam
      print*, 'optical depth to H=3.5 at rmind',taunoam
      
      outillum=1                !for now
c     set vger position in disk
c     if (iveeg.eq.1) then
      if (outillum.eq.0) then
         xvg=1.*autors
         yvg=0.
         hnoam=z1*xvg**b
         zvg=3.44*hnoam
         c0r2h=c0*kappa*(xvg)**(-a)*2.*hnoam
         if(c0.gt.0.) then
            taunoam=c0r2h*0.0003d0
            xnoam=1.d0-1.d0/c0r2h
         end if
         write(12,*) 'vger position, x,y,z ',xvg,yvg,zvg
         write(12,*) 'optical depth from top of disk to vger ',taunoam
      else
c     set position of star emitting outside photons
c     xs=rmax*1.4d0
c     ys=0.0d0*rmax
c     zs=rmax*1.20
         xs=rmax*xsrc0
         ys=rmax*ysrc0
         zs=rmax*zsrc0
         rs=dsqrt(xs**2+ys**2+zs**2)
c     set vger position in outer envelope (for outside illumination)
c     xvg=rmax*0.4d0
c     yvg=0.d0
c     zvg=-rmax*0.2d0
         xvg=rmax*xvg0
         yvg=rmax*yvg0
         zvg=rmax*zvg0
         rvg=dsqrt(xvg**2+yvg**2+zvg**2)
         
c     calculate direction cosines from star to vger
         deltax=xs-xvg
         deltay=ys-yvg
         deltaz=zs-zvg
         stot2=deltax**2+deltay**2+deltaz**2
         stot=dsqrt(stot2)
         cost=deltaz/stot
         cost=check(cost,routine)
         sint=dsqrt(1.-cost**2) !theta ranges from 0-pi so sint>0
         phi=datan2(deltay,deltax)
         if (phi.lt.0) phi=phi+r2p
         cosp=dcos(phi)
         sinp=dsin(phi)      
         ux=sint*cosp
         uy=sint*sinp
c     calculate x,y,z position of direct light from star to 
c     vger once it hits the envelope
         call Rdist(ifound,t,ux,uy,cost,xvg,yvg,zvg,rmax)
         fact=0.9999d0
 333     xse=xvg+ux*t*fact
         yse=yvg+uy*t*fact
         zse=zvg+cost*t*fact
         rpse=dsqrt(xse**2+yse**2)
         rtmp=dsqrt(xse**2+yse**2+zse**2)	 
         if (rtmp.ge.rmax) then
            fact=fact*0.99d0
            go to 333
         endif  
      endif
c     endif
      
      
      if (iwav.eq.1) then
         write(12,*) 'nri,nzi,nfinei',nri,nzi,nfinei
         write(12,*) 'nfined',nfined
         write(12,*) 'nx',nx
         write(12,*) 'infall rate (solar m per yr),ihole',rate,ihole
         write(12,*) 'streamline angle for wind',thetmu0
         write(12,*) 'idust',idust
         write(12,*) 'limb,occult',limb,occult
         write(12,*) 'z=z1*r**b. b = ',b
         write(12,*) 'z1 ',z1
         write(12,*) 'rho(z=0)=c0*r**(-a). a = ',a
         write(12,*) 'c0 ',c0
         if (itherm.eq.1) then
            write(12,*) 'includes thermal emission'
            write(12,*) 'tstar = ',tstar
         else
            write(12,*) 'does not include thermal emission'
         endif
         write(12,*) 'density at rmind,30AU,rmaxd',rhomin,rho30,rhomax
         write(12,*) 'oa (tan(oa)=zmax/rmaxd)',oadeg
         write(12,*) 'mass in disk, mass of core ',massd,massc
         write(12,*) 'density at midplane of disk, rmin',rho0
      endif
      write(12,*) 'tauzmin,tauz30,tauzmax ',tauzmin,tauz30,tauzmax
      
      
      return
      end
      
