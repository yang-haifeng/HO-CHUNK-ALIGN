      subroutine gridset

c  set up grid:  3-D spherical grid with variable spacing in r and 
c  theta (set by exponents, rexp, texp)
c  make a bunch of arrays necessary for find_wall subroutine

c history:
c 00/09/06 (baw): write subroutine
	
      implicit none
      include 'tts.txt'
      include 'out.txt'
      include 'grid.txt'
      include 'opacin.txt'

      real*8 dr,dt,dp,rad,phi,thet,pi,r2p,densd,dense,vol
     $ ,cost,sint,pihalf,eps,tau,rming1,rming2,rming3,tiny
     $     ,rmincgs,taud,taue,mu0,tauv,taufmin,mu,thetatm
     $     ,thetmin,tauzrat,sigmu,rsonr,Tdisk,taudisk
     $     ,mudisk,c1,a1,tauw,rtop,rfac,tauRr,taumid,tauRth
     $     ,tauRos,res,thettop
      real xyarr(200,200),x,y,xmax,dx,r
      real chiR,ierfc
      integer ir,it,ip,nttmp,gflag,n,ix,iy

      real*8 linterp

      pi=4.d0*datan(1.d0)
      r2p=2.d0*pi
      pihalf=pi/2.d0
      tiny=1.d-15
      rmincgs=rstar*rsol

c     convert opacity to units of cm^2/gm*rstar(cgs) because distance
c     is in units if 1/rstar, and dtau=kapd*rho*ds
c     no, kappa is a variable now, so make rho=rho*rmincgs
c     and leave kappa in cm^2/gm units.
c     kapd=kapd*rmincgs

c     calculate some constants for the TSC envelope density
      call envset

c     make grid include minimum and maximum values.

      print*, ' '
      print*, 'grid setup'

c     rgrid
      gflag=0
      rarr(1)=rmin
c     since there's empty space below rmine and rmind....
      if (rddust.lt.rmine) then
         gflag=1
         rarr(2)=rddust
         rming1=rddust
         rming2=rmine
         irbeg=2
      else if (rddust.gt.rmine) then
         gflag=2
         rarr(2)=rmine
         rming1=rmine
         rming2=rddust
         irbeg=3
      else if (rddust.eq.rmine) then
         gflag=3
         rarr(2)=rmine
         rming1=rmine
         irbeg=2
      endif

c     since there's empty space below rmine and rddust....
      dr=(rmax-rming1)/(dfloat(nrg-2))**rexp
      print*,'rexp,dr,rmin',rexp,dr/autors,rmin/autors
      do ir=3,nrg
         rarr(ir)=rming1+dr*(dfloat(ir-2))**rexp
c     print*,'rarr ',ir,rarr(ir)/au
      enddo

c     set rmine to a grid location
      if (gflag.eq.1) then
         call locate(rarr,nrg,rmine,ir)
         rarr(ir+1)=rmine
      endif

c     r-squared array
      do ir=1,nrg
         r2arr(ir)=rarr(ir)**2
c     print first few points of rarr
         if (ir.lt.10) then
            print*,'rarr/rmin,ir ',rarr(ir)/rmin,ir
         endif
      enddo

      open(unit=15,file='rarr.dat',status='unknown')
      write(15,*) nrg,' = number of grid points in r'
      write(15,*) '      index     r/rstar         r(au)'
      do ir=1,nrg
         write(15,*) ir,sngl(rarr(ir)),sngl(rarr(ir)/autors)
      enddo
      close(15)

      print*,'rarr(nrg),r2arr(nrg)',rarr(nrg),r2arr(nrg)

      if (ntg.gt.1) then
c     set up a tmp array going from 0 - 90 from equ. to pole
c     make theta=0 bin 5 degrees wide (otherwise, too much noise,
c     nothing happens there anyway)
         thettop=5.*pi/180.
         if (massd.eq.0.d0) then

            ntatm=1
            thetatm=0.d0
            tmptharr(1)=0.d0
            nttmp=(ntg+1)/2

         else

c     taufmin=min(0.1,0.1*kappaf*tauv)
            taufmin=min(.001,0.001*kappaf*tauv)
            
            print*,'kappaf*tauv',kappaf*tauv
            
            mu=min(mu0*sqrt(2.*dlog(tauv*kappaf/taufmin)),0.5)
            print*,'tauv,kappaf,taufmin',tauv,kappaf,taufmin

c     ************
c     for disk-only model, let mu=0.5 to sample entire disk height
c     at high-res
c     mu=0.5
c     ***********

            thetatm=pihalf-acos(mu)
c     test
            print*,'thetatm',thetatm*180./pi
            nttmp=(ntg+1)/2

c     *********************
c     change this for disk or envelope runs
c     for envelope with 1 degree resolution in polar region
c     except for first bin (theta=5).
            res=1.
c     for disk with nothing in the polar region
c     res = size of angle bin in poles, approximately 
c     res=10.
c     *********************

            ntatm=nttmp-(pihalf-thetatm-thettop)*90./pihalf/res

            thetmin=min(0.1*asin(mu0),thetatm/ntatm)
            
            tmptharr(1)=0.d0
            do it=2,ntatm
               tmptharr(it)=thetmin *
     1              (thetatm/thetmin)**((it-2.d0)/(ntatm-2.d0))
c     print*,'tmptharr ',tmptharr(it)*180.d0/pi
            enddo
            
         endif
         do it=ntatm+1,nttmp-1
            tmptharr(it)=thetatm+
     1           (it-ntatm)*(pihalf-thetatm-thettop)/(nttmp-ntatm)
            print*,'tmptharr(it)',tmptharr(it)
         end do
         tmptharr(nttmp)=pihalf
         
         thetarr(nttmp)=pihalf
         do it=2,nttmp-1
            thetarr(nttmp+1-it)=pihalf-tmptharr(it)
            thetarr(nttmp-1+it)=pihalf+tmptharr(it)
         enddo
         thetarr(1)=0.d0
         thetarr(ntg)=pi
         
c     this is not the eps used in the rest of the code!  see
c     newdisktherm for that
         eps=1.d-8
         do it=1,ntg
            if (thetarr(it).lt.0) then
               thetarr(it)=0.d0
               costarr(it)=1.d0
               sintarr(it)=0.d0
            else if (thetarr(it).gt.pi) then
               thetarr(it)=pi
               costarr(it)=-1.d0
               sintarr(it)=0.d0
            else
               costarr(it)=cos(thetarr(it))
               sintarr(it)=sin(thetarr(it))
            endif
            if(thetarr(it).gt.(pihalf-eps).and.thetarr(it).lt.
     1           (pihalf+eps)) 
     1           then
               tan2arr(it)=-1.d0
            else
               tan2arr(it)=tan(thetarr(it))**2
            endif
c            print*,'thetarr,costarr,tan2arr '
c            print*,'thetarr',thetarr(it)*180.d0/pi
         enddo
      else
         thetarr(1)=0.d0
      endif

      open(unit=15,file='tharr.dat',status='unknown')
      write(15,*) ntg,' = number of grid points in theta'
      write(15,*) 
     1     'index  thet(rad)  thet(deg)    cost     tan**2(thet)'
      do it=1,ntg
         write(15,*) it,sngl(thetarr(it)),sngl(thetarr(it)*180.d0/pi),
     1        sngl(costarr(it)),sngl(tan2arr(it))
      enddo
      close(15)

      if (npg.gt.1) then
         dp=r2p/dfloat(npg-1)
         do ip=1,npg
            phiarr(ip)=dp*(dfloat(ip-1))
            aarr(ip)=dsin(phiarr(ip))
            barr(ip)=-dcos(phiarr(ip))
c     print*,'phiarr, aarr, barr ',phiarr(ip),aarr(ip),barr(ip)
c     carr(ip)=
c     darr(ip)=
         enddo
      else
         dp=r2p
         ip=1
         phiarr(ip)=0.d0
         aarr(ip)=dsin(phiarr(ip))
         barr(ip)=-dcos(phiarr(ip))
      endif

      open(unit=15,file='phiarr.dat',status='unknown')
      write(15,*) npg
      do ip=1,npg
         write(15,*) phiarr(ip)*180.d0/pi
      enddo
      close(15)

c     calculate density in grid
      massenv=0.d0
      massdisk=0.d0
      do ir=1,nrg-1
         rad=0.5d0*(rarr(ir)+rarr(ir+1))
         dr=rarr(ir+1)-rarr(ir)
         if (ntg.gt.1) then     !2- or 3-D atmosphere
            do it=1,ntg-1
               thet=0.5d0*(thetarr(it)+thetarr(it+1))
               cost=cos(thet)
               sint=sin(thet)
               dt=thetarr(it+1)-thetarr(it)
               if (npg.gt.1) then
                  do ip=1,npg-1   !3-D atmosphere
                     phi=0.5d0*(phiarr(ip)+phiarr(ip+1))
                     dp=phiarr(ip+1)-phiarr(ip)
                     call densenv(rad,thet,cost,sint,pi,pihalf,phi
     1                    ,dense)
                     vol=rad**2*sint*dt*dp*dr*rmincgs**3
                     massenv=massenv+dense*vol/rmincgs 
                     call densdisk(rad,sint,cost,phi,densd)
                     massdisk=massdisk+densd*vol
                     densarr(ir,it,ip)=dense+densd
                     massarr(ir,it,ip)=(dense+densd)*vol
                  end do
               else           !2-D atmosphere
                  ip=1       
                  phi=0.d0
                  dp=r2p
                  call densenv(rad,thet,cost,sint,pi,pihalf,phi
     1                 ,dense)
                  vol=rad**2*sint*dt*dp*dr*rmincgs**3
                  massenv=massenv+dense*vol/rmincgs
                  call densdisk(rad,sint,cost,phi,densd)
                  massdisk=massdisk+densd*vol
                  densarr(ir,it,ip)=dense+densd
                  massarr(ir,it,ip)=(dense+densd)*vol
               endif
            enddo
         else              !1-D atmosphere
            phi=0.d0
            ip=1
            it=1
            thet=0.d0
            cost=1.d0
            sint=0.d0
            call densenv(rad,thet,cost,sint,pi,pihalf,phi,dense)
            vol=rad**2*dr*2.d0*r2p*rmincgs**3
            massenv=massenv+dense*vol/rmincgs
            call densdisk(rad,sint,cost,phi,densd)
            massdisk=massdisk+densd*vol
            densarr(ir,it,ip)=dense+densd
            massarr(ir,it,ip)=(dense+densd)*vol
c     print*,'densarr',densarr(ir,it,ip)
         endif
      end do
      massenv=massenv/1.989d33/rmincgs    !extra factor in dense
      massdisk=massdisk/1.989d33/rmincgs
c     important note:  densarr(ir,it,ip) is the density halfway between
c     ir and ir+1, it and it+1, ip and ip+1; it is NOT the density at
c     ir,it,ip.

      print*, 'massenv (in solar units) ',massenv
      print*, 'massdisk from grid, compared to input ',massdisk,
     $     massd/1.989d33
      print*, ' '

c     write out rarr,thetarr,densarr at phi=0 so you can read them into
c     IDL.


c     integrate optical depth along the theta angles
c     at phi=0.
      open(unit=15,file='tau.dat',status='unknown')
      write(15,*) 
     1 '    theta          tau_env        tau_disk       tau_tot'
      do it=1,ntg
         tau=0.d0
         taud=0.d0
         taue=0.d0
         if (ntg.eq.1) then
            thet=60.d0*pi/180.d0
         else
            thet=thetarr(it)
         endif
         if (ntg.gt.1.and.it.eq.(ntg/2)+1) then
c     TSC blows up at thet=90
            print*,'thet',thet*180./pi
            thet=89.999d0*pi/180.d0
         endif
         cost=cos(thet)
         sint=sin(thet)
         phi=0.d0
         do ir=1,nrg-1
            dr=rarr(ir+1)-rarr(ir)
            rad=0.5d0*(rarr(ir)+rarr(ir+1))
            call densenv(rad,thet,cost,sint,pi,pihalf,phi,dense)
            call densdisk(rad,sint,cost,phi,densd)
            call locate(cosbarr,MXTHETB,nthetab,cost,iinc)
            kappa=linterp(cosbarr,kapd,MXTHETB,cost,iinc)
            taud=taud+kappa*(densd)*dr
            taue=taue+kappa*(dense)*dr
         end do
         write(15,*) sngl(thet*180.d0/pi),sngl(taue)
     1        ,sngl(taud),sngl(taue+taud)
      end do
      close(15)

c     calculate Av along ethet,ephi directions (peeled image)
c     assume Kappa_V = 224.8
c      open(unit=15,file='A_v.dat',status='unknown')
c      write(15,*) 
c     1 '    theta           Av_env         Av_disk        Av_tot'
c      thet=thete
c      cost=cos(thet)
c      sint=sin(thet)
c      phi=0.d0
c      do ir=1,nrg-1
c         dr=rarr(ir+1)-rarr(ir)
c         rad=0.5d0*(rarr(ir)+rarr(ir+1))
c         call densenv(rad,thet,cost,sint,pi,pihalf,phi,dense)
c         call densdisk(rad,sint,cost,phi,densd)
c         taud=taud+224.8*rmincgs*(densd)*dr
c         taue=taue+224.8*rmincgs*(dense)*dr
c      end do
c      write(15,*) sngl(thet*180.d0/pi),sngl(taue*1.086)
c     1     ,sngl(taud*1.086),sngl((taue+taud)*1.086)
c      close(15)

c     write out density array to read in IDL
c     make the cartesian grid here (fortran is faster than IDL)
c     this makes a grid in the x-z plane
      n=200
c     xmax=sngl(rarr(nrg))
c     xmax=rmaxd
      xmax=sngl(xymaxdens)
      print*,'xmax',xmax
      dx=xmax/float(n)
      do ix=1,n
         x=dx*(float(ix)-0.5)
      do iy=1,n
         y=dx*(float(iy)-0.5)
         if (iy.eq.1.and.ix.eq.1) print*,'x0,y0',x,y
         r=sqrt(x**2+y**2)
         thet=atan2(y,x)
         call locate(rarr,nrg,dble(r),ir)
         if (ntg.gt.1) then
            call locate(thetarr,ntg,(thet),it)
         else
            it=1
         endif
         xyarr(ix,iy)=sngl(densarr(ir,it,1))
      end do
      end do
      open(unit=15,file='densxz.dat',status='unknown')
      do ix=1,n
         write(15,*) (xyarr(ix,iy),iy=1,n)
      enddo
      close(15)

c     write out array in x-y plane at z=.5*rmax to see phi dependence
      thet=pi/2.d0
      if (ntg.gt.1) then
         call locate(thetarr,ntg,(thet),it)
c     print*,'it',it
      else
         it=1
      endif
      do ix=1,n
         x=dx*(float(ix)-0.5)
      do iy=1,n
         y=dx*(float(iy)-0.5)
         if (iy.eq.1.and.ix.eq.1) print*,'x0,y0',x,y
         r=sqrt(x**2+y**2)
         phi=atan2(y,x)
         call locate(rarr,nrg,dble(r),ir)
c     print*,'ir',ir
         if (npg.gt.1) then
            call locate(phiarr,npg,(phi),ip)
         else
            ip=1
         endif
         xyarr(ix,iy)=sngl(densarr(ir,it,ip))
      end do
      end do
      open(unit=15,file='densxy.dat',status='unknown')
      do ix=1,n
         write(15,*) (xyarr(ix,iy),iy=1,n)
      enddo
      close(15)

      end


      






