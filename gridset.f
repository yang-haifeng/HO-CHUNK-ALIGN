      subroutine gridset

c  set up grid:  3-D spherical grid with variable spacing in r and 
c  theta (set by exponents, rexp, texp)
c  make a bunch of arrays necessary for find_wall subroutine

c history:
c     2/22/02 BAW: arbritrary B-field orientation
c 00/09/06 (baw): write subroutine
	
      implicit none
      include 'tts.txt'
      include 'out.txt'
      include 'grid.txt'
      include 'opacin.txt'
      include 'tab.txt'

      real*8 dr,dt,dp,rad,phi,thet,pi,r2p,densd,dense
     $ ,cost,sint,pihalf,eps,tau,rming1,rming2,rming3,tiny
     $     ,rmincgs,taud,taue,rfac,thettop,thetatm,mu,res,xp,yp,zp
     $     ,thetmin,kappa,mu0,sinpbz,cospbz,phibz,sintbz,costbz
     1     ,cospb,sinpb,phib,sintb,costb,sup,sqp,cosp,sinp
     1     ,ga,gb,gc,gs,pwr
      real xyarr(200,200),x,y,xmax,dx,r
      integer ir,it,ip,nttmp,gflag,n,ix,iy,irbeg,nrwall,ntatm,iinc,flag

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

c     read in galli-shu-toroid constants
c     open(unit=15,file='gs.in',status='unknown')
c       read(15,*) gs,ga,gb,gc,pwr
c     close(15)

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

      if (massd.eq.0.d0) then
c     since there's empty space below rmine and rddust....
         dr=(rmax-rming1)/(dfloat(nrg-2))**rexp
         print*,'rexp,dr,rmin',rexp,dr/autors,rmin/autors
         do ir=3,nrg
            rarr(ir)=rming1+dr*(dfloat(ir-2))**rexp
c     print*,'rarr ',ir,rarr(ir)/au
         enddo

      else
         
c     finely sample disk to radial optical depth taux (computed
c     in setup_wave.f)
         nrwall=irbeg+int(taux/1.d0)+1     !tau=1 steps
         print*,'nrwall',nrwall
         if (nrwall.gt.nrg) then
            write(*,*) 'too many points in rgrid'
            stop
         endif
         print*,'rdx,rddust',rdx,rddust
         rfac=(dlog(rdx)/dlog(rddust))**(1.d0/float(nrwall-irbeg))
         print*,'rfac',rfac
         rarr(irbeg)=rddust
         do ir=irbeg,nrwall-1
            rarr(ir+1)=rarr(ir)**rfac
         end do
         print*,'irbeg,nrwall',irbeg,nrwall

         rming1=rarr(nrwall)
         dr=(rmax-rming1)/(dfloat(nrg-nrwall))**rexp
         print*,'rexp,dr,rmin',rexp,dr/autors,rming1
         do ir=nrwall+1,nrg
            rarr(ir)=rming1+dr*(dfloat(ir-nrwall))**rexp
         enddo

      endif

c     set rmine to a grid location
      if (gflag.eq.1) then
         call locate(rarr,nrg,nrg,rmine,ir)
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

c     r-ave array
      do ir=1,nrg-1
         ravearr(ir)=0.5*(rarr(ir)+rarr(ir+1))
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
c     make spacing fine for, say first 10 scale heights
            mu=10.d0*z1*rddust**b/rddust
            thetatm=asin(mu)
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
            print*,'ntatm',ntatm

            mu0=z1*rddust**b/rddust
            print*,'thet of 1 scale height at disk',
     1           sngl(asin(mu0)*180./pi),' degrees'
            thetmin=min(0.1*asin(mu0),thetatm/float(ntatm))

c     jon's log spacing.
c            tmptharr(1)=0.d0
c            do it=2,ntatm
c               tmptharr(it)=thetmin *
c     1              (thetatm/thetmin)**((it-2.d0)/(ntatm-2.d0))
cc     print*,'tmptharr ',tmptharr(it)*180.d0/pi
c            enddo

c     barb's equal spacing
            do it=2,ntatm
               tmptharr(it)=float(it-1)*thetatm/float(ntatm-1) 
            enddo
            
         endif
         do it=ntatm+1,nttmp-1
            tmptharr(it)=thetatm+
     1           (it-ntatm)*(pihalf-thetatm-thettop)/(nttmp-ntatm)
c     print*,'tmptharr(it)',tmptharr(it)
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

c     theta-ave array
      do it=1,ntg-1
         thetavearr(it)=0.5*(thetarr(it)+thetarr(it+1))
      enddo

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

c     phi-ave array
      do ip=1,npg-1
         phiavearr(ip)=0.5*(phiarr(ip)+phiarr(ip+1))
      enddo

      open(unit=15,file='phiarr.dat',status='unknown')
      write(15,*) npg
      do ip=1,npg
         write(15,*) phiarr(ip)*180.d0/pi
      enddo
      close(15)

c     calculate density and B-field in grid
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
                  do ip=1,npg-1 !3-D atmosphere
                     phi=0.5d0*(phiarr(ip)+phiarr(ip+1))
                     dp=phiarr(ip+1)-phiarr(ip)
                     call densenv(rad,thet,cost,sint,pi,pihalf,phi
     1                    ,dense)
                     volarr(ir,it,ip)=rad**2*sint*dt*dp*dr*rmincgs**3
                     massenv=massenv+dense*volarr(ir,it,ip) 
                     call densdisk(rad,sint,cost,phi,densd)
                     massdisk=massdisk+densd*volarr(ir,it,ip)
                     densarr(ir,it,ip)=(dense+densd)*rmincgs
                     massarr(ir,it,ip)=(dense+densd)*volarr(ir,it,ip)
c                    call bset(rad/rmax,rc/rmax,cost,sint,phi,r2p,pi,
c    1                    costbz,sintbz,phibz,cospbz,sinpbz,gs,ga,gb,
c    1                    gc,pwr)
c                    costbzarr(ir,it,ip)=costbz
c                    sintbzarr(ir,it,ip)=sintbz
c                    phibzarr(ir,it,ip)=phibz
c                    cospbzarr(ir,it,ip)=cospbz
c                    sinpbzarr(ir,it,ip)=sinpbz
                     costbzarr(ir,it,ip)=1.
                     sintbzarr(ir,it,ip)=0.
                     phibzarr(ir,it,ip)=0.
                     cospbzarr(ir,it,ip)=1.
                     sinpbzarr(ir,it,ip)=0.
                  end do
               else             !2-D atmosphere
                  ip=1       
                  phi=0.d0
                  dp=r2p
                  call densenv(rad,thet,cost,sint,pi,pihalf,phi
     1                 ,dense)
                  volarr(ir,it,ip)=rad**2*sint*dt*dp*dr*rmincgs**3
                  massenv=massenv+dense*volarr(ir,it,ip)
                  call densdisk(rad,sint,cost,phi,densd)
                  massdisk=massdisk+densd*volarr(ir,it,ip)
                  densarr(ir,it,ip)=(dense+densd)*rmincgs
                  massarr(ir,it,ip)=(dense+densd)*volarr(ir,it,ip)
                  call bset(rad/rmax,rc/rmax,cost,sint,phi,r2p,pi,
     1                 costbz,sintbz,phibz,cospbz,sinpbz,gs,ga,gb,
     1                 gc,pwr)
                  costbzarr(ir,it,ip)=costbz
                  sintbzarr(ir,it,ip)=sintbz
                  phibzarr(ir,it,ip)=phibz
                  cospbzarr(ir,it,ip)=cospbz
                  sinpbzarr(ir,it,ip)=sinpbz               
               endif
            enddo
         else                   !1-D atmosphere
            phi=0.d0
            ip=1
            it=1
            thet=0.d0
            cost=1.d0
            sint=0.d0
            call densenv(rad,thet,cost,sint,pi,pihalf,phi,dense)
            volarr(ir,it,ip)=rad**2*dr*2.d0*r2p*rmincgs**3
            massenv=massenv+dense*volarr(ir,it,ip)
            call densdisk(rad,sint,cost,phi,densd)
            massdisk=massdisk+densd*volarr(ir,it,ip)
            densarr(ir,it,ip)=(dense+densd)*rmincgs
c     units of distance in code are in rmincgs
            massarr(ir,it,ip)=(dense+densd)*volarr(ir,it,ip)
c     print*,'densarr',densarr(ir,it,ip)
            call bset(rad/rmax,rc/rmax,cost,sint,phi,r2p,pi,
     1           costbz,sintbz,phibz,cospbz,sinpbz,gs,ga,gb,gc
     1           ,pwr)
            costbzarr(ir,it,ip)=costbz
            sintbzarr(ir,it,ip)=sintbz
            phibzarr(ir,it,ip)=phibz
            cospbzarr(ir,it,ip)=cospbz
            sinpbzarr(ir,it,ip)=sinpbz            
         endif
      end do
      massenv=massenv/1.989d33
      massdisk=massdisk/1.989d33
c     important note:  densarr(ir,it,ip) is the density halfway between
c     ir and ir+1, it and it+1, ip and ip+1; it is NOT the density at
c     ir,it,ip.
      
      print*, 'massenv (in solar units) ',massenv
      print*, 'massdisk from grid, compared to input ',massdisk,
     $     massd
      print*, ' '
      write(12,*) 'massenv (in solar units) ',massenv
      write(12,*) 'massdisk from grid, compared to input ',massdisk,
     $     massd
      
c     write out rarr,thetarr,densarr at phi=0 so you can read them into
c     IDL.
      
      
c     integrate optical depth along the theta angles
c     at phi=0.
      open(unit=15,file='tau.dat',status='unknown')
      write(15,*) 
     1     '   theta     tau_env     tau_disk     tau_tot'
      ip=1
      sqp=0.
      sup=0.
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
            flag=1
         else
            flag=0
         endif
         cost=cos(thet)
         sint=sin(thet)
         phi=phiarr(ip)
         cosp=cos(phi)
         sinp=sin(phi)
         do ir=1,nrg-1
            dr=rarr(ir+1)-rarr(ir)
            rad=0.5d0*(rarr(ir)+rarr(ir+1))
            call densenv(rad,thet,cost,sint,pi,pihalf,phi,dense)
            call densdisk(rad,sint,cost,phi,densd)
c     rotate to B-field frame, calc. opacity in this frame
            xp=rad*sint*cosp
            yp=rad*sint*sinp
            zp=rad*cost
            
            call rotate(ir,it,ip,pi,r2p,sqp,sup,sqp,sup
     1           ,cost,sint,phi,cosp,sinp
     1           ,costb,sintb,phib,cospb,sinpb)
            call locate(cosbarr,MXTHETB,nthetab,costb,iinc)
            kappa=linterp(cosbarr,kapd,MXTHETB,costb,iinc)
            taud=taud+kappa*(densd)*dr*rmincgs
            taue=taue+kappa*(dense)*dr*rmincgs
c     print this out to see dtaud---log(R) spacing works!
c            if (flag.eq.1) then
c               print*,'r,dtaud in equator',rarr(ir),
c     1              kappa*(densd)*dr*rmincgs
c            endif
         end do
         write(15,901) sngl(thet*180.d0/pi),sngl(taue)
     1        ,sngl(taud),sngl(taue+taud)
      end do
      close(15)
901     format((f10.5),3(1x,1pe12.5))


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
         call locate(rarr,nrg,nrg,dble(r),ir)
         if (ntg.gt.1) then
            call locate(thetarr,ntg,ntg,(thet),it)
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
         call locate(thetarr,ntg,ntg,(thet),it)
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
         call locate(rarr,nrg,nrg,dble(r),ir)
c     print*,'ir',ir
         if (npg.gt.1) then
            call locate(phiarr,npg,npg,(phi),ip)
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

c     write out density array
      open(unit=15,file='dens.unf',status='unknown',form='unformatted')
      write(15) nrg-1,ntg-1,npg-1
      write(15) (ravearr(ir)/autors,ir=1,nrg-1)
      write(15) (thetavearr(it),it=1,ntg-1)
      write(15) (phiavearr(ip),ip=1,npg-1)
      write(15) (((densarr(ir,it,ip)/rmincgs,ir=1,nrg-1),it=1,ntg-1)
     1     ,ip=1,npg-1)
      close(15)

c     stop
      end


      







