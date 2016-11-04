      subroutine gridset

c  set up grid:  3-D spherical grid with variable spacing in r and theta
c  (set by exponents, rexp, texp
c  make a bunch of arrays necessary for find_wall subroutine

c history:
c 00/09/06 (baw): write subroutine
	
      implicit none
      include 'mars.txt'
      include 'taunum.txt'
      include 'grid.txt'

      real*8 tmpdens(nsig),dr,rad,thet,phi,dense,t,bnu,dv,km2cm
     $ ,cost,sint,pihalf,pi,r2p,rmin,rmax,dt,eps,dp,massatm,massco2
      integer ir,it,ip,nttmp,is

      pi=4.d0*datan(1.d0)
      r2p=2.d0*pi
      pihalf=pi/2.d0
      km2cm=1.d5
      massco2=14.d0*2.d0*16.d0*1.6735344d-24

c     make grid include minimum and maximum values.

      print*, ' '
      print*, 'grid setup'

c     moved rarr initialization to readatm.f (to calculate optical depths
c     consistently).

c     r-squared array
      do ir=1,nrg
         r2arr(ir)=rarr(ir)**2
c     print first few points of rarr
c     if (ir.lt.10) then
c     print*,'rarr/rmin,ir ',rarr(ir)/rmin,ir
c     endif
      enddo

      if (ntg.gt.1) then
c     set up a tmp array going from 0 - 90 from equ. to pole
c     if rexp > 1, spacing will increase from equ. to pole
         nttmp=(ntg+1)/2
         pihalf=pi/2.d0
         dt=pihalf/(dfloat(nttmp-1))*texp
c     check: say ntg=19, nttmp=10, if texp=1, dt=10, correct
         do it=1,nttmp
            tmptharr(it)=dt*(dfloat(it-1))**texp
c     print*,'tmptharr ',tmptharr(it)
         enddo
         do it=1,nttmp
            thetarr(it)=pihalf-tmptharr(nttmp-it+1)
       if ((it+nttmp).le.ntg) thetarr(it+nttmp)=tmptharr(it)+pihalf+dt
         enddo
         eps=pihalf/100.
         do it=1,ntg
            costarr(it)=cos(thetarr(it))
       if(thetarr(it).gt.(pihalf-eps).and.thetarr(it).lt.(pihalf+eps)) 
     1           then
               tan2arr(it)=-1.d0
            else
               tan2arr(it)=tan(thetarr(it))**2
            endif
c     print*,'thetarr,costarr,tan2arr '
c     print*,thetarr(it),costarr(it),tan2arr(it)
         enddo
      endif

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
         ip=1.d0
         phiarr(ip)=0.d0
         aarr(ip)=dsin(phiarr(ip))
         barr(ip)=-dcos(phiarr(ip))
      endif

c     now calculate density in grid
      massatm=0.d0
      do ir=1,nrg-1
         rad=0.5d0*(rarr(ir)+rarr(ir+1))
         dr=(rarr(ir+1)-rarr(ir))*km2cm
         if (ntg.gt.1) then              !2- or 3-D atmosphere
            do it=1,ntg-1
               thet=0.5d0*(thetarr(it)+thetarr(it+1))
               cost=cos(thet)
               sint=sin(thet)
               dt=thetarr(it+1)-thetarr(it)
               if (npg.gt.1) then             !3-D atmosphere
                  do ip=1,npg-1
                     phi=0.5d0*(phiarr(ip)+phiarr(ip+1))
                     dp=phiarr(ip+1)-phiarr(ip)
                     call densatm(rad,tmpdens) !not set up for 3-D
                     do is=1,nsig
                        dense=tmpdens(is)
                        densarr(ir,it,ip,is)=dense
                     enddo
                     dv=(rad*km2cm)**2*dr*sint*dt*dp 
                     massatm=massatm+densarr(ir,it,ip,5)*dv*massco2
                  end do
               else
                  ip=1
                  dp=r2p
                  call densatm(rad,tmpdens) !not set up for 2-D
                  do is=1,nsig
                     dense=tmpdens(is)
                     densarr(ir,it,ip,is)=dense
                  enddo
                  dv=(rad*km2cm)**2*dr*sint*dt*dp 
                  massatm=massatm+densarr(ir,it,ip,5)*dv*massco2
               endif
            enddo
         else
            ip=1
            it=1
            call densatm(rad,tmpdens)
            do is=1,nsig
               dense=tmpdens(is)
               densarr(ir,it,ip,is)=dense
            enddo
            dv=(rad*km2cm)**2*dr*r2p*2.d0
c     print*,'rad,dr,dv ',rad*km2cm,dr,dv
c     print*,'r2p',r2p
            massatm=massatm+densarr(ir,it,ip,5)*dv*massco2
         endif
      end do
      massatm=massatm/1.989d33
c     important note:  densarr(ir,it,ip) is the density halfway between
c     ir and ir+1, it and it+1, ip and ip+1; it is NOT the density at
c     ir,it,ip.

      print*, 'mass of CO2 in atm (in solar units) ',massatm

c     write out rarr,thetarr,densarr at phi=0 so you can read them into
c     IDL.

c     integrate optical depths along normal to atmosphere
c     at phi=0.
      tdust=0.d0
      tco2=0.d0
      tcloud=0.d0
      to3=0.d0
      tnot=0.d0
      it=1
      ip=1
      do ir=1,nrg-1
         dr=(rarr(ir+1)-rarr(ir))   
         rad=0.5d0*(rarr(ir)+rarr(ir+1))
         tdust=tdust+(sigma(1)+sigma(2))*densarr(ir,it,ip,1)*dr
         tcloud=tcloud+(sigma(3)+sigma(4))*densarr(ir,it,ip,3)*dr
         to3=to3+(sigma(7)+sigma(8))*densarr(ir,it,ip,7)*dr
         tco2=tco2+(sigma(5)+sigma(6))*densarr(ir,it,ip,5)*dr
      end do

      print*,'in subroutine gridset'
      print*,'dust optical depth ',tdust
      print*,'CO2 optical depth ',tco2
      print*,'cloud optical depth ',tcloud
      print*,'O3 optical depth ',to3
      tnot=tdust+tcloud+tco2+to3
      print*,'total optical depth ',tnot

c     now calculate energy emitted by each grid cell
      if (itherm.eq.'Y') then
         do ir=1,nrg-1
            rad=0.5d0*(rarr(ir)+rarr(ir+1))
            dr=(rarr(ir+1)-rarr(ir))*km2cm
            t=tarr(ir)
c     if at some future time, t is 3-D, modify this
            if (ntg.gt.1) then  !2- or 3-D atmosphere
               do it=1,ntg-1
                  thet=0.5d0*(thetarr(it)+thetarr(it+1))
                  sint=sin(thet)
                  dt=thetarr(it+1)-thetarr(it)
                  if (npg.gt.1) then !3-D atmosphere
                     do ip=1,npg-1
                        phi=0.5d0*(phiarr(ip)+phiarr(ip+1))
                        dp=phiarr(ip+1)-phiarr(ip)
                        dv=(rad*km2cm)**2*dr*sint*dt*dp
                        call plancknu(t,wnu,bnu)
                        earr(ir,it,ip)=bnu*dv*(densarr(ir,it,ip,2)*
     1                       sigma(2)+densarr(ir,it,ip,4)*sigma(4)
     2                       +densarr(ir,it,ip,6)*sigma(6)
     3                       +densarr(ir,it,ip,8)*sigma(8))
c     that's the absorption opacities
                     end do
                  else
                     dp=r2p
                     ip=1
                     dp=phiarr(ip+1)-phiarr(ip)
                     dv=(rad*km2cm)**2*dr*sint*dt*dp
                     call plancknu(t,wnu,bnu)
                     earr(ir,it,ip)=bnu*dv*(densarr(ir,it,ip,2)*
     1                    sigma(2)+densarr(ir,it,ip,4)*sigma(4)
     2                    +densarr(ir,it,ip,6)*sigma(6)
     3                    +densarr(ir,it,ip,8)*sigma(8))
                  endif
               enddo
            else
               ip=1
               it=1
               dv=(rad*km2cm)**2*dr*r2p*2.d0
               call plancknu(t,wnu,bnu)
               earr(ir,it,ip)=bnu*dv*(densarr(ir,it,ip,2)*
     1              sigma(2)+densarr(ir,it,ip,4)*sigma(4)
     2              +densarr(ir,it,ip,6)*sigma(6)
     3              +densarr(ir,it,ip,8)*sigma(8))/km2cm 
c     note, sigma is in units of 1/km, hence the factor of km2cm
               if (ir.eq.1) print*,'bnu of lowest level atm',bnu
     1              ,'rad,dr',rad,dr
c               print*,'ir,bnu,dv,densarr,earr',ir,bnu,dv,
c     1              densarr(ir,it,ip,2),earr(ir,it,ip)
            endif
         end do
      endif

c     stop
      end


      







