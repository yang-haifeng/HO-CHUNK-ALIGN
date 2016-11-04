c     *********************************************************

      subroutine propagate(iphot,sini,xmax,fractx,dum,eps)

      implicit none

      include 'stokes.txt'
      include 'tts.txt'
      include 'opacin.txt'
      include 'taunum.txt'
      include 'random.txt'
      include 'vger.txt'
      include 'tab.txt'

      real*8 xran,xpold,ypold,rpold,rho2,zpold,x,y,xmax,fractx
     1   ,sipold,tsum,pol2,eps,alb,coskapd,kap,kappa_q,kap_ext
     1   ,kapsI,kapsQ,kaps
      real*8 sini(15)
      real ran2
      integer ii(3),dum(3),iphot,iscat,k,ia,ithet,ix,iy
     1   ,icount,i,ir,it,ip,iinc
      real*8 linterp

      xpold=xp
      ypold=yp
c     rpold=rp
c     rhopold=dsqrt(rp**2+zp**2)
      zpold=zp
c     already calculated rsq and rtot before call to this routine.

      do i=1,3
         ii(i)=dum(i)
c     print*,'ii ',ii(i)
      enddo

c     iflag is set to 1 if photon scatters
      iflag=0
      exitflag=0
      aflag=0
      iscat=0

c     sample tau
      xran=ran2(i1)
      tau=-dlog(xran)	
c     integrate over distance until the optical depth equals tau.
      call tauint(iphot,tsum,ii,r2p,sip,sqp,sup,svp,eps)

      if(exitflag.eq.1) go to 300
      if(aflag.eq.1) then
         write(6,*) 'shouldnt be here, aflag=1'
         go to 400
      end if
      tot=tot+1.d0
      
c     photon scatters until exit exits disk
 30   continue
      iflag=1
c     sip=sip*rlam
c     sqp=sqp*rlam
c     sup=sup*rlam
c     svp=svp*rlam
      xran=ran2(i1)
      coskapd=(cost)
      call locate(cosbarr,MXTHETB,nthetab,coskapd,iinc)
      kap=linterp(cosbarr,kapd,MXTHETB,uz,iinc)
      kappa_q=linterp(cosbarr,kapd_q,MXTHETB,uz,iinc)
      kap_ext=kap+kappa_q*sqp/sip
      kapsI=linterp(cosbarr,kaps_i,MXTHETB,uz,iinc)
      kapsQ=linterp(cosbarr,kaps_q,MXTHETB,uz,iinc)
      kaps=kapsI+kapsQ*sqp/sip
      alb=kaps/kap_ext
c      print*,'alb,kap_ext,Q/I,cost',alb,kap_ext,sqp/sip,cost
      if(xran.le.alb) then
          
         iscat=iscat+1
c     if (sip.lt.(1.0d-3*rlam)) then
c     aflux=aflux+1
c     go to 400
c     endif
         sipold=sip
         if (ipeel.eq.1) then
             call peeloff(xp,yp,zp,sip,sqp,sup,svp
     1           ,cost,sint,cosp,sinp,phi
     1           ,pi,r2p,hit,htot,rsq,rtot
     1           ,peak,tsum,ii,idust,iflag,iphot,eps,kaps)
         endif
c     NOTE**** need to write vger_3d.f subroutine*******
         if(iveeg.eq.1) then
            call vger(xp,yp,zp,sip,sqp,sup,svp
     1           ,cost,sint,cosp,sinp,phi
     1           ,pi,r2p,hit,htot,rsq,rtot
     1           ,peak,tsum,ii,idust,iflag,iphot,eps)
         endif
         
         pol2=sqp**2+sup**2+svp**2
         if (pol2 .gt. sip**2) then
            print*,'error, P^2, I^2 ', pol2,sip**2
            continue
         endif
         
c     xran=ran2(i1)
c     if(xran.le.rlam) then
c     if (idust.lt.2) then
c     call stokes
c     else
         call stokes16
c     end if
         pol2=sqp**2+sup**2+svp**2
         if (pol2 .gt. sip**2) then
            print*,'error, P^2, I^2 ', pol2,sip**2
            continue
         endif
c     iscat=iscat+1
      else
c     photon absorbed, reemitted at longer wavelengths which we won't
c     consider here
         if(alb.eq.1.) write(6,*)'error! rlam=1, yet photon absorbed'
         aflux=aflux+1.
         go to 400
c     call therm
      end if
c     sample tau
      xran=ran2(i1)
      tau=-dlog(xran)
c     integrate distance until optical depth equals tau.
      xpold=xp
      ypold=yp
      zpold=zp
c     not sure if we need to do this here, might already have this
c      rsq=xp**2+yp**2+zp**2
c      rtot=dsqrt(rsq)
c     rpold=rp
c     rhopold=dsqrt(rp**2+zp**2)
      ux=sint*cosp
      uy=sint*sinp
      uz=cost
      call tauint(iphot,tsum,ii,r2p,sip,sqp,sup,svp,eps)
      ir=ii(1)
      it=ii(2)
      ip=ii(3)
      if(exitflag.eq.1) go to 300
      if (aflag.eq.1) then
c     write(6,*) 'photon absorbed by star, iphot',iphot
         go to 400
      end if 
      go to 30
      
 300  continue
c     photon exits
c     bin angle
c     NOTE:  assuming axisymmetric, and z=-z.  
c     k=int(float(nmu-1)*dabs(cost)+1.5d0)
c     NOTE:  no longer assuming z=-z with toroidal B-field
      k=int(float(nmu)*(1.d0-cost)/2.d0)+1
c      k=min(int(float(nmu)*dabs(cost)+1),nmu)
      if(k.lt.0.or.k.gt.nmu) then
         print*,'cost binning error, k, nmu',k,nmu
      endif

      x = ypold*cosp-xpold*sinp
      y = zpold*sini(k)-ypold*u(k)*sinp
     &     -xpold*u(k)*cosp
      
c     first, sum fluxes
      rho2=(x**2+y**2)
      if (rho2.lt.aperture2(3)) then
         flux=flux+1
         if (iflag.eq.1) then
            sflux=sflux+1
            nscat=nscat+iscat
         endif
      endif
      do ia=1,3
         if (rho2.lt.aperture2(ia)) then
c     flux=flux+1
            si(k,ia)=si(k,ia)+sip
            sq(k,ia)=sq(k,ia)+sqp
            su(k,ia)=su(k,ia)+sup
            sv(k,ia)=sv(k,ia)+svp
            si2(k,ia)=si2(k,ia)+sip*sip
            sq2(k,ia)=sq2(k,ia)+sqp*sqp
            su2(k,ia)=su2(k,ia)+sup*sup
            sv2(k,ia)=sv2(k,ia)+svp*svp
            nums(k,ia)=nums(k,ia)+1
            if (iflag.eq.1) then
               spi(k,ia)=spi(k,ia)+sip
               spq(k,ia)=spq(k,ia)+sqp
               spu(k,ia)=spu(k,ia)+sup
               spv(k,ia)=spv(k,ia)+svp
            end if
         end if
      end do
      
c     ithet=int(acos(cost)*fractt+0.5)
c      ithet=int(float(nmu-1)*dabs(cost)+1.5d0)
c      if(ithet.lt.1.or.ithet.gt.nmu) then
c         print*, 'theta wrong, ithet',ithet
c         print*, 'cost,sint',cost,sint
c         print*, 'iphot',iphot
c      end if
      
c      if (cost.lt.0.d0) then
c         x = -ypold*cosp+xpold*sinp
c         y = -zpold*sini(ithet)-ypold*u(ithet)*sinp
c     &        -xpold*u(ithet)*cosp
c      else
c         x = ypold*cosp-xpold*sinp
c         y = zpold*sini(ithet)-ypold*u(ithet)*sinp
c     &        -xpold*u(ithet)*cosp
c      endif
      
c     testing**************
c     if (cost.ge.0) then
      
      ix=int((x+xmax)/fractx)+1
      iy=int((y+xmax)/fractx)+1
      if (ix.le.nx.and.iy.le.nx.and.ix.gt.0.d0.and.iy.gt.0.d0) then
         image(ix,iy,k)=image(ix,iy,k)+(sip)
         imagei(ix,iy,k)=imagei(ix,iy,k)+(sip)
         imageq(ix,iy,k)=imageq(ix,iy,k)+(sqp)
         imageu(ix,iy,k)=imageu(ix,iy,k)+(sup)
         imagev(ix,iy,k)=imagev(ix,iy,k)+(svp)
         image2(ix,iy,k)=image2(ix,iy,k)+sip*sip
         numi(ix,iy,k)=numi(ix,iy,k)+1
      else
         icount=icount+1
      endif
      
      if(iflag.eq.0) then
         star(k)=star(k)+sip
         starq(k)=starq(k)+sqp
         staru(k)=staru(k)+sup
         starv(k)=starv(k)+svp
      end if
      
 400  continue
      return
      end
      
      
c     *********************************************************
