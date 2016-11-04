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
     1   ,kapsI,kapsQ,kaps,supb,sqpb,sinpb,cospb,phib,sintb,costb
      real*8 sini(nmu)
      real ran2
      integer ii(3),dum(3),iphot,iscat,k,ia,ix,iy
     1   ,icount,i,ir,it,ip,iinc

      real*8 linterp

      xpold=xp
      ypold=yp
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
      call tauint(iphot,tsum,ii,r2p,sip,sqp,sup,svp,eps
     1     ,cost,sint,phi,cosp,sinp,pi
     1     ,costb,sintb,phib,cospb,sinpb,sqpb,supb)
      if(exitflag.eq.1) go to 300
      if(aflag.eq.1) then
         write(6,*) 'shouldnt be here, aflag=1'
         go to 400
      end if
      tot=tot+1.d0
      
c     photon scatters until exit exits disk
 30   continue
      iflag=1

      xran=ran2(i1)
c     use rotated stokes and B-angles calculated in tauint

      call locate(cosbarr,MXTHETB,nthetab,costb,iinc)
      kap=linterp(cosbarr,kapd,MXTHETB,costb,iinc)
      kappa_q=linterp(cosbarr,kapd_q,MXTHETB,costb,iinc)
      kap_ext=kap+kappa_q*sqpb/sip
      kapsI=linterp(cosbarr,kaps_i,MXTHETB,costb,iinc)
      kapsQ=linterp(cosbarr,kaps_q,MXTHETB,costb,iinc)
      kaps=kapsI+kapsQ*sqpb/sip
      alb=kaps/kap_ext
c      print*,'alb,kap_ext,Q/I,cost',alb,kap_ext,sqp/sip,cost

c     nokill option
c      sip=sip*alb
c      sqp=sqp*alb
c      sup=sup*alb
c      sqpb=sqpb*alb
c      supb=supb*alb
c      svp=svp*alb

c     kill option
      if(xran.le.alb) then
         
         iscat=iscat+1

c     nokill
c         if (sip.lt.(1.0d-3*alb)) then
c            aflux=aflux+1
c            go to 400
c         endif

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
         
c     shouldn't need to rotate since we already know sqpb,supb,
c     and costb,sintb, etc.
         call rotate(ii(1),ii(2),ii(3),pi,r2p,sqp,sup,sqpb,supb
     1        ,cost,sint,phi,cosp,sinp
     1        ,costb,sintb,phib,cospb,sinpb)
         call stokes16(costb,sintb,cospb,sinpb,phib,sqpb,supb)
         call rotate_back(ii(1),ii(2),ii(3),pi,r2p,sqp,sup,sqpb,supb
     1        ,cost,sint,phi,cosp,sinp
     1        ,costb,sintb,phib,cospb,sinpb)
         pol2=sqpb**2+supb**2+svp**2
         if (pol2 .gt. sip**2) then
            print*,'error, P^2, I^2 ', pol2,sip**2
            continue
         endif

c     kill
      else
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
      ux=sint*cosp
      uy=sint*sinp
      uz=cost
c      print*,''
c      print*,'i,q,u,v before tauint ',sip,sqp,sup,svp
      call tauint(iphot,tsum,ii,r2p,sip,sqp,sup,svp,eps
     1     ,cost,sint,phi,cosp,sinp,pi
     1     ,costb,sintb,phib,cospb,sinpb,sqpb,supb)
c      print*,'i,q,u,v after tauint ',sip,sqp,sup,svp
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
c      k=int(float(nmu-1)*dabs(cost)+1.5d0)
c     NOTE:  no longer assuming z=-z with toroidal B-field
      k=int(float(nmu)*(1.d0-cost)/2.d0)+1      
      if(k.lt.0.or.k.gt.nmu) then
         print*,'cost binning error, k, nmu',k,nmu
         print*, 'cost,sint',cost,sint
         print*, 'iphot',iphot
      endif

      x = ypold*cosp-xpold*sinp
      y = zpold*sini(k)-ypold*u(k)*sinp
     &     -xpold*u(k)*cosp

c     first, sum fluxes
c      rhopold=dsqrt(xpold**2+ypold**2+zpold**2)
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
      
c      if (cost.lt.0.d0) then
c         x = -ypold*cosp+xpold*sinp
c         y = -zpold*sini(k)-ypold*u(k)*sinp
c     &        -xpold*u(k)*cosp
c      else
c         x = ypold*cosp-xpold*sinp
c         y = zpold*sini(k)-ypold*u(k)*sinp
c     &        -xpold*u(k)*cosp
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
