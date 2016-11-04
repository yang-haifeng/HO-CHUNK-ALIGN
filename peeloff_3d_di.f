      subroutine peeloff(dum1,dum2,dum3,dum4,dum5,dum6,dum7
     1   ,dum9,dum10,dum11,dum12,dum13
     1   ,pi,r2p,hit,htot,dum15,dum16
     1   ,peak,tsum,dum14,dum8,iflag,iphot,eps,kaps)

c     don't want this routine to change any of the input variables!!!!
c     fortran 90 would make this easier....

c     integrate optical depth until photon exits in a specified
c     direction, or hits the star.

      implicit none

      real*8 s,dsmax,xp1,yp1,zp1,ux1,uy1,uz1
     1   ,pi,r2p,sip,sqp,sup,svp,rsq1
     $   ,bmu,sinbm,yimage,ximage
     $   ,costp,sintp,cospp,sinpp,phip
     $   ,taumax,dtau,eps,deps,dtau_eps
      real*8 tsum,rtau,kap,kappa_q,kappa
     1   ,xpnew,ypnew,zpnew,opac,dscoste,rsqnew
     1   ,dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum9,dum10
     1   ,dum11,dum12,dum13,dum3a,dum15,dum16,rtot1
     1   ,phinew,ri1,hit,htot,peak,xold,yold,zold
     1   ,rhopold,factor,dt,lp1,beta,coskapd,dtaui,dtauq
     1     ,sipold,sqpold,supold,svpold,c1,d1,kap_uv,a1,b1
     1     ,kapsI,kapsQ,norm,kaps

      integer dum14(3),ii(3),ioffset(3)
      integer idust,ix,iy,irho,dum8,iphot,iflag,ir,it,ip
     1 ,iinc

      character routine*20

      real*8 linterp

      include 'opacin.txt'
      include 'tts.txt'
      include 'grid.txt'
      include 'tab.txt'

      routine='peeloff_3d_di'

c     copy input variables so as not to change them
      xp1=dum1
      yp1=dum2
      zp1=dum3
      xold=xp1
      yold=yp1
      zold=zp1
      sip=dum4
      sqp=dum5
      sup=dum6
      svp=dum7
      idust=dum8
      costp=dum9
      sintp=dum10
      cospp=dum11
      sinpp=dum12
      phip=dum13
      ii(1)=dum14(1)
      ii(2)=dum14(2)
      ii(3)=dum14(3)
      rsq1=dum15
      rtot1=dum16

c     print*,'input variables'
c     print*,xp1,yp1,zp1,sip,sqp,sup,svp,idust,costp,sintp,cospp,
c     1     sinpp,phip,ii(1),ii(2),ii(3)

      tsum=0.d0
      taumax=25.d0

c     need new value of sqp for opacity calculation so call
c     stokespeel here instead of down below like in spherical dust
c     codes
      sipold=sip
      sqpold=sqp
      supold=sup
      svpold=svp
      if (iflag.eq.1) then
         call stokespeel16(sip,sqp,sup,svp,costp,coste
     1        ,hit,htot,pi,r2p,peak,phip,phie)
c         coskapd=costp
c         call locate(cosbarr,MXTHETB,nthetab,coskapd,iinc)
c         kapsI=linterp(cosbarr,kaps_i,MXTHETB,coskapd,iinc)
c         kapsQ=linterp(cosbarr,kaps_q,MXTHETB,coskapd,iinc)
c         norm=(kapsI*sipold+kapsQ*sqpold)/sipold
         sipold=sip/kaps*5.19209325d16
         sqpold=sqp/kaps*5.19209325d16
         supold=sup/kaps*5.19209325d16
         svpold=svp/kaps*5.19209325d16
c     factor is (XH/mH/1.d8)*4*pi  to get it back to integrated phase function
      endif

c     find direction-dependent opacity
      coskapd=(coste)
      call locate(cosbarr,MXTHETB,nthetab,coskapd,iinc)
      kap=linterp(cosbarr,kapd,MXTHETB,coskapd,iinc)
      kappa_q=linterp(cosbarr,kapd_q,MXTHETB,coskapd,iinc)
      kap_uv=linterp(cosbarr,Ccpol,MXTHETB,coskapd,iinc)
c     kappa=kap+kappa_q*sqp/sip
c     print*,'in peeloff_3d_di,kap,kappa,coste',kap,kappa,coste

      rhopold=rtot1
      if (rhopold.lt.rarr(1)) then
         print*,'wo!  peeloff, error, r ',rhopold
         print*,'rarr(1)',rarr(1)
      endif
      ux1=sinte*cospe
      uy1=sinte*sinpe
      uz1=coste

c     print*,'ux1,uy1,uz1',ux1,uy1,uz1
c     print*,'rsq',rsq1,xp1*xp1+yp1*yp1+zp1*zp1
c     print*,'rtot',rtot1,dsqrt(xp1*xp1+yp1*yp1+zp1*zp1)

c      print*,''
c      print*,'before loop,sip,sqpold,sqp',sip,sqpold,sqp
      if (sqp/sip .gt.1.d0.or.sqp/sip.lt.-1.d0) then
         print*,''
         print*,'in peeloff, before integration, sqp>sip'
         print*,'sqp,sip,sqp/sip',sqp,sip,sqp/sip
      endif

10    continue
      call testgrid(xp1,yp1,zp1,rsq1,rtot1,ii,routine,iphot,r2p)
      ir=ii(1)
      it=ii(2)
      ip=ii(3)

      if (npg.gt.1.and.ntg.gt.1) then
         call find_wall(dt,ioffset,aarr(ip),aarr(ip+1),barr(ip)
     &        ,barr(ip+1),0.d0,0.d0,0.d0,0.d0,r2arr(ir),r2arr(ir+1)
     &        ,tan2arr(it),tan2arr(it+1),ux1,uy1,uz1,xp1,yp1,zp1)
c     print*,'dt ',dt
c     print*,'ioffsets ',ioffset(1),ioffset(2),ioffset(3)
      else if (ntg.gt.1) then
         call find_wall_2D(dt,ioffset,r2arr(ir),r2arr(ir+1)
     &        ,tan2arr(it),tan2arr(it+1),ux1,uy1,uz1,xp1,yp1,zp1)
      else
         call find_wall_1D(dt,ioffset,r2arr(ir),r2arr(ir+1)
     &        ,ux1,uy1,uz1,xp1,yp1,zp1)
      endif

      if (dt.lt.0.d0) then
         print*,'problem with find_wall! did not find surface! '
         print*,'in peeloff,ir,rtot2,r2arr(ir)',ir,rsq1,r2arr(ir)
         print*,'ux,uy,uz,iphot',ux1,uy1,uz1,iphot
         print*,'xtot,ytot,ztot',xp1,yp1,zp1
         dt=eps*(rarr(ir+1)-rarr(ir))
      endif

      dtaui=(kap-kappa_q)*dt*densarr(ir,it,ip)
      dtauq=(kap+kappa_q)*dt*densarr(ir,it,ip)
      if (dtaui.gt.100.d0) then
         sip=0.d0
         sqp=sip*0.999999d0
         sup=0.d0
         svp=0.d0
         kappa=kap+kappa_q
      else
         a1=(sipold-sqpold)*exp(-dtaui)
         b1=(sipold+sqpold)*exp(-dtauq)
         if (a1.lt.0.d0) then
            print*,'a1 < 0',a1
            a1=1.d-3
         endif
c     play
         sip=0.5*(a1+b1)
         sqp=0.5*(b1-a1)
         c1=exp(-kap*dt*densarr(ir,it,ip))
         d1=kap_uv*dt*densarr(ir,it,ip)
         sup=c1*(supold*cos(d1)-svpold*sin(d1))
         svp=c1*(supold*sin(d1)+svpold*cos(d1))
c     print*,'sip,sqpold,sqp',sip,sqpold,sqp
         kappa=kap+kappa_q*0.5*(sqpold/sipold+sqp/sip)
      endif
      if (sip.gt.0.d0) then
      if (sqp/sip .gt.1.d0.or.sqp/sip.lt.-1.d0) then
         print*,'ERROR, in peeloff_3d; sqp/sip',sqp/sip
         print*,'sqpold,sqp,sip',sqpold,sqp,sip
         print*,'a,b'
         print*,a1,b1
         print*,'kap,kappa_q,kap_uv,kappa',kap,kappa_q,kap_uv,kappa
         print*,'dt*densarr(ir,it,ip)',dt*densarr(ir,it,ip)
         print*,'dt,densarr(ir,it,ip)',dt,densarr(ir,it,ip)
         print*,'rtot,zp,dt',rtot1,zp1,dt
         print*,'tsum',tsum
         print*,'exponentials'
         print*,exp(-(kap-kappa_q)*dt*densarr(ir,it,ip)),
     1        exp(-(kap+kappa_q)*dt*densarr(ir,it,ip))
         print*,'dtaui,dtauq',dtaui,dtauq
         stop
      endif  
      endif
      sipold=sip
      sqpold=sqp
      supold=sup
      svpold=svp
      opac=kappa*densarr(ir,it,ip)
      deps=(rarr(ir+1)-rarr(ir))*eps
      dtau_eps=opac*deps
      if (dtau_eps.gt.0.001) then
         print*,''
         print*,'WARNING, in peeloff_3d'
         print*,'deps is too big, opac*deps',dtau_eps
         print*,''
      endif
      dtau=opac*dt+dtau_eps      !optical depth to edge of cell
      tsum=tsum+dtau
      if ((tsum).gt.taumax) then
c        photon weight will be too small to matter
         go to 300
      endif

      ii(1)=ii(1)+ioffset(1)
      ii(2)=ii(2)+ioffset(2)
      ii(3)=ii(3)+ioffset(3)
c     account for crossing the phi=0 plane:
      if (ii(3).eq.0) ii(3)=npg-1
      if (ii(3).eq.npg) ii(3)=1

 50   xpnew=xp1+(dt+deps)*ux1
      ypnew=yp1+(dt+deps)*uy1
      zpnew=zp1+(dt+deps)*uz1
      rsqnew=xpnew*xpnew+ypnew*ypnew+zpnew*zpnew
      if ((rsqnew-rsq1).eq.0.d0) then
         print*,''
         print*,'ERROR, in peeloff_3d'
         print*,'rsqnew=rsq, step size too small'
         print*,'increasing deps'
         deps=deps*2.d0
         go to 50
      endif
      xp1=xpnew
      yp1=ypnew
      zp1=zpnew
      rsq1=rsqnew
      rtot1=dsqrt(rsq1)

c     test to see if photon exits system
      if (rsq1.gt.r2arr(nrg)) then
c        photon exits system
         go to 250
      endif

      if (rsq1.lt.r2arr(1)) then
         go to 300
      endif

      go to 10
250   continue
c     write(6,*) 'tsum.gt.tau',tsum,tau

      
c     if you are here, you left the atmosphere and are heading to HST
c     weight photon intensity by optical depth and
c     account for probability of intercepting solid angle of scattered light.
c     sfactor is 1/4pi
c      if (iflag.eq.1) then
c         factor=exp(-tsum)*sfactor
c         factor=dexp(-tsum)
c         sip=sip*factor
c         sqp=sqp*factor
c         sup=sup*factor
c         svp=svp*factor
c         if (sip.lt.1.e-10) go to 300
c         call findangle3(pi,r2p,costp,sintp,phip,cospp,sinpp
c     1   ,coste,sinte,phie,cospe,sinpe,bmu,sinbm,phinew,ri1)
c     don't care about phinew
c     if (idust.lt.2) then
c     call stokespeel(idust,sip,sqp,sup,svp,coste,sinte
c     1      ,hit,htot,pi,r2p,peak,ri1,bmu,costp,sintp)
c     else
c         call stokespeel6(sip,sqp,sup,svp,coste,sinte
c     1        ,hit,htot,pi,r2p,peak,ri1,bmu,costp,sintp)
c     endif
c      else
c         factor=dexp(-tsum)
c         sip=sip*factor
c         sqp=sqp*factor
c         sup=sup*factor
c         svp=svp*factor
c      endif

      yimage = zold*sinte-yold*coste*sinpe
     &   -xold*coste*cospe
      ximage = yold*cospe-xold*sinpe
c     print*,'ximage,rmaxi,fractxh ',ximage,rmaxi,fractxh
      ix = int((ximage+rmaxi)/fractxh)+1
      iy = int((yimage+rmaxi)/fractxh)+1
      if (ix.le.nxhst.and.iy.le.nxhst.and.ix.gt.0.and.iy.gt.0) then
         tihst(ix,iy) = tihst(ix,iy) + (sip)
         tqhst(ix,iy) = tqhst(ix,iy) + (sqp)
         tuhst(ix,iy) = tuhst(ix,iy) + (sup)
         tvhst(ix,iy) = tvhst(ix,iy) + (svp)
      endif

      if (rhopold.lt.aperture(1)) then
         ti(1)=ti(1)+(sip)
         tq(1)=tq(1)+(sqp)
         tu(1)=tu(1)+(sup)
         tv(1)=tv(1)+(svp)
         numt(1)=numt(1)+1
         ti2(1)=ti2(1)+sip*sip
         tq2(1)=tq2(1)+sqp*sqp
         tu2(1)=tu2(1)+sup*sup
         tv2(1)=tv2(1)+svp*svp
         if (iflag.eq.1) then
            tis(1)=tis(1)+(sip)
            tqs(1)=tqs(1)+(sqp)
            tus(1)=tus(1)+(sup)
            tvs(1)=tvs(1)+(svp)
         endif
      endif
      if (rhopold.lt.aperture(2)) then
         ti(2)=ti(2)+(sip)
         tq(2)=tq(2)+(sqp)
         tu(2)=tu(2)+(sup)
         tv(2)=tv(2)+(svp)
         numt(2)=numt(2)+1
         ti2(2)=ti2(2)+sip*sip
         tq2(2)=tq2(2)+sqp*sqp
         tu2(2)=tu2(2)+sup*sup
         tv2(2)=tv2(2)+svp*svp
         if (iflag.eq.1) then
            tis(2)=tis(2)+(sip)
            tqs(2)=tqs(2)+(sqp)
            tus(2)=tus(2)+(sup)
            tvs(2)=tvs(2)+(svp)
         endif
      endif
      if (rhopold.lt.aperture(3)) then
         ti(3)=ti(3)+(sip)
         tq(3)=tq(3)+(sqp)
         tu(3)=tu(3)+(sup)
         tv(3)=tv(3)+(svp)
         numt(3)=numt(3)+1
         ti2(3)=ti2(3)+sip*sip
         tq2(3)=tq2(3)+sqp*sqp
         tu2(3)=tu2(3)+sup*sup
         tv2(3)=tv2(3)+svp*svp
         if (iflag.eq.1) then
            tis(3)=tis(3)+(sip)
            tqs(3)=tqs(3)+(sqp)
            tus(3)=tus(3)+(sup)
            tvs(3)=tvs(3)+(svp)
         endif
      endif

      if (iflag.eq.0) then
         tid=tid+sip
         tqd=tqd+sqp
         tud=tud+sup
         tvd=tvd+svp
      endif

 300  continue
      return
      end
      
