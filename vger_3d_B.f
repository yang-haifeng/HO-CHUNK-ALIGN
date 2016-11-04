      subroutine vger(dum1,dum2,dum3,dum4,dum5,dum6,dum7
     1   ,dum9,dum10,dum11,dum12,dum13
     1   ,pi,r2p,hit,htot,dum15,dum16
     1   ,peak,tsum,dum14,dum8,iflag,iphot,eps)

c     don't want this routine to change any of the input variables!!!!
c     fortran 90 would make this easier....

c     WARNING, HAVE NOT IMPLEMENTED CHANGES FOR ARBRITRARY B-FIELD
C     AND ALIGNED GRAINS.

      implicit none

      real*8 s,dsmax,xp1,yp1,zp1,rp1,ux1,uy1,uz1
     1   ,pi,r2p,sip,sqp,sup,svp,rsq1
     $   ,bmu,sinbm
     $   ,costp,sintp,cospp,sinpp,phip
     $   ,taumax,dtau,eps,deps
      real*8 tsum
     1   ,xpnew,ypnew,zpnew,opac,dscost,rsqnew
     1   ,dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum9,dum10
     1   ,dum11,dum12,dum13,dum3a,stot2,stot,phi
     1   ,phinew,ri1,hit,htot,peak,xold,yold,zold
     1   ,rhopold,factor,deltax,deltay,deltaz,sint,cost,sinp,cosp
     1   ,phiin,sinpin,costin,cospin,sintin,thetin,rnum,kappa,dtau_eps
     1   ,dt,lp1,beta,newdt,newdtau,dum15,dum16,rtot1,sipold,coskapd
     1   ,kap,kappa_q,kap_uv,sqpold,supold,svpold,c1,d1,a1,b1


      integer dum14(3),ii(3),ioffset(3),l,m,iinc
      integer idust,dum8,iphot,iflag,ir,it,ip

      real*8 check
      character*20 routine
      real*8 linterp

      include 'opacin.txt'
      include 'tts.txt'
      include 'vger.txt'
      include 'grid.txt'
      include 'tab.txt'

      routine='vger_3d_di'
      
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

      tsum=0.d0
      s=0.d0
      taumax=25.d0

      deltax=xvg-xold
      deltay=yvg-yold
      deltaz=zvg-zold
      stot2=deltax**2+deltay**2+deltaz**2
      stot=dsqrt(stot2)
      cost=deltaz/stot
      cost=check(cost,routine)
      sint=dsqrt(1.-cost**2) !theta ranges from 0-pi so sint>0
      phi=datan2(deltay,deltax)
      if (phi.lt.0) phi=phi+r2p
      cosp=dcos(phi)
      sinp=dsin(phi)
c     rsq1=xp1*xp1+yp1*yp1+zp1*zp1
c     rhopold=dsqrt(rsq1) 
      rhopold=rtot1
      ux1=sint*cosp
      uy1=sint*sinp
      uz1=cost

c     need new value of sqp for opacity calculation so call
c     stokespeel here instead of down below like in spherical dust
c     codes
      if (iflag.eq.1) then
         call stokespeel16(sip,sqp,sup,svp,costp,cost
     1        ,hit,htot,pi,r2p,peak,phip,phi)
      endif
      sipold=sip
      sqpold=sqp
      supold=sup
      svpold=svp

c     find direction-dependent opacity
      coskapd=(costp)
      call locate(cosbarr,MXTHETB,nthetab,coskapd,iinc)
      kap=linterp(cosbarr,kapd,MXTHETB,coskapd,iinc)
      kappa_q=linterp(cosbarr,kapd_q,MXTHETB,coskapd,iinc)
      kap_uv=linterp(cosbarr,Ccpol,MXTHETB,coskapd,iinc)
c     kappa=kap+kappa_q*sqp/sip
c     print*,'in peeloff_3d_di,kap,kappa,coste',kap,kappa,coste


      do while (s.lt.stot)
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
         dt=eps*(rarr(ir+1)-rarr(ir))
      endif

      a1=(sipold-sqpold)*exp(-(kap-kappa_q)*dt*densarr(ir,it,ip))
      b1=(sipold+sqpold)*exp(-(kap+kappa_q)*dt*densarr(ir,it,ip))
      if (a1.lt.0.d0) then
         print*,'a < 0',a1
         a1=1.d-3
      endif
c     play
      sip=0.5*(a1+b1)
      sqp=0.5*(b1-a1)
c     sqp=(b-a)/(a+b)*sip
      c1=exp(-kap*dt*densarr(ir,it,ip))
      d1=kap_uv*dt*densarr(ir,it,ip)
c      sup=c*(supold*cos(d)-svpold*sin(d))/(a+b)*sip
c      svp=c*(supold*sin(d)+svpold*cos(d))/(a+b)*sip
      sup=c1*(supold*cos(d1)-svpold*sin(d1))
      svp=c1*(supold*sin(d1)+svpold*cos(d1))
c      print*,'sip,sqpold,sqp',sip,sqpold,sqp
      kappa=kap+kappa_q*0.5*(sqpold/sipold+sqp/sip)
      if (sqp/sip .gt.1.d0.or.sqp/sip.lt.-1.d0) then
         print*,'ERROR, in peeloff_3d; sqp/sip',sqp/sip
         print*,'sqpold,sqp,sip',sqpold,sqp,sip
         print*,'a,b'
         print*,a,b
         print*,'exponentials'
         print*,exp(-(kap-kappa_q)*dt*densarr(ir,it,ip)),
     1        exp(-(kap+kappa_q)*dt*densarr(ir,it,ip))
      endif  
      opac=kappa*densarr(ir,it,ip)
      deps=(rarr(ir+1)-rarr(ir))*eps
      dtau_eps=opac*deps
      if (dtau_eps.gt.0.001) then
         print*,''
         print*,'WARNING, in peeloff_3d'
         print*,'deps is too big, opac*deps',dtau_eps
         print*,''
      endif
      dtau=opac*dt+dtau_eps  
      tsum=tsum+dtau
      if ((tsum).gt.taumax) then
c        photon weight will be too small to matter
         go to 300
      endif
      sipold=sip
      sqpold=sqp
      supold=sup
      svpold=svp

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
      s=s+dt

c     testing
      if (rsq1.gt.r2arr(nrg)) then
c        photon exits system
         go to 250
      endif

      if (rsq1.lt.r2arr(1)) then
         go to 300
      endif

      if (s.ge.stot) then
         go to 250
      endif

      enddo
 
250   continue
c     write(6,*) 'tsum.gt.tau',tsum,tau
      
c     if you are here, you reached vger
c     but probably went a little beyond stot. backup to stot.
c     calculate newdt = distance from last wall to pathfinder
      newdt=s-stot
c     check, newdt should be smaller than dt
      if (newdt.gt.dt) then
         print*,'shooot, error in vger'
         print*,'newdt should be < dt',newdt,dt
c     newdtau is optical depth corresponding to that distance
         newdtau=newdt*opac
c     correct tsum
         tsum=tsum-dtau+newdtau
c     correct stokes parameters
         a1=(sip-sqpold)*exp(-(kap-kappa_q)*newdt*densarr(ir,it,ip))
         b1=(sip+sqpold)*exp(-(kap+kappa_q)*newdt*densarr(ir,it,ip))
         c1=exp(-kap*newdt*densarr(ir,it,ip))
         d1=kap_uv*newdt*densarr(ir,it,ip)
         if (a1.lt.0.d0) then
            print*,'a < 0',a1
            a1=1.d-3
         endif
         sip=0.5*(a1+b1)
         sqp=0.5*(b1-a1)
         sup=c1*(supold*cos(d1)-svpold*sin(d1))
         svp=c1*(supold*sin(d1)+svpold*cos(d1))
         kappa=kap+kappa_q*0.5*(sqpold+sqp)/sip
c         sqpold=sqp
c         supold=sup
c         svpold=svp
      endif

c     weight photon intensity by optical depth and
c     account for probability of intercepting solid angle of 
c     scattered light.
c     vgfactor is rstar^2 in cgs.  because, normalizing to stellar
c     flux=N/(4*pi*rstar^2).  and also dividing by 4*pi, so they
c     cancel out.
c      if (iflag.eq.1) then
c         factor=dexp(-tsum)*vgfactor
c         sip=sip*factor
c         sqp=sqp*factor
c         sup=sup*factor
c         svp=svp*factor
cc        if (sip.lt.1.e-10) go to 300
c	 call findangle3(pi,r2p,costp,sintp,phip,cospp,sinpp
c     1   ,cost,sint,phi,cosp,sinp,bmu,sinbm,phinew,ri1)
cc        don't care about phinew
c         call stokespeel6(sip,sqp,sup,svp,cost,sint
c     1   ,hit,htot,pi,r2p,peak,ri1,bmu,costp,sintp)
c      else
c         factor=dexp(-tsum)*vgfactor
c         sip=sip*factor
c         sqp=sqp*factor
c         sup=sup*factor
c         svp=svp*factor
c      endif

      phiin=phi+pi               !reverse direction to get incoming
      if (phiin .gt. r2p) phiin = phiin - r2p
      cospin=-cosp
      sinpin=-sinp
      costin=-cost               !reverse direction to get incoming
      sintin=sint
      thetin=acos(costin)
      countvg=countvg+1
c     vger binning
      rnum=binpvg*phiin/r2p
      m=int(rnum+0.5)+1
      if (m.eq.(int(binpvg)+1)) m=1
                                   !this is the phi=0 bin
                                   !     binning in theta
      rnum=bincvg*thetin/pi
      l=int(rnum)+1
      if (m.gt.0.and.m.le.npvg.and.l.gt.0.and.l.le.ncvg) then
         ivg(l,m)=ivg(l,m)+sip
         qvg(l,m)=qvg(l,m)+sqp
         uvg(l,m)=uvg(l,m)+sup
         vvg(l,m)=vvg(l,m)+svp
         ivgflux=ivgflux+sip
         qvgflux=qvgflux+sqp
         uvgflux=uvgflux+sup
         vvgflux=vvgflux+svp
         ivgerr=ivgerr+sip*sip
         qvgerr=qvgerr+sqp*sqp
         uvgerr=uvgerr+sup*sup
         vvgerr=vvgerr+svp*svp
      else
        print*,'array out of bounds, pathfinder '
        print*,'m,l,npvg,ncvg ',m,l,npvg,ncvg
      endif


 300  continue
      return
      end
      
