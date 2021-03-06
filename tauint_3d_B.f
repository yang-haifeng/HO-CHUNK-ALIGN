c     **********************************************************

      subroutine tauint(iphot,tsum,ii,r2p,sip,sqp,sup,svp,eps
     1     ,cost,sint,phi,cosp,sinp,pi
     1     ,costb,sintb,phib,cospb,sinpb,sqpb,supb)

      implicit none

      include 'taunum.txt'
      include 'tts.txt'
      include 'opacin.txt'
      include 'grid.txt'
      include 'tab.txt'

      real*8 tsum,dt,opac,dtau,ds,beta,lp1,r2p,eps,deps
     1     ,rsqnew,xpnew,ypnew,zpnew,dtau_eps,dum1,sqp,sup,svp
     1     ,kappa,kap,sqpold,supold,svpold,c1,d1,kap_uv,kappa_q
     1     ,a1,b1,sip,dtaui,dtauq,supb,sqpb,sinpb,cospb,phib
     1     ,sintb,costb,pi,sinp,cosp,phi,sint,cost
      integer ioffset(3),ii(3),iphot,it,ir,ip,ir1,it1,ip1,iinc

      real*8 linterp

      character routine*20
      routine='tauint_3d_B'
      tsum=0.d0

      tsum=0.d0

c     print*,' '
c     print*,'in tauint_3d',sip,sqp,sup,svp

      if (sqp/sip .gt.1.d0.or.sqp/sip.lt.-1.d0) then
         print*,''
         print*,'in TAUINT, before integration, sqp>sip'
         print*,'sqp,sip,sqp/sip',sqp,sip,sqp/sip
      endif

10    continue

      call testgrid(xp,yp,zp,rsq,rtot,ii,routine,iphot,r2p)

      ir=ii(1)
      it=ii(2)
      ip=ii(3)

      if (npg.gt.1.and.ntg.gt.1) then
         call find_wall(dt,ioffset,aarr(ip),aarr(ip+1),barr(ip)
     &        ,barr(ip+1),0.d0,0.d0,0.d0,0.d0,r2arr(ir),r2arr(ir+1)
     &        ,tan2arr(it),tan2arr(it+1),ux,uy,uz,xp,yp,zp)
c     print*,'dt ',dt
c     print*,'ioffsets ',ioffset(1),ioffset(2),ioffset(3)
      else if (ntg.gt.1) then
         call find_wall_2D(dt,ioffset,r2arr(ir),r2arr(ir+1)
     &        ,tan2arr(it),tan2arr(it+1),ux,uy,uz,xp,yp,zp)
      else
         call find_wall_1D(dt,ioffset,r2arr(ir),r2arr(ir+1)
     &        ,ux,uy,uz,xp,yp,zp)
      endif

      if (dt.lt.0.d0) then
         print*,'problem with find_wall! did not find surface! '
         print*,'in tauint,ir,rtot2,r2arr(ir)',ir,rsq,r2arr(ir)
         print*,'ux,uy,uz,iphot',ux,uy,uz,iphot
         print*,'xtot,ytot,ztot',xp,yp,zp
         dt=eps*(rarr(ir+1)-rarr(ir))
      endif

c     rotate stokes to B-field frame, calc. opacity in this frame
c     print*,'tauint loop'
      call rotate(ir,it,ip,pi,r2p,sqp,sup,sqpb,supb
     1     ,cost,sint,phi,cosp,sinp
     1     ,costb,sintb,phib,cospb,sinpb)

      sqpold=sqpb
      supold=supb
      svpold=svp

c     find angle-dependent opacity
      call locate(cosbarr,MXTHETB,nthetab,costb,iinc)
      kap=linterp(cosbarr,kapd,MXTHETB,costb,iinc)
      kappa_q=linterp(cosbarr,kapd_q,MXTHETB,costb,iinc)
      kap_uv=linterp(cosbarr,Ccpol,MXTHETB,costb,iinc)

c     optical depth to edge of cell (add small amount because we
c     will increment distance by a small amount)

      dtaui=(kap-kappa_q)*dt*densarr(ir,it,ip)
      dtauq=(kap+kappa_q)*dt*densarr(ir,it,ip)
      if (dtaui.gt.100.d0) then
         sqpb=sip*0.99999d0
         supb=0.d0
         svp=0.d0
         kappa=kap+kappa_q
         opac=kappa*densarr(ir,it,ip)
         dtau=opac*dt
      else
         a1=(sip-sqpold)*exp(-dtaui)
         b1=(sip+sqpold)*exp(-dtauq)
         if (a1.lt.0.d0) then
            print*,'a < 0',a1
            a1=1.d-3
         endif
         sqpb=(b1-a1)/(a1+b1)*sip
         c1=exp(-kap*dt*densarr(ir,it,ip))
         d1=kap_uv*dt*densarr(ir,it,ip)
         supb=c1*(supold*cos(d1)-svpold*sin(d1))/(a1+b1)*sip*2.d0
         svp=c1*(supold*sin(d1)+svpold*cos(d1))/(a1+b1)*sip*2.d0
c         kappa=kap+kappa_q*0.5*(sqpold+sqp)/sip
         dtau=-dlog(0.5*(a1+b1)/sip)
         opac=dtau/dt
      endif
      deps=(rarr(ir+1)-rarr(ir))*eps
      if ((tsum+dtau).gt.tau) then
         ds=(tau-tsum)/opac
         pathI(ir,it,ip)=pathI(ir,it,ip)+ds
         pathQ(ir,it,ip)=pathQ(ir,it,ip)+sqp*ds
         pathU(ir,it,ip)=pathU(ir,it,ip)+sup*ds
         pathV(ir,it,ip)=pathV(ir,it,ip)+svp*ds
         xp=xp+ds*ux
         yp=yp+ds*uy
         zp=zp+ds*uz
         rsq=xp**2+yp**2+zp**2
         rtot=dsqrt(rsq)
c     recalculate q,u,v
         a1=(sip-sqpold)*exp(-(kap-kappa_q)*ds*densarr(ir,it,ip))
         b1=(sip+sqpold)*exp(-(kap+kappa_q)*ds*densarr(ir,it,ip))
         if (a1.lt.0.d0) then
            print*,'a < 0',a1
            a1=1.d-3
         endif
         sqpb=(b1-a1)/(b1+a1)*sip
         c1=exp(-kap*ds*densarr(ir,it,ip))
         d1=kap_uv*ds*densarr(ir,it,ip)
         supb=c1*(supold*cos(d1)-svpold*sin(d1))/(a1+b1)*sip*2.d0
         svp=c1*(supold*sin(d1)+svpold*cos(d1))/(a1+b1)*sip*2.d0    
c         kappa=kap+kappa_q*0.5*(sqpold+sqp)/sip
         sqpold=sqpb
         supold=supb
         svpold=svp
c        we are done with optical depth integration
         go to 300
      endif

      if (sqpb/sip .gt.1.d0.or.sqpb/sip.lt.-1.d0) then
         print*,'ERROR, in tauint; sqp/sip,iphot',sqpb/sip,iphot
         print*,'sqpold,sqp,sip',sqpold,sqpb,sip
         print*,'a,b'
         print*,a1,b1
         print*,'kap,kappa_q,kap_uv,kappa',kap,kappa_q,kap_uv,kappa
         print*,'dt*densarr(ir,it,ip)',dt*densarr(ir,it,ip)
         print*,'ds*densarr(ir,it,ip)',ds*densarr(ir,it,ip)
         print*,'dt,densarr(ir,it,ip)',dt,densarr(ir,it,ip)
         print*,'rtot,zp,dt',rtot,zp,dt
         print*,'tsum,tau',tsum,tau
         print*,'exponentials'
         print*,exp(-(kap-kappa_q)*dt*densarr(ir,it,ip)),
     1        exp(-(kap+kappa_q)*dt*densarr(ir,it,ip))
         print*,'dtaui,dtauq',dtaui,dtauq
         stop
      endif
c      sqpold=sqp
c      supold=sup
c      svpold=svp

      ii(1)=ii(1)+ioffset(1)
      ii(2)=ii(2)+ioffset(2)
      ii(3)=ii(3)+ioffset(3)
c     account for crossing the phi=0 plane:
      if (ii(3).eq.0) ii(3)=npg-1
      if (ii(3).eq.npg) ii(3)=1

 50   xpnew=xp+(dt+deps)*ux
      ypnew=yp+(dt+deps)*uy
      zpnew=zp+(dt+deps)*uz
      pathI(ir,it,ip)=pathI(ir,it,ip)+dt
      pathQ(ir,it,ip)=pathQ(ir,it,ip)+sqp*dt
      pathU(ir,it,ip)=pathU(ir,it,ip)+sup*dt
      pathV(ir,it,ip)=pathV(ir,it,ip)+svp*dt
      rsqnew=xpnew*xpnew+ypnew*ypnew+zpnew*zpnew
      if ((rsqnew-rsq).eq.0.d0) then
         print*,''
         print*,'ERROR, in tauint_3d'
         print*,'rsqnew=rsq, step size too small'
         print*,'increasing deps'
         deps=deps*2.d0
         go to 50
      endif
      xp=xpnew
      yp=ypnew
      zp=zpnew
      rsq=rsqnew
      rtot=dsqrt(rsq)
      tsum=tsum+dtau

c     testing to see if photon exits
      if (rsq.gt.r2arr(nrg)) then
c        photon exits system
         exitflag=1
         go to 300
      endif

      if (rsq.lt.r2arr(1)) then
c     print*,'rsq,rminsq ',rsq,r2arr(1)
         aflux=aflux+1.d0
         aflag=1
         go to 300
      endif

      call rotate_back(ir,it,ip,pi,r2p,sqp,sup,sqpb,supb
     1     ,cost,sint,phi,cosp,sinp
     1     ,costb,sintb,phib,cospb,sinpb)

      go to 10

300   continue

      call rotate_back(ir,it,ip,pi,r2p,sqp,sup,sqpb,supb
     1     ,cost,sint,phi,cosp,sinp
     1     ,costb,sintb,phib,cospb,sinpb)
c     print*,'sqpb,supb in tauint',sqpb,supb

      return

      end


