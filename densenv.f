c     ***********************************************************
      subroutine densenv(rad,thetin,costin,sina,pi,pihalf,phi
     1     ,dens)
      
c     calculates density in envelope. called during grid setup
c     so doesn't need to be optimized.

c  history:
c  00/09/06 (baw): write subroutine, modified from opacenv.f
c  01/02/17 (baw): combine opacitybub subroutine into here.

      implicit none
      include 'tts.txt'
      include 'opacin.txt' 

      real*8 rad,dens,cosa,rado,radi,phi,thet,thetin,pihalf,y
     $     ,pi,r2,costin,m,n,fact,sina,xmu,rx,rp,zup,xmu0,factor
     $     ,rx2,xmu0new,zlo,zp

      integer iflag

c ... begin g77
      real*8 dcosd,dsind,dtand
      external dcosd,dsind,dtand
c     ... end g77


c     for now, assuming +z = -z
      thet=thetin
      cosa=abs(costin)
      if (thet.gt.pihalf) thet=pi-thet
      rp=rad*sina
      zp=rad*cosa

      zup=c1e*rp**ex1+z01
      xmu=zp/rad
      if(rad.gt.rmine) then
         rx=rad/rc
         call zerod(xmu,rx,xmu0,iflag)
c         if (iflag.eq.1) then
c            print*,'thetin,pi',thetin*180.d0/pi,pi
c         endif
c     if (iflag.eq.1) print*,'xmu,rx,xmu0,iflag',xmu,rx,xmu0,iflag
         factor=1.d0/((xmu/xmu0)+2.d0*(xmu0**2)/rx)
c        dens=rhoe0*(rx**(-1.5))*factor/(dsqrt(1.d0+xmu/xmu0))
	 dens=1.d-18 ! set to constant envelope for now.
         rx2=rad/rchole         
         call zerod(xmu,rx2,xmu0new,iflag)
c     if (iflag.eq.1) print*,'xmu,rx2,xmu0new,iflag',xmu,rx2,
c     1        xmu0new,iflag
      else
         dens=0.
      end if         

      if(ihole.eq.1) then
         if(istream.eq.1.and.xmu0new.gt.windmu0) then
            if(zp.gt.zup) then
               dens=rhoconst1*rad**(-exf)
            else
               dens=rhoconst2*rad**(-exf)
            end if
            if(rad.lt.zflowmin) dens=rhoamb
         end if
         if(ipoly.eq.1) then
            zlo=c2e*rp**ex2+z02
c     baw: 2/5/99 new.  doing this so holes can intersect and create 
c     different shapes!  
            if (zp.gt.zlo) then
               dens=rhoconst2*rad**(-exf)
               if(rad.lt.zflowmin) dens=rhoamb
            endif
            if(zp.gt.zup) then
               dens=rhoconst1*rad**(-exf)
               if(rad.lt.zflowmin) dens=rhoamb
            end if
         end if
         if(ibub.eq.1) then
            if (rad.lt.zbub2*xmu**nbub) then
               if (rad.lt.zbub1*xmu**nbub) then
                  dens=rhoconst1
               else
                  dens=rhoconst2
               endif
            endif
c     if(rp.lt.roa.and.rad.gt..8d0*zbub1) dens=rhoconst1
            if(xmu.gt.cosbuboa) dens=rhoconst1
         endif
      end if
c     make density stop at r>rmax
      if (rad. gt. rmax) dens=0.d0
c     compare to ambient density; if rho lt rhoamb, set rho to rhoamb.
      if(dens.lt.rhoamb) then
         dens=rhoamb
      endif
      
c     calculate vphi---velocity in azimuthal direction.
c     xsint=rp/rad
c     xsint0=sqrt(1.-xmu0**2)
c     308.72=sqrt(G*M/rstar/rsol) so rad is in correct units
c     vphi=308.72/sqrt(rad)*sqrt(1.-xmu/xmu0)*xsint0/xsint

      return 
      end

