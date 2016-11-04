c     **********************************************************

      subroutine diskflux
c
c     calculates flux from disk
c
c history:
c 2000/07/28 (mjw):  change filter on rtmp to .gt.1.d0
c

      implicit none

      include 'tts.txt'
      include 'stokes.txt'
      include 'vger.txt'
      include 'grid.txt'

      real*8 rnorm,fimage,rnp,hstnorm,rtmp,xtmp,ytmp,norm
     1  ,thet

      integer i,ia,ix,iy,it,icenter,l,m,ir,ip

c     total flux from disk
      rnp=dfloat(np)
c      rnorm=rnp/float(nmu-1)/normstar
c      rnorm=rnp/dfloat(nmu-1)
      rnorm=rnp/dfloat(nmu)
c      hstnorm=rnp/normstar
      hstnorm=rnp
      do ia=1,3
      do i=1,nmu
	   rtmp=dfloat(nums(i,ia))
	   if (rtmp.gt.1.d0)  then
c           errors, first, standard deviations
            si2(i,ia)=dsqrt((si2(i,ia)-(si(i,ia)**2)/rtmp)/(rtmp-1.d0))
            sq2(i,ia)=dsqrt((sq2(i,ia)-(sq(i,ia)**2)/rtmp)/(rtmp-1.d0))
            su2(i,ia)=dsqrt((su2(i,ia)-(su(i,ia)**2)/rtmp)/(rtmp-1.d0))
            sv2(i,ia)=dsqrt((sv2(i,ia)-(sv(i,ia)**2)/rtmp)/(rtmp-1.d0))
c           divide by average intensity per photon (assumes error
c           in I is 0)
            rtmp=si(i,ia)/dfloat(nums(i,ia))
	      si2(i,ia)=si2(i,ia)/rtmp
	      sq2(i,ia)=sq2(i,ia)/rtmp
	      su2(i,ia)=su2(i,ia)/rtmp
	      sv2(i,ia)=sv2(i,ia)/rtmp
c           somewhat poisson in that errors should go down with more photons
            rtmp=dsqrt(dfloat(nums(i,ia)))
	      si2(i,ia)=si2(i,ia)/rtmp+1.d0/rtmp
	      sq2(i,ia)=sq2(i,ia)/rtmp
	      su2(i,ia)=su2(i,ia)/rtmp
	      sv2(i,ia)=sv2(i,ia)/rtmp	   
	   endif
c     if(i.eq.1.or.i.eq.nmu) then
c     si(i,ia)=si(i,ia)/rnorm*2.d0
c     sq(i,ia)=sq(i,ia)/rnorm*2.d0
c     su(i,ia)=su(i,ia)/rnorm*2.d0
c     sv(i,ia)=sv(i,ia)/rnorm*2.d0
c     bins 1 and nmu are half the size as the others
c     else
           si(i,ia)=si(i,ia)/rnorm
           sq(i,ia)=sq(i,ia)/rnorm
           su(i,ia)=su(i,ia)/rnorm
           sv(i,ia)=sv(i,ia)/rnorm
c     end if
      end do
	rtmp=dfloat(numt(ia))
	if (rtmp.gt.1.d0) then
c        standard_deviation
         ti2(ia)=dsqrt((ti2(ia)-(ti(ia)**2)/rtmp)/(rtmp-1.d0))
         tq2(ia)=dsqrt((tq2(ia)-(tq(ia)**2)/rtmp)/(rtmp-1.d0))
         tu2(ia)=dsqrt((tu2(ia)-(tu(ia)**2)/rtmp)/(rtmp-1.d0))
         tv2(ia)=dsqrt((tv2(ia)-(tv(ia)**2)/rtmp)/(rtmp-1.d0))
c        divide by average i per photon
         rtmp=ti(ia)/dfloat(numt(ia))
	   ti2(ia)=ti2(ia)/rtmp
	   tq2(ia)=tq2(ia)/rtmp
	   tu2(ia)=tu2(ia)/rtmp
	   tv2(ia)=tv2(ia)/rtmp
c        account for somewhat poisson nature
         rtmp=dsqrt(dfloat(numt(ia)))
	   ti2(ia)=ti2(ia)/rtmp+1.d0/rtmp
	   tq2(ia)=tq2(ia)/rtmp
	   tu2(ia)=tu2(ia)/rtmp
	   tv2(ia)=tv2(ia)/rtmp
      endif

      ti(ia)=ti(ia)/hstnorm
      tq(ia)=tq(ia)/hstnorm
      tu(ia)=tu(ia)/hstnorm
      tv(ia)=tv(ia)/hstnorm
      tis(ia)=tis(ia)/hstnorm
      tqs(ia)=tqs(ia)/hstnorm
      tus(ia)=tus(ia)/hstnorm
      tvs(ia)=tvs(ia)/hstnorm
      end do
      tid=tid/hstnorm
      tqd=tqd/hstnorm
      tud=tud/hstnorm
      tvd=tvd/hstnorm

c     check constants
c     rnorm has factor of 2*pi because you made all the photons
c     come out at one phi.  
c     rnorm=r2p*delnt*float(np)
c     that was the value for the correct rnorm.  using that
c     value, the flux at a given angle from the star, if there
c     is not absorbing disk, is 1/4pi.  We will multiply the
c     flux by 4pi so that it is normalized to the stellar flux.
c      rnorm=0.5*dmu*float(np)/normstar
      rnorm=0.5d0*dmu*dfloat(np)

c     multiply rnorm by 2 if you are doubling photons by
c     symmetry through z axis (in diskim)
      rnorm=rnorm*1.d0
      fimage=0.d0
      print*,'normalizing flux image (not pol images) by ',rnorm
c     whenever we deal with bins 1 and nmu, multiply by 2 because
c     bin is half size of others.

      icenter=nx/2+1
      do it=1,nmu
c     if(it.eq.1.or.it.eq.nmu) then
c     star(it)=star(it)/(rnorm)*2.d0
c     starq(it)=starq(it)/(rnorm)*2.d0
c     staru(it)=staru(it)/(rnorm)*2.d0
c     starv(it)=starv(it)/(rnorm)*2.d0
c     else
         star(it)=star(it)/(rnorm)
         starq(it)=starq(it)/(rnorm)
         staru(it)=staru(it)/(rnorm)
         starv(it)=starv(it)/(rnorm)
c     end if
      end do
      do it=1,nmu
      do ix=1,nx
      do iy=1,nx
c        standard deviation of imagei array
	   rtmp=dfloat(numi(ix,iy,it))
	   if (image2(ix,iy,it).gt.1.d0) then
c           print*,'rtmp',rtmp
            image2(ix,iy,it)=(image2(ix,iy,it)-imagei(ix,iy,it)**2/rtmp)
     1      /rtmp/(rtmp-1.d0)
	   endif
         if(ix.eq.icenter.and.iy.eq.icenter) then
            write(6,*) 'image before norm'
     1      ,image(icenter,icenter,it),icenter,it
         end if
         if(image(ix,iy,it).lt.0.d0) then
            write(6,*) 'before normalizing'
            write(6,*) 'image lt 0',image(ix,iy,it),ix,iy,it
         end if
c         if(it.eq.1.or.it.eq.nmu) then
c            image(ix,iy,it)=image(ix,iy,it)/(rnorm)*2.d0
c         else
            image(ix,iy,it)=image(ix,iy,it)/(rnorm)
c         end if
         fimage=fimage+(image(ix,iy,it))/dfloat(nmu)
         if(ix.eq.icenter.and.iy.eq.icenter) then
            write(6,*) 'image(icenter,icenter,it),icenter,it'
     1     ,image(icenter,icenter,it),icenter,it
         end if
         if(image(ix,iy,it).lt.0.d0) then
            write(6,*) 'image lt 0',image(ix,iy,it),ix,iy,it
         end if
      end do
      end do
      end do

      do ix=1,nxhst
      do iy=1,nxhst
         tihst(ix,iy)=tihst(ix,iy)/(hstnorm)
         tqhst(ix,iy)=tqhst(ix,iy)/(hstnorm)
         tuhst(ix,iy)=tuhst(ix,iy)/(hstnorm)
         tvhst(ix,iy)=tvhst(ix,iy)/(hstnorm)
      enddo
      enddo

      norm=pi/bincvg*r2p/binpvg*dfloat(np)
      do l=1,ncvg
         thet=pi/bincvg*(dfloat(l)-0.5)
c     print*,'thet',thet
         do m=1,npvg
            ivg(l,m)=ivg(l,m)/norm/dsin(thet)
            qvg(l,m)=qvg(l,m)/norm/dsin(thet)
            uvg(l,m)=uvg(l,m)/norm/dsin(thet)
            vvg(l,m)=vvg(l,m)/norm/dsin(thet)
         enddo
      enddo
c     errors;
      if (countvg.gt.0) then
         rtmp=countvg
c        first, standard deviations
         ivgerr=dsqrt((ivgerr-ivgflux**2/rtmp)/(rtmp-1.d0))
         qvgerr=dsqrt((qvgerr-qvgflux**2/rtmp)/(rtmp-1.d0))
         uvgerr=dsqrt((uvgerr-uvgflux**2/rtmp)/(rtmp-1.d0))
         vvgerr=dsqrt((vvgerr-vvgflux**2/rtmp)/(rtmp-1.d0))
c        we've now calculated average q,u,v st.dev. per collected photon
c        want error in q/i, u/i, v/i
c        assume i has no error.  divide by average i per photon
         rtmp=ivgflux/countvg
	   ivgerr=ivgerr/rtmp
	   qvgerr=qvgerr/rtmp
	   uvgerr=uvgerr/rtmp
	   vvgerr=vvgerr/rtmp
c        now, it's not exactly poisson since each photon has a different
c        intensity, but the errors do decrease with increasing number
c        of photons.  so you want to decrease the errors by some
c        sqrt(1./n) type of number.  instead of using total number
c        of photons, we might vary it, make it smaller
         rtmp=dsqrt(1.d0*dfloat(countvg))
	   ivgerr=ivgerr/rtmp+1./rtmp   !in case st.dev.=0.
	   qvgerr=qvgerr/rtmp
	   uvgerr=uvgerr/rtmp
	   vvgerr=vvgerr/rtmp
      endif
c     normalize 
	ivgflux=ivgflux/rnp
	qvgflux=qvgflux/rnp
	uvgflux=uvgflux/rnp
	vvgflux=vvgflux/rnp

      write(6,*)'flux from image',fimage

c     calculate average intensity in grid for four stokes vectors
      do ir=1,nrg-1
      do it=1,ntg-1
      do ip=1,npg-1
c     ds is in units of rstar so multiply by rmincgs, right?
c     will have to multiply output by stellar luminosity at this
c     wavelength, if it is known
         pathI(ir,it,ip)=pathI(ir,it,ip)/r2p*volarr(ir,it,ip)
     1      *rstar*rsol
         pathQ(ir,it,ip)=pathQ(ir,it,ip)/r2p*volarr(ir,it,ip)
     1      *rstar*rsol
         pathU(ir,it,ip)=pathU(ir,it,ip)/r2p*volarr(ir,it,ip)
     1      *rstar*rsol
         pathV(ir,it,ip)=pathV(ir,it,ip)/r2p*volarr(ir,it,ip)
     1      *rstar*rsol
      enddo
      enddo
      enddo

      return
      end



c     *********************************************************
