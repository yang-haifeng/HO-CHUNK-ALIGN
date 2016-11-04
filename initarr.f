      subroutine initarr

      implicit none

      include 'tts.txt'
      include 'vger.txt'
      include 'grid.txt'

      integer i,j,k,ir,it,ip
	
c     set up angle array
c      dmu=1.d0/(float(nmu)-1.d0)
c      dmuhalf=.5d0*dmu
c      do i=1,nmu
c         u(i)=(float(i)-1.d0)*dmu
c         write(6,*) 'mu',u(i)
c      end do

c     theta goes from 0-pi
      dmu=2.d0/(float(nmu))
      do i=1,nmu
         u(i)=1.d0-(dfloat(i)-0.5d0)*dmu
c     th(i)=dacos(u(i)
c     sini(i)=dsqrt(1.-u(i)**2)
c     print*,'u(i) ',u(i)
         print*,'mu(i) ',u(i)
      end do

c     initialize arrays
      do k=1,nmu
         star(k)=0.0d0
         starq(k)=0.0d0
         staru(k)=0.0d0
         starv(k)=0.0d0
      end do

      do k=1,nmu
      do i=1,3
         si(k,i)=0.0d0
         sq(k,i)=0.0d0
         su(k,i)=0.0d0
         sv(k,i)=0.0d0
         si2(k,i)=0.0d0
         sq2(k,i)=0.0d0
         su2(k,i)=0.0d0
         sv2(k,i)=0.0d0
         nums(k,i)=0
         spi(k,i)=0.0d0
         spq(k,i)=0.0d0
         spu(k,i)=0.0d0
         spv(k,i)=0.0d0
      enddo
      enddo

      do i=1,3
         ti(i)=0.0d0
         tq(i)=0.0d0
         tu(i)=0.0d0
         tv(i)=0.0d0
         ti2(i)=0.0d0
         tq2(i)=0.0d0
         tu2(i)=0.0d0
         tv2(i)=0.0d0
	   numt(i)=0
         tis(i)=0.0d0
         tqs(i)=0.0d0
         tus(i)=0.0d0
         tvs(i)=0.0d0
      end do

      do i=1,nx
      do j=1,nx
      do k=1,nmu
         image(i,j,k)=0.0d0
         imagei(i,j,k)=0.0d0
         imageq(i,j,k)=0.0d0
         imageu(i,j,k)=0.0d0
         imagev(i,j,k)=0.0d0
         image2(i,j,k)=0.0d0
         numi(i,j,k)=0
      end do
      end do
      end do

      do i=1,nxhst
      do j=1,nxhst
         tihst(i,j)=0.0d0
         tqhst(i,j)=0.0d0
         tuhst(i,j)=0.0d0
         tvhst(i,j)=0.0d0
      end do
      end do

      binpvg=npvg
      bincvg=ncvg
      do i=1,ncvg
         do j=1,npvg
            ivg(i,j)=0.d0
            qvg(i,j)=0.d0
            uvg(i,j)=0.d0
            vvg(i,j)=0.d0
         enddo
      enddo
	ivgflux=0.d0
	qvgflux=0.d0
	uvgflux=0.d0
	vvgflux=0.d0
	ivgerr=0.d0
	qvgerr=0.d0
	uvgerr=0.d0
	vvgerr=0.d0


      do ir=1,nrg
      do it=1,ntg
      do ip=1,npg
         densarr(ir,it,ip)=0.d0
         massarr(ir,it,ip)=0.d0
         volarr(ir,it,ip)=0.d0
         pathI(ir,it,ip)=0.d0
         pathQ(ir,it,ip)=0.d0
         pathU(ir,it,ip)=0.d0
         pathV(ir,it,ip)=0.d0
      enddo
      enddo
      enddo


      return
      end


c     ********************************************************
