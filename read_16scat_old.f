      subroutine read_16scat(wave,filename)
 
c
c This program reads in the unformatted binary file written by the
c IDL routine load_16ampld.pro.  
c
c m.j. wolff and b. a. whitney, Space Science Institute, 2000/03/26
c
c history:
c 00/04/12 (baw) sorry, put include file to save memory
c        mike will hate this
c 00/04/11 (baw) get rid of wavelength index in arrays
c 00/07/28 (mjw): remove extraneous commas that drove lf95 crazy
c                   (occurrences of "write(*,*),something")
c 00/07/31 (mjw): add NWAV to read statement
c 00/08/01 (baw): fix diagnostic output
c 00/10/05 (mjw): explicitly define PEAKNEW
c 00/11/27 (mjw): update for new format of files (add DPHI, 
c                   g is no longer read in...)
c                 scattering quantities are real*4 (NOT angles, though)
c 00/11/29 (mjw): add normalization, which occurs if |integral-0.5| > epsilon
c                 update diagnostic output to reflect inverse order of
c		  cosine arrays (i.e., rlam(90 degrees) = rlam(1)
c
      implicit none

c ... scalar arguments
      character filename*(*)
      real *8 wave


c ... include files/definitions
      include 'tab.txt'
      include 'stokes.txt'

c ... local scalar
      integer nhalf,ipeak,ibpeak,ippeak,ibppeak
      integer i,iredo,j,k,l,m,nwav
      real*4 sum,epsilon
      real*8 XH,mH,avealb,avekapi,avekapq,tmp,ppeak,
     &     peaknew

c ... local arrays
      real*4 fx(MXDPHI),
     &     fxy(MXDPHI),sin_th,dh

c ... functions
      real*4 simp
      external simp

      data XH,mH,epsilon /0.69d0,1.67d-24,1.e-4/
c     XH = 0.69d0  (mass fraction of hydrogen)
c     mH = 1.67d-24 (grams of mH)
      pi = 4.d0*datan(1.d0)


      print*,' '
      write(12,*) ' '
      print*,' reading in dust properties from file ',filename
      write(12,*) ' reading in dust properties from file ',filename

      open(unit=66,file=filename,form='unformatted')

c ... get file array size info
      read(66) nwav,nthetab,ntheta,ndphi
      print*,' NWAV = ',nwav
      print*,' NTHETAB = ',nthetab
      print*,' NTHETA = ',ntheta
      print*,' NDPHI = ',ndphi

c ... check array sizes
      if (nwav.gt.1) stop 'nwav .gt. 1 NOT allowed'
      if (nthetab.gt.MXTHETB) stop 'increase MXTHETB'
      if (ntheta.gt.MXTHET) stop 'increase MXTHET'
      if (ndphi.gt.MXDPHI) stop 'increase MXDPHI'

c ... score the data
      read(66) wave
      print*,'wavelength',wave
      read(66) (cosbarr(i),i=1,nthetab)
      print*,'cosbarr(1),cosbarr(nthetab)',cosbarr(1),cosbarr(nthetab)
      read(66) (cosarr(i), i=1,ntheta)
      print*,'cosarr(1),cosarr(ntheta)',cosarr(1),cosarr(ntheta)
      read(66) (dphi(i), i=1,ndphi)
      print*,'dphi(1),dphi(ndphi)',dphi(1),dphi(ndphi)
      read(66) (kapd(i),i=1,nthetab)      !Cext
      read(66) (rlam(i),i=1,nthetab)      !Csca
      read(66) (kapd_q(i),i=1,nthetab)    !Cpol
      read(66) (Ccpol(i),i=1,nthetab)     !Ccpol
      read(66) (((((Fmat(j,k,l,m,i),j=1,4),k=1,4),
     &     l=1,nthetab),m=1,ntheta),i=1,ndphi)

      close(66)

      avealb=0.d0
      avekapi=0.d0
      avekapq=0.d0
      do i=1,nthetab
c        kapd is initially cext, rlam is initially cscat, kapd_q is Cpol
c        convert to quantities of interest
         rlam(i)=rlam(i)/kapd(i)
         kapd(i)=kapd(i)*XH/mH/1.d8
         kapd_q(i)=kapd_q(i)*XH/mH/1.d8
         Ccpol(i)=Ccpol(i)*XH/mH/1.d8
         avealb=avealb+rlam(i)	
         avekapi=avekapi+kapd(i)	
         avekapq=avekapq+kapd_q(i)	
c     print*,'cosbarr(i),i',cosbarr(i),i
      enddo
      avealb=avealb/float(nthetab)
      avekapi=avekapi/float(nthetab)
      avekapq=avekapq/float(nthetab)

      nhalf=nthetab/2
      print*,'average alb, alb(0deg), alb(90), alb(180deg) '
      print*, avealb
      print*, rlam(nthetab),rlam(nhalf),rlam(1)
      print*,'average kapd, kapd(0deg), kapd(90), kapd(180deg) '
      print*, avekapi
      print*, kapd(nthetab),kapd(nhalf),kapd(1)
      print*,'average kapd_q, kapd_q(0deg), kapd_q(90), kapd_q(180deg) '
      print*, avekapq
      print*, kapd_q(nthetab),kapd_q(nhalf),kapd_q(1)
      print*, ' '

      write(12,*) 'average alb, alb(0deg), alb(90), alb(180deg) '
      write(12,*) avealb
      write(12,*) rlam(nthetab),rlam(nhalf),rlam(1)
      write(12,*)'average kapd, kapd(0deg), kapd(90), kapd(180deg) '
      write(12,*) avekapi
      write(12,*) kapd(nthetab),kapd(nhalf),kapd(1)
      write(12,*)'avg. kapd_q, kapd_q(0deg), kapd_q(90), kapd_q(180deg)'
      write(12,*) avekapq
      write(12,*) kapd_q(nthetab),kapd_q(nhalf),kapd_q(1)
      write(12,*) ' '
 
 0100 iredo = 0

c  check that integral of F_11 over 4 pi is 1.  note we
c  are only integrating of 2 pi here...so integral should = 0.5
c  to precision epsilon
      do l=1,nthetab
         do i=1,ndphi
            do m=1,ntheta
               sin_th = sngl(dsqrt(1.d0-cosarr(m)**2))
               fxy(m) = Fmat(1,1,l,m,i)*sin_th
            end do
            dh = dacos(cosarr(1))-dacos(cosarr(2))
            fx(i) = simp(ntheta,MXDPHI,dh,fxy)
         enddo
         dh = dphi(2) - dphi(1)
         sum = simp(ndphi,MXDPHI,dh,fx)
         write(*,'(f10.5,a,f7.2)') sum,
     &        ' = integral of phs fn. over 2pi for thetab = ',
     &        dacos(cosbarr(l)) * 180.d0 / pi
         write(12,'(f10.5,a,f7.2)') sum,
     &        ' = integral of phs fn. over 2pi for thetab = ',
     &        dacos(cosbarr(l)) * 180.d0 / pi
         
c     ... renormalize, if necessary (e.g., in angular resolution)
c     ... baw, renormalize anyway, to divide by 4*pi
c     if (abs(1.0-sum).gt.epsilon) then
c     write(*,*) 'Renormalizing Scattering Matrix (F, not Z)'
c     write(12,*) 'Renormalizing Scattering Matrix (F, not Z)'
c     iredo = 1
         do i=1,ndphi
            do m=1,ntheta
               do j=1,4
                  do k=1,4
                     Fmat(j,k,l,m,i)=Fmat(j,k,l,m,i)/sum*pi*4.0
                  end do
               end do
            end do
         end do
c     endif
      enddo
c      if (iredo.eq.1) goto 0100

c     find peak of s11
      peak=0.
      ppeak=0.
      ippeak=0
      ibppeak=0
      ipeak=0
      ibpeak=0
      do l=1,nthetab
         do m=1,ntheta
            do i=1,ndphi
               if (Fmat(1,1,l,m,i).ge.peak) then
                  peak=Fmat(1,1,l,m,i)
                  ipeak=m
                  ibpeak=l
               endif
               tmp=-Fmat(1,2,l,m,i)/Fmat(1,1,l,m,i)
               if(tmp.ge.ppeak) then
                  ppeak=tmp
                  ippeak=m
                  ibppeak=l
               endif
            end do
         enddo
      enddo

      write(6,*) 'peak of s11',peak
      print*,'occurs at cos(thetb),cos(thet) = ',
     $ cosbarr(ibpeak),cosarr(ipeak)
      print*, 'peak of linear polarization ',ppeak
      print*,'occurs at cos(thetb),cos(thet) = ',
     $ cosbarr(ibppeak),cosarr(ippeak)

c     increase peak to account for polarized radiation
c      peaknew=-2.85*peak*ppeak*kapd_q(1)/kapd(1)
c      peaknew=peak*ppeak*(kapd(nhalf-1)-kapd_q(nhalf-1))/
c     $(kapd(nhalf-1))
c      peaknew=peak*ppeak*(kapd(nhalf+3)-kapd_q(nhalf+3))/
c     $(kapd(nhalf+3)+kapd_q(nhalf+3))
cc      peaknew=1.34*peak*ppeak
c     wild guess, just increase in proportion to the opacity variations.
      peaknew=peak*(abs(kapd_q(nhalf))+kapd(nhalf))/kapd(nhalf)
      print*,'peaknew ',peaknew
      if (peaknew.gt.peak) peak=peaknew
      print*,'peak increased to ',peak
      print*, ' '

      return
      end

      real function simp(ndata,nmax,h,f)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     this routine computes int( f(x) ) from a to b with mesh
c     being uniform.  the method employed is the alternate extended Simpson's 
c     rule (p. 117, Numerical Recipes in C)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none

      integer ndata,i,nmax
      real h,f(nmax),coeff,integral
 
      if(ndata.lt.8) stop '_simp requires at least 8 points' 

      integral = 0.
      
c ... compute coefficients for quadrature sum

      do 10 i=0,ndata-1
         if( (i.eq.0). or .(i.eq.(ndata-1)) ) then
            coeff = 17. / 48.
         else if( (i.eq.1). or .(i.eq.(ndata-2)) ) then
            coeff = 59. / 48.
         else if( (i.eq.2). or .(i.eq.(ndata-3)) ) then
            coeff = 43. / 48.
         else if( (i.eq.3). or .(i.eq.(ndata-4)) ) then
            coeff = 49. / 48.
         else
            coeff = 1.
         endif

         integral = integral + coeff*f(i+1)

 10   continue

      simp = h*integral

      return
      end
      


      
