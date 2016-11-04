      subroutine reatab(cfile,rlam,kapd,wave)

      implicit none

C     .. Scalar Arguments ..
      character cfile*50

C     .. Local Scalars ..
      integer iopar, narrin, nar, i
      character cline*70
      real*8 sum, thet, r, v, rlam, rjunk, kapd, rjunk2, ppeak, tmp,wave

      parameter(nar=181)

C     .. Local Arrays ..
      real*8 s11in(nar),s12in(nar),s33in(nar),s34in(nar),s22in(nar),
     1   s44in(nar),cosarrin(nar),s(nar),q(nar),
     1   s11tmp(nar),s12tmp(nar),s33tmp(nar),s34tmp(nar),s22tmp(nar),
     1   s44tmp(nar),cosarrtmp(nar)

C     .. External Subroutines ..
        external errmsg, wrimsg

C     .. Intrinsic Functions ..

C     .. Data statements ..
      data iopar /18/

c     .. Include statements ..
      include 'tab.txt'
      include 'stokes.txt'

c     testing
      
      open(unit=iopar,file=cfile,status='old')
	print*,' '

      narrin=0
      if (nelems.eq.4) then
         READ(IOPAR,FMT=5,ERR=1000)CLINE
         CALL WRIMSG(' ',CLINE)
c	  read(iopar,5)cline
 5       format(a70)
         do i=1,1000
            read(iopar,10,END=999,ERR=1000) thet,
     &      s11tmp(i),s12tmp(i),s33tmp(i),s34tmp(i)
c           planning to run stokes6, with 6-element matrix
c           so set redundant elements
            s22tmp(i) = s11tmp(i)
            s44tmp(i) = s33tmp(i)
c	      tharr(i)=thet*pi/180.
	      cosarrtmp(i)=dcos(thet*pi/180.d0)
	      narrin=narrin+1 
c	      sum=sum+s11arr(i)*sin(tharr(i))
	   end do
 10      format(1x,5(1pe11.4,1x))
      else if (nelems.eq.6) then
	   print*,'reading in 6-element dust scattering matrix'
	   do i=1,2
	      read(iopar,fmt=15,err=1000)cline
 15         format(a70)
            print*,cline
         enddo
	   read(iopar,*,err=1000) wave
	   print*,'wavelength: ',wave
	   read(iopar,*,err=1000) rjunk
	   print*, 'integral (Cext*n(a)da)/n_H (micron^2/H) ',rjunk
	   read(iopar,*,err=1000) rjunk2
	   print*, 'integral (Csca*n(a)da)/n_H (micron^2/H) ',rjunk2
	   rlam=rjunk2/rjunk
	   print*, 'albedo ',rlam
	   read(iopar,*,err=1000) kapd
	   print*, 'opacity (cm^2/g of gas) ',kapd
	   read(iopar,*,err=1000) rjunk
	   print*, '<cos(theta)> ',rjunk
	   write(12,*) 'wavelength, opacity, albedo, g '
         write(12,*) wave,rlam,kapd,rjunk
	   do i=1,2
	      read(iopar,fmt=15,err=1000)cline
         enddo
	   do i=1,1000
            read(iopar,*,END=999,ERR=1000) thet,
     &      s11tmp(i),s22tmp(i),s33tmp(i),s44tmp(i),
     &      s12tmp(i),s34tmp(i)
c           if file tabulates polarization instead of s12, do the following:
c	      s12tmp(i)=-s12tmp(i)*s11tmp(i)
	      cosarrtmp(i)=cos(thet*pi/180.d0)
	      narrin=narrin+1 
	   end do
      endif
 
 999  continue

	close(iopar)

c     now have to reverse order of tables
      do i=1,narrin
	   cosarrin(i)=cosarrtmp(narrin-i+1)
	   s11in(i)=s11tmp(narrin-i+1)
	   s22in(i)=s22tmp(narrin-i+1)
	   s33in(i)=s33tmp(narrin-i+1)
	   s44in(i)=s44tmp(narrin-i+1)
	   s12in(i)=s12tmp(narrin-i+1)
	   s34in(i)=s34tmp(narrin-i+1)
	enddo

	dct=2.d0/dfloat(narr)
	do i=1,narr
	   cosarr(i)=-1.+(dfloat(i)-0.5d0)*dct    !monotonically increasing
	                                       !for spline fits
	enddo

      if (nar.lt.narrin) print*,'messed up, reatab, nar, narrin', 
     1  nar,narrin


c     not very elegant....
	call onespl (cosarrin,s11in,s,q,narrin,nar)
	sum=0
	do i=1,narr
	   r=cosarr(i)
	   call spline(r,cosarrin,s11in,s,narrin,nar,v)
	   s11arr(i)=v
c	   sum=sum+v*dsin(dacos(r))*dct
	   sum=sum+0.5*v*dct
	enddo
	call onespl (cosarrin,s22in,s,q,narrin,nar)
	do i=1,narr
	   r=cosarr(i)
	   call spline(r,cosarrin,s22in,s,narrin,nar,v)
	   s22arr(i)=v
	enddo
	call onespl (cosarrin,s33in,s,q,narrin,nar)
	do i=1,narr
	   r=cosarr(i)
	   call spline(r,cosarrin,s33in,s,narrin,nar,v)
	   s33arr(i)=v
	enddo
	call onespl (cosarrin,s44in,s,q,narrin,nar)
	do i=1,narr
	   r=cosarr(i)
	   call spline(r,cosarrin,s44in,s,narrin,nar,v)
	   s44arr(i)=v
	enddo
	call onespl (cosarrin,s12in,s,q,narrin,nar)
	do i=1,narr
	   r=cosarr(i)
	   call spline(r,cosarrin,s12in,s,narrin,nar,v)
	   s12arr(i)=v
	enddo
	call onespl (cosarrin,s34in,s,q,narrin,nar)
	do i=1,narr
	   r=cosarr(i)
	   call spline(r,cosarrin,s34in,s,narrin,nar,v)
	   s34arr(i)=v
	enddo

	open(unit=15,status='unknown',file='dustfine.dat')
	do i=1,narr
	   write(15,30) cosarr(narr-i+1),s11arr(narr-i+1),s22arr(narr-i+1)
     1   ,s33arr(narr-i+1),s44arr(narr-i+1),s12arr(narr-i+1)
	1   ,s34arr(narr-i+1)
	enddo
	close(15)
 30      format(1x,7(1pe11.4,1x))

	peak=0.
	ppeak=0.
	do i=1,narr
	   s11arr(i)=s11arr(i)/sum
	   s12arr(i)=s12arr(i)/sum
	   s22arr(i)=s22arr(i)/sum
	   s33arr(i)=s33arr(i)/sum
	   s34arr(i)=s34arr(i)/sum
	   s44arr(i)=s44arr(i)/sum
	   if(s11arr(i).ge.peak) peak=s11arr(i)
	   tmp=-s12arr(i)/s11arr(i)
	   if(tmp.ge.ppeak) ppeak=tmp
	end do

	write(6,*) 'integral of s11',sum
	write(6,*) 'peak of s11',peak
	print*, 'peak of linear polarization ',ppeak
	write(6,*)'narr',narr

      print*, ' ' 
	return
 1000 CONTINUE
      CALL ERRMSG('FATAL','REATAB',' Error reading .par file')
 1010	return
	end
