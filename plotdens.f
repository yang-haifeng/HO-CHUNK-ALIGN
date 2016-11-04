      program plotdens4

      implicit none

      integer nrg,ntg,npg,nrgmax,ntgmax,npgmax

      parameter (nrgmax=400,ntgmax=499,npgmax=1)
      
      real*8 densarr(nrgmax,ntgmax,npgmax),rarr(nrgmax)
     $     ,phiarr(npgmax),thetarr(ntgmax)

      real xyarr(200,200),x1,y1,xmax(4),dx,r1,thet1,r2,thet2,x2,y2,t

      integer nrg1,ntg1,npg1,nx,ix,iy,ir2,it2,ip,n,i,ir1,it1,ir,it,
     1     ntot,tmp1,tmp2

      character*20 filename

      open(unit=15,file='dens.unf',status='unknown'
     1 ,form='unformatted')
      read(15) nrg,ntg,npg
	print*,nrg,ntg,npg
c     do test on this another day
      read(15) (rarr(ir1),ir1=1,nrg)
      read(15) (thetarr(it1),it1=1,ntg)
      read(15) (phiarr(ip),ip=1,npg)
      read(15) (((densarr(ir1,it1,ip),ir1=1,nrg),it1=1,ntg),ip=1,npg)
      close(15)      

      print*,'done reading in the data'
      n=200
c      print*,'enter 4 values xmax'
c      read*,xmax(1),xmax(2),xmax(3),xmax(4)
      xmax(1)=.625
      xmax(2)=12.5
      xmax(3)=250.
      xmax(4)=5000.
      xmax(1)=1.
      xmax(2)=5.
      xmax(3)=500.
      xmax(4)=5000.

      do i=1,4
      dx=xmax(i)/float(n)
      if (i.eq.1) filename='densxz1.dat'
      if (i.eq.2) filename='densxz2.dat'
      if (i.eq.3) filename='densxz3.dat'
      if (i.eq.4) filename='densxz4.dat'
      do ix=1,n
c     x=dx*(float(ix)-0.5)
         x1=dx*(float(ix)-1)
         x2=dx*(float(ix))
         do iy=1,n
            y1=dx*(float(iy)-1)
            y2=dx*(float(iy))
c     y=dx*(float(iy)-0.5)
c     if (iy.eq.1.and.ix.eq.1) print*,'x0,y0',x,y
         r1=sqrt(x1**2+y1**2)
         r2=sqrt(x2**2+y2**2)
         thet1=atan2(x1,y1)
         thet2=atan2(x2,y2)
c     if(x1.eq.0) print*,'x1=0,thet1',thet1
         call locate(rarr,nrg,dble(r1),ir1)
         call locate(rarr,nrg,dble(r2),ir2)
         if (ntg.gt.1) then
            call locate(thetarr,ntg,dble(thet1),it1)
            call locate(thetarr,ntg,dble(thet2),it2)
         else
            it1=1
            it2=1
         endif
         if (thet1.eq.0) it1=1
         if (it1.gt.it2) then
            tmp2=it2
            tmp1=it1
            it1=tmp2
            it2=tmp1
c     print*,'shit, it2,it1',it2,it1
c     print*,'x1,x2,y1,y2,it1,it2,ir1,ir2'
c     print* ,x1,x2,y1,y2,it1,it2,ir1,ir2
         endif
         if (ir1.gt.ir2) then
            print*,'shit, ir1,ir2',ir1,ir2
            print*,'x1,x2,y1,y2,it1,it2,ir1,ir2'
            print* ,x1,x2,y1,y2,it1,it2,ir1,ir2
         endif
c         if (x1.eq.0) print*,'x1=0,it1',it1
c         if (it2.eq.1) it2=2
c         if (it1.eq.1) it1=2
         t=0.d0
         ntot=0
         do it=it1,it2
            do ir=ir1,ir2
c     if (tdustarr(ir,it,1).eq.1000.d0) tdustarr(ir,it,1)=0.d0
               t=t+densarr(ir,it,1)
               ntot=ntot+1
            enddo
         enddo
         if (ntot.eq.0) print*,'ntot=0,screwed up,x,y',x1,x2,y1,y2
         xyarr(ix,iy)=t/float(ntot)
      end do
      end do
      open(unit=15,file=filename,status='unknown')
      do iy=1,n
         write(15,*) (xyarr(ix,iy),ix=1,n)
      enddo
      close(15)
      enddo

      print*,'n=',n
      print*,thetarr(1),thetarr(2)

      end

c     ***********************************************

      SUBROUTINE locate(xx,n,x,j)
c     from numerical recipes
c     searches an ordered table, using bisection
      INTEGER j,n
      REAL*8 x,xx(n)
      INTEGER jl,jm,ju
      jl=0
      ju=n+1
10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      goto 10
      endif
      j=jl
      return
      END
      
      
