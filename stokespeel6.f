      subroutine stokespeel6(sip,sqp,sup,svp,cost,sint
     1 ,hit,htot,pi,r2p,peak,ri1,bmu,costp,sintp)

c history:
c 99/03/22 (baw) put in 6-element matrix, use rejection method only
c 95/01/17 (mjw): add calls to ERRMSG,WRIMSG
c 95/12/08 (baw): sample from tables instead of analytic functions

c     This subroutine
c     calculates scattering using tabular phase functions.
c     this is for randomly oriented spheres or spheroids, which require 
c     four to six unique matrix elements
c     placed the following way:
c        |  S11 S12 0   0  | 
c        |  S12 S22 0   0  |
c        |  0  0   S33 S34 |
c        |  0  0  -S34 S44 |
c     Bohren&Huffman eqn (13.21).
c     note:  white and bh have sign of S34 switched.  I'm using B-H's 
c     convention here, but in the original Stokes routines, I used
c     White's.  

c     This matrix is rotated into and out of the scattering frame to
c     get stokes parameters in observer's frame
 
      implicit none

      character cmsger*70
c      include 'mrncoef.txt'
      include 'tab.txt'
      include 'random.txt'
 
      integer itab

      real*8 sip,sqp,sup,svp,cost,sint,pi,r2p
     1 ,hit,htot,peak,ri1,bmu,costp,sintp

      real*8 xran,cosi3,sini3,cosb2,ri3
     1 ,b,sinbt,sini2,bott,cosi2,sin2i3,sin2i2,cos2i3,cos2i2
     1 ,sin2,cos2,sin2cos1,cos2sin1,s11,s12,s22,s33,s34,s44
     1 ,a11,a12,a13,a21,a22
     1 ,a23,a24,a31,a32,a33,a34,a42,a43,a44,cosi1,sini1,sin2i1
     1 ,cos2i1,a,si,sq,su,sv,rprob

C     External Subroutines ..
      external errmsg, wrimsg
      
      real ran2

c     uses rejection method for sampling from
c     exact phase function; can be inefficient for
c     forward-peaked dust.
      
c5     continue
c      htot=htot+1.d0
      itab=int((1.d0+bmu)/dct)+1
      if (itab.lt.1.or.itab.gt.narr) print*,'whoops,itab',itab
      s11=s11arr(itab)
      s12=s12arr(itab)
      cosb2=bmu**2
      b=cosb2-1.d0

      if(dabs(bmu).gt.1.d0) then
         write(cmsger,'(a,f)')'bmu eq.1',bmu
         call ERRMSG('WARNING','STOKES',cmsger)
         if(bmu.gt.1.d0) then
            bmu=1.d0
            cosb2=1.d0
            b=0.d0
         else
            bmu=-1.d0
            cosb2=1.d0
            b=0.d0
         end if
      end if
      
      sinbt=dsqrt(1.d0-cosb2)
      
      
c     **** ri1 gt pi ****
c     don't try to combine this with ri1 lt pi.  there are several
c     slight differences and it will require more than one if statement
      
      if (ri1.gt.pi) then
         ri3=r2p-ri1
         cosi3=dcos(ri3)
         sini3=dsin(ri3)
         sin2i3=2.d0*sini3*cosi3
         cos2i3=2.d0*cosi3**2-1.d0
         a11=s11
         a12=s12*cos2i3
         a13=s12*sin2i3
c        rejection test
         rprob=(a11*sip+a12*sqp+a13*sup)/sip
         if(rprob.gt.peak) then
c        this shouldn't be the case since we measured peak of p1
            write(cmsger,'(a,f14.11,f14.11)')'rprob gt peak!',
     $		rprob,peak
            call ERRMSG('WARNING','STOKES',cmsger)
            peak=rprob
         end if
c         xran=ran2(i1)         
c         if(peak*xran.gt.rprob) go to 5
c         hit=hit+1
         a=rprob
         
         if(bmu.eq.1.d0) then
            write(cmsger,'(a,f)')'bmu eq.1',bmu
            call ERRMSG('WARNING','STOKES',cmsger)
            go to 10
         else if(bmu.eq.-1.d0) then
               write(cmsger,'(a,f)')'bmu eq.-1',bmu
               call ERRMSG('WARNING','STOKES',cmsger)
               sini3=1.d0
               sini2=1.d0
               cosi3=0.d0
               cosi2=0.d0
c               cost=-costp
c               sint=sintp
         else
c           anything with p on the end is the primed angle, which is
c           the angle from whence it came. oops--except cosbwp.
c            cost=costp*bmu+sintp*sinbt*cosi3
c           write(6,*)'costp,bmu,sintp,sinbt,cosi3,cost'
c           write(6,*)costp,bmu,sintp,sinbt,cosi3,cost
c           everytime we do an arcsin or arccos we have to test to see
c           if roundoff error has caused cos or sin to be slightly
c           greater than 1.  then we have to set them to 1.
            if (dabs(cost).lt.1.d0) then
c               sint=abs(sqrt(1-cost**2))
               sini2=sini3*sintp/sint
c              if (sini2.lt.0) write(6,*)'sini2 lt 0',sini2
               bott=sint*sinbt
               cosi2=costp/bott-cost*bmu/bott
               if (dabs(cosi2).gt.1.0001d0) then
                  write(cmsger,'(a,f)')'cosi2 big',cosi2
                  call ERRMSG('WARNING','STOKES',cmsger)
               endif
            else
               sint=0.d0
               call ERRMSG('WARNING','STOKESPEEL','sint,sini2=0')
               sini2=0.0d0
               if (cost.ge.1.d0) then
                  cosi2=-1.d0
                  cost=1.d0
               endif
               if (cost.le.-1.d0) then
                  cosi2=1.d0
                  cost=-1.d0
               endif
            end if
         endif

c         cosdph=-cosi2*cosi3+sini2*sini3*bmu
c         if (abs(cosdph).gt.1) then
c            if(abs(cosdph).gt.1.001) then
c               write(cmsger,'(a,f)')'abs(cosdph).gt.1.001',cosdph
c               call ERRMSG('WARNING','STOKES',cmsger)
c            end if
c            if(cosdph.gt.1) then
c               cosdph=1.0
c            else 
c               cosdph=-1.0
c            end if
c         end if
c         phi=phip+acos(cosdph)
c         if(phi.gt.r2p) phi=phi-r2p
c         if(phi.lt.0.) phi=phi+r2p
      
         sin2i2=2.d0*sini2*cosi2
         cos2i2=2.d0*cosi2**2-1.d0
         sin2=sin2i2*sin2i3
         cos2=cos2i2*cos2i3
         sin2cos1=sin2i2*cos2i3
         cos2sin1=cos2i2*sin2i3

         s22=s22arr(itab)
         s33=s33arr(itab)
         s34=s34arr(itab)
         s44=s44arr(itab)

         a21=s12*cos2i2
         a22=s22*cos2-s33*sin2
         a23=s22*cos2sin1+s33*sin2cos1
         a24=s34*sin2i2
         a31=-s12*sin2i2
         a32=-s22*sin2cos1-s33*cos2sin1
         a33=-s22*sin2+s33*cos2
         a34=s34*cos2i2
         a42=s34*sin2i3
         a43=-s34*cos2i3
         a44=s44
 
c       **** ri1 lt pi ****
      else   
 
         cosi1=dcos(ri1)
         sini1=dsin(ri1)
         sin2i1=2.d0*sini1*cosi1
         cos2i1=2.d0*cosi1**2-1.d0
         a11=s11
         a12=s12*cos2i1
         a13=-s12*sin2i1
c        rejection test
         rprob=(a11*sip+a12*sqp+a13*sup)/sip
         if(rprob.gt.peak) then
c           this shouldn't be the case since we measured peak of p1
            write(cmsger,'(a,f14.11,f14.11)')'rprob gt peak!',
     $		rprob,peak
            call ERRMSG('WARNING','STOKES',cmsger)
            peak=rprob
         end if
c     xran=ran2(i1)         
c     if(peak*xran.gt.rprob) go to 5
c     hit=hit+1
         a=rprob
         
         if(bmu.eq.1.d0) then
            write(cmsger,'(a,f)')'bmu eq.1',bmu
            call ERRMSG('WARNING','STOKES',cmsger)
            go to 10
         else if(bmu.eq.-1.d0) then
            write(cmsger,'(a,f)')'bmu eq.-1',bmu
            call ERRMSG('WARNING','STOKES',cmsger)
            sini1=1.d0
            sini2=1.d0
            cosi3=0.d0
            cosi2=0.d0
c            cost=-costp
c            sint=sintp
         else
         
c           anything with p on the end is the primed angle, which is
c           the angle from whence it came. oops--except cosbwp.
 
c            cost=costp*bmu+sintp*sinbt*cosi1
c           write(6,*)'costp,bmu,sintp,sinbt,cosi1,cost'
c           write(6,*)costp,bmu,sintp,sinbt,cosi1,cost
c           everytime we do an arcsin or arccos we have to test to see
c           if roundoff error has caused cos or sin to be slightly
c           greater than 1.  then we have to set them to 1.
            if (dabs(cost).lt.1.d0) then
c               sint=abs(sqrt(1-cost**2))
               sini2=sini1*sintp/sint
c              if (sini2.lt.0) write(6,*)'sini2 lt 0',sini2
               bott=sint*sinbt
               cosi2=costp/bott-cost*bmu/bott
               if (dabs(cosi2).gt.1.0001d0) then
                  write(cmsger,'(a,f)')'cosi2 big',cosi2
                  call ERRMSG('WARNING','STOKES',cmsger)
               endif
            else   
c               sint=0
               call ERRMSG('WARNING','STOKESPEEL','sint,sini2=0')
               sini2=0.0d0
               if (cost.ge.1.d0) then
                  cosi2=-1.d0
c                  cost=1.
               endif
               if (cost.le.-1.d0) then
                  cosi2=1.d0
c                  cost=-1.
               endif
            end if
         endif
c         cosdph=-cosi1*cosi2+sini1*sini2*bmu
c         if (abs(cosdph).gt.1) then
c            if(abs(cosdph).gt.1.001)then
c               write(cmsger,'(a,f)')'abs(cosdph).gt.1.001',cosdph
c               call ERRMSG('WARNING','STOKES',cmsger)
c            end if
c            if(cosdph.gt.1) then
c               cosdph=1.0
c            else
c               cosdph=-1.0
c            end if
c         end if
         
c         phi=phip-acos(cosdph)
c         if(phi.gt.r2p) phi=phi-r2p
c         if(phi.lt.0.) phi=phi+r2p
         
         sin2i2=2.d0*sini2*cosi2
         cos2i2=2.d0*cosi2**2-1.d0
         sin2=sin2i2*sin2i1
         cos2=cos2i2*cos2i1
         sin2cos1=sin2i2*cos2i1
         cos2sin1=cos2i2*sin2i1

         s22=s22arr(itab)
         s33=s33arr(itab)
         s34=s34arr(itab)
         s44=s44arr(itab)

         a21=s12*cos2i2
         a22=s22*cos2-s33*sin2
         a23=-s22*cos2sin1-s33*sin2cos1
         a24=-s34*sin2i2
         a31=s12*sin2i2
         a32=s22*sin2cos1+s33*cos2sin1
         a33=-s22*sin2+s33*cos2
         a34=s34*cos2i2
         a42=-s34*sin2i1
         a43=-s34*cos2i1
         a44=s44
 
      end if
 
      si=(a11*sip+a12*sqp+a13*sup)
      sq=(a21*sip+a22*sqp+a23*sup+a24*svp)
      su=(a31*sip+a32*sqp+a33*sup+a34*svp)
      sv=(a42*sqp+a43*sup+a44*svp)

c     write(6,*)'si,sq,su',si,sq,su,nscat
 
      sip=si
      sqp=sq
      sup=su
      svp=sv
c      cosp=cos(phi)
c      sinp=sin(phi)
      
      
 10   continue
      return
      end
 
c     ***********************************************************************
