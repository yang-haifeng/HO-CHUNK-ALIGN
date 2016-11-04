      subroutine stokespeel16(sip,sqp,sup,svp,cost,coste
     1 ,hit,htot,pi,r2p,peak,phip,phie)

c history:
c 00/12/11 (baw) new tables have 3 dimensions, theta_prime, theta,phi,
c     so use trilinear interpolation to sample!  this is a major
c     revision.  for now, assume field is in z-direction.  in this
c     case, there is no rotation in and out of scattering frame,
c     just sampling for new theta,phi given theta_prime.  the next
c     modification will have to take into account field in different
c     direction.
c 00/04/11 (baw) put in 16-element matrix
c 99/03/22 (baw) put in 6-element matrix, use rejection method only
c 95/01/17 (mjw): add calls to ERRMSG,WRIMSG
c 95/12/08 (baw): sample from tables instead of analytic functions

c     This subroutine
c     calculates scattering using tabular phase functions.
c     This is required for dichroic scattering or non-spherical 
c     particles
c        |  S11 S12 S13 S14 | 
c        |  S12 S22 S23 S34 |
c        |  S31 S32 S33 S34 |
c        |  S41 S42 S43 S44 |

      implicit none

      character cmsger*70
c     include 'mrncoef.txt'
      include 'tab.txt'
      include 'random.txt'
 
      integer iinc,isca,iphi

      real*8 sip,sqp,sup,svp,pi,r2p
     1 ,hit,htot,peak,costp,sintp,coste,phie
     1 ,phip,delphi,cost

      real*8 xran
     1 ,a11,a12,a13,a21,a22
     1 ,a41,a14
     1 ,a23,a24,a31,a32,a33,a34,a42,a43,a44
     1 ,a,si,sq,su,sv,rprob

C     External Subroutines ..
      external errmsg, wrimsg
      
      real ran2
      real*4 trilinfmat

      delphi=phie-phip
      if (delphi.lt.0.d0) delphi=delphi+r2p
      if (delphi.gt.r2p) delphi=delphi-r2p
c     if (delphi.gt.pi)  delphi=r2p-delphi

      costp=(cost)
c     if (cost.lt.0) coste=-coste

      call locate(cosbarr,MXTHETB,nthetab,costp,iinc)
      call locate(cosarr,MXTHET,ntheta,coste,isca)
      call locate(dphi,MXDPHI,ndphi,delphi,iphi)
      a11=trilinfmat(1,1,iinc,isca,iphi,costp,coste,delphi)
      a12=trilinfmat(1,2,iinc,isca,iphi,costp,coste,delphi)
      a13=trilinfmat(1,3,iinc,isca,iphi,costp,coste,delphi)
      a14=trilinfmat(1,4,iinc,isca,iphi,costp,coste,delphi)

c      rprob=(a11*sip+a12*sqp+a13*sup+a14*svp)/sip
c      if(rprob.gt.peak) then
c     this shouldn't be the case since we measured peak of p1
c         write(cmsger,'(a,f14.11,f14.11)')'rprob gt peak!',
c     $        rprob,peak
c         call ERRMSG('WARNING','STOKES',cmsger)
c         peak=rprob
c      end if

      a21=trilinfmat(2,1,iinc,isca,iphi,costp,coste,delphi)
      a22=trilinfmat(2,2,iinc,isca,iphi,costp,coste,delphi)
      a23=trilinfmat(2,3,iinc,isca,iphi,costp,coste,delphi)
      a24=trilinfmat(2,4,iinc,isca,iphi,costp,coste,delphi)
      a31=trilinfmat(3,1,iinc,isca,iphi,costp,coste,delphi)
      a32=trilinfmat(3,2,iinc,isca,iphi,costp,coste,delphi)
      a33=trilinfmat(3,3,iinc,isca,iphi,costp,coste,delphi)
      a34=trilinfmat(3,4,iinc,isca,iphi,costp,coste,delphi)
      a41=trilinfmat(4,1,iinc,isca,iphi,costp,coste,delphi)
      a42=trilinfmat(4,2,iinc,isca,iphi,costp,coste,delphi)
      a43=trilinfmat(4,3,iinc,isca,iphi,costp,coste,delphi)
      a44=trilinfmat(4,4,iinc,isca,iphi,costp,coste,delphi)         

      si=(a11*sip+a12*sqp+a13*sup+a14*svp)
      sq=(a21*sip+a22*sqp+a23*sup+a24*svp)
      su=(a31*sip+a32*sqp+a33*sup+a34*svp)
      sv=(a41*sip+a42*sqp+a43*sup+a44*svp)

c     write(6,*)'si,sq,su',si,sq,su,nscat
 
      sip=si
      sqp=sq
      sup=su
      svp=sv

c     already know what theta and phi are so don't need to calculate them 
c     as in regular stokes routine
      
 10   continue
      return
      end
 
