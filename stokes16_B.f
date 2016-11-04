      subroutine stokes16(costb,sintb,cospb,sinpb,phib,sqpb,supb)


c history:
c 00/12/11 (baw) new tables have 3 dimensions, theta_prime, theta,phi,
c     so use trilinear interpolation to sample!  this is a major
c     revision.  for now, assume field is in z-direction.  in this
c     case, there is no rotation in and out of scattering frame,
c     just sampling for new theta,phi given theta_prime.  the next
c     modification will have to take into account field in different
c     direction.  
c 00/04/10 (baw) sampling from tables is different from before. use
c     bilinear interpolation of table which depends on incident angle
c     w.r.t B-field and scattering angle.
c 00/03/25 (baw) put in 16-element matrix
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
      include 'stokes.txt'
      include 'tab.txt'
      include 'random.txt'
 
      integer iinc,isca,iphi

      real*8 xran,bmu,costp,sintp,phip,delphi,delph
     1 ,a11,a12,a13,a14,a21,a22,a23,a24,a31,a32,a33,a34
     1 ,a41,a42,a43,a44
     1 ,a,si,sq,su,sv,rprob
     1 ,costb,sintb,cospb,sinpb,phib,sqpb,supb

C     External Subroutines ..
      external errmsg, wrimsg
      
      real ran2
      real*8 trilinfmat,check

      character*20 routine

      routine='stokes16_B'

c     costp=dabs(cost)
      costp=costb
      phip=phib
      call locate(cosbarr,MXTHETB,nthetab,costp,iinc)

c     uses rejection method for sampling from
c     exact phase function; can be inefficient for
c     forward-peaked dust.
      
5     continue
      htot=htot+1.d0
      xran=ran2(i1)         
      bmu=1.d0-2.d0*xran
c     note:  bmu is not scattering angle, as in previous version of
c     stokes.  it is new theta measured from B-field axis.
      delphi=r2p*ran2(i1)
      if (delphi.gt.r2p) delphi=r2p
      if (delphi.lt.0.d0) delphi=0.d0
      call locate(cosarr,MXTHET,ntheta,bmu,isca)
      call locate(dphi,MXDPHI,ndphi,delphi,iphi)
      a11=trilinfmat(1,1,iinc,isca,iphi,costp,bmu,delphi)
      a12=trilinfmat(1,2,iinc,isca,iphi,costp,bmu,delphi)
      a13=trilinfmat(1,3,iinc,isca,iphi,costp,bmu,delphi)
      a14=trilinfmat(1,4,iinc,isca,iphi,costp,bmu,delphi)

c     rejection test
      rprob=(a11*sip+a12*sqpb+a13*supb+a14*svp)/sip
      if(rprob.gt.peak) then
c     this shouldn't be the case since we measured peak of p1
         write(cmsger,'(a,f14.11,f14.11)')'rprob gt peak!',
     $        rprob,peak
         call ERRMSG('WARNING','STOKES',cmsger)
c         print*,'sip,sqp,sup,svp',sip,sqp,sup,svp
         peak=rprob
      end if
      xran=ran2(i1)         
      if(peak*xran.gt.rprob) go to 5
      hit=hit+1
      a=rprob
         
      a21=trilinfmat(2,1,iinc,isca,iphi,costp,bmu,delphi)
      a22=trilinfmat(2,2,iinc,isca,iphi,costp,bmu,delphi)
      a23=trilinfmat(2,3,iinc,isca,iphi,costp,bmu,delphi)
      a24=trilinfmat(2,4,iinc,isca,iphi,costp,bmu,delphi)
      a31=trilinfmat(3,1,iinc,isca,iphi,costp,bmu,delphi)
      a32=trilinfmat(3,2,iinc,isca,iphi,costp,bmu,delphi)
      a33=trilinfmat(3,3,iinc,isca,iphi,costp,bmu,delphi)
      a34=trilinfmat(3,4,iinc,isca,iphi,costp,bmu,delphi)
      a41=trilinfmat(4,1,iinc,isca,iphi,costp,bmu,delphi)
      a42=trilinfmat(4,2,iinc,isca,iphi,costp,bmu,delphi)
      a43=trilinfmat(4,3,iinc,isca,iphi,costp,bmu,delphi)
      a44=trilinfmat(4,4,iinc,isca,iphi,costp,bmu,delphi)

c     sampled from a so divide by a.
      si=(a11*sip+a12*sqpb+a13*supb+a14*svp)/a
      sq=(a21*sip+a22*sqpb+a23*supb+a24*svp)/a
      su=(a31*sip+a32*sqpb+a33*supb+a34*svp)/a
      sv=(a41*sip+a42*sqpb+a43*supb+a44*svp)/a
      
c     write(6,*)'si,sq,su',si,sq,su,nscat
 
      sip=si
      sqpb=sq
      supb=su
      svp=sv
c     print*,'in stokes16_B,sip',sip

      phib=delphi+phip
      if(phib.gt.r2p) phib=phib-r2p
      if(phib.lt.0.d0) phib=phib+r2p
      cospb=dcos(phib)
      sinpb=dsin(phib)

      costb=check(bmu,routine)
      sintb=dsqrt(1.d0-costb**2)
      
 10   continue
      return
      end
 

