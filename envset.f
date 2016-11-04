      subroutine envset

c     2001/02/17 (baw) constants for envelope assuming TSC
c     sets cavity shape parameters also.
c     2003/12/26 set c1e,c2e when hole opening angle=90.

      implicit none

      include 'tts.txt'
      include 'opacin.txt'

      real*8 pi,const,sphmass,rd,rdcgs,rmincgs,z1max,r1max,rtmp,r2max

c ... begin g77
        real*8 dcosd,dsind,dtand
        external dcosd,dsind,dtand
c ... end g77

      pi=4.d0*datan(1.d0)

c     infall;  Rd is disk radius, units of rmin 
c     BAW 7/5/98 make rd = rchole = rc
c     rd = rchole
c     BAW 12/15/99 make rd=rc
      rd=rc
      rdcgs=rd*rstar*rsol
c     rhoe0 is factor in front of density distribution, given by
c     infall calculation.
c     rhoe0=(infallrate)/(4*pi)/sqrt(G*Mcore)/rd**1.5
      const=1.d0/3.15576d7*dsqrt(msol)/(4.d0*pi)/dsqrt(6.67259d-8)
      write(6,*)'const',const
      rhoe0=const*rate/dsqrt(massc)/rdcgs**1.5
      sphmass=2.d0/3.d0/3.15576d7/dsqrt(msol)*rate/dsqrt(massc)/
     1     dsqrt(2.d0*6.67259d-8)
      sphmass=sphmass*(rmax*rstar*rsol)**1.5
      write(6,*) 'rhoe0 of envelope',rhoe0
      rmincgs=rstar*rsol
      write(6,*) 'rmincgs',rmincgs
      windmu0=dcosd(thetmu0)
      print*,'windmu0 ',windmu0

c     ambient density
c     rhoamb=0.d0

c     hole in bubble
c     roa=rmax*dtand(buboa)
      cosbuboa=dcosd(buboa)
      
c     1999
c     outflow boundaries
c     z1max is height where opening angle is measured.  take as outer
c     bound.
      z1max=rmax
      if (thet1.lt.89.999) then
         r1max=z1max*dtand(thet1)
c     z=a+b*x**beta
c     c1=b
c     z01=a
c     ex1=beta
c     r=x  (r is cylindrical radius)
         c1e=(z1max-z01)/r1max**ex1
         print*,'z1max,z01,r1max,ex1,c1e'
         print*,z1max,z01,r1max,ex1,c1e
c     z01 is input
         if(z01.lt.0.d0) then
            rtmp=(-z01/c1e)**(1.d0/ex1)
            print*,'hole 1 intersects disk at ',rtmp/autors, ' AU'
         endif
         if(ipoly.eq.1) then
            r2max=z1max*dtand(thet2)
            c2e=(z1max-z01)/r2max**ex2
            print*,'r2max,ex2,c2e',r2max,ex2,c2e
c     z02 is input
            if(z02.lt.0.d0) then
               rtmp=(-z02/c2e)**(1.d0/ex2)
               print*,'hole 2 intersects disk at ',rtmp/autors, ' AU'
            endif
         end if
      else
         c1e=0.
         c2e=0.
      endif
      
      return
      end
