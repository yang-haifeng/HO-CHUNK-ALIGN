      subroutine plancknu(t,nu,bnu)
 
c     calculate planck function, 
c     for given temperature, frequency
      
      implicit none

      real*8 nu,h,k,bnu,t,xpn,hnkt,c
         
      data c /2.99792458d10/
      data h /6.6260755d-27/
      data k /1.380658d-16/

      hnkt=h/k*nu/t
      if(hnkt.gt.170.d0) hnkt=170.d0
      xpn=dexp(hnkt)
      bnu=2.d0*h*nu*nu/c*nu/c/(xpn-1.d0)
	
      return
      end

