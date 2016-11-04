c     ************************************************************
      subroutine densdisk(rad,sint,cost,phi,dens)

c     calculates density in disk with Gaussian scale height. 
c     called during grid setup so doesn't need to be optimized.

c  history:
c  00/09/06 (baw): write subroutine, modified from opacdiskg.f

      implicit none
      include 'tts.txt' 
      real*8 zp,rp,opacd,rad,sint,phi,dens,cost,gexp,zmr

      rp=rad*sint
      zp=rad*cost
      zmr=z1*(rp/rmin)**b

      if(rad.gt.rmaxd.or.rhod0.eq.0.d0.or.rad.lt.rmind.
     $     or.zmr.lt.1.d-10) then
         dens=0.d0
c     print*,'rad,rhod0 ',rad,rhod0
         return
      end if

c     print*,'sint,cost',sint,cost

      gexp=(zp/zmr)**2
      if(gexp.gt.1.d-15) then
         if(gexp.lt.120.d0) then
	    dens=rhod0*(rp/rmin)**(-a)*dexp(-0.5d0*gexp)
         else
	    dens=0.d0
         end if
      else
         dens=rhod0*(rp/rmin)**(-a)
      end if
      
      return
      end
      
c     **********************************************************
