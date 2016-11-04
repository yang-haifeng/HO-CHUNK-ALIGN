      subroutine thermset

c     set thermal constants

      implicit none

      include 'tts.txt'
      include 'grid.txt'
      include 'taunum.txt'

      if (itherm.eq.0) then
         ns=np+npout
         nd=0
         ne=0
      else
c     calculate luminosity of star
         nu=2.99792458d14/wave  !in microns
         call plancknu(tstar,nu,bnu)
         print*,'stellar Bnu ',bnu
c     bnu=bnu*4*pi     !convert from energy/s/cm2/Sr to energy/s/cm2
         rstarcgs=rstar*rsol
         area=pi*(rstarcgs)**2
         normstar=area*bnu      !this is an energy/s; normalization 
                                !for output
         print*,'stellar energy at this wavelength at earth',normstar
         print*,'luminosity ',4.d0*area/3.826d33*tstar**4*5.669d-5 
                                !energy/s
         normstar1=normstar*4.d0*pi !flux=pi*bnu area=4*pi*rstar**2
         print*,'luminosity at this wavelength ',normstar1
         ns=np+npout
      endif



