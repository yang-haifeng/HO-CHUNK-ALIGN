c     *******************************************************

      subroutine setup(cwav,cfldust,ctherm,nwav,MXWAV)

      implicit none

      integer MXWAV,nwav
      character cwav(MXWAV)*(*),cfldust(MXWAV)*(*),ctherm(MXWAV)*(*)
	

      include 'tab.txt'
      include 'tts.txt'
      include 'stokes.txt'
      include 'tabl.txt'
      include 'opacin.txt'
      include 'out.txt'
c      include 'mrncoef.txt'
      include 'random.txt'

c      character cfile*50
      character cflpar*50
      integer iopar,i

c     data
      data rsol /6.9598d10/
      data msol /1.989d33/
      data rmin /1.d0/
      data aflux,flux,sflux /3*0.d0/
      data cflpar /'mctherm.par'/,iopar /11/

c     read in variables

      call reapar(cflpar,iopar)
c occult not in current documentation
      occult = 1

c     convert units from AU to rstar.
      autors=214.94d0/rstar
      rmax=rmax*autors
      rmaxd=rmaxd*autors
      rmaxi=rmaxi*autors      !image size
      rchole=rchole*autors
	rc=rc*autors
      zbub1=zbub1*autors
      zbub2=zbub2*autors
c     note: aperture in radius, not diameter
      aperture2(1)=aperture2(1)*autors
      aperture2(2)=aperture2(2)*autors
      aperture2(3)=aperture2(3)*autors
      write(12,*) 'aperture sizes in Rstar',(aperture2(i),i=1,3)
      aperture2(1)=aperture2(1)**2
      aperture2(2)=aperture2(2)**2
      aperture2(3)=aperture2(3)**2
      z01=z01*autors
      z02=z02*autors
      zflowmin=zflowmin*autors
      npsav=np+npout
      print*,'npsav',npsav
      
      idust=2
c     nelems=6
      
      write(12,*) 'nx,nxhst ',nx,nxhst
      write(12,*) 'envelope inner dust radius, rmine ',rmine
      write(12,*) 'disk inner radius, rmind ',rmind
      write(12,*) 'disk inner dust radius, rddust ',rddust
      
      call waveset(cwav,cfldust,ctherm,nwav,MXWAV)
                                ! read in from wavein.dat

c     checks
      if (ihole.eq.1.and.istream.eq.1.and.ipoly.eq.1) then
         write(6,*) 'ERROR, cant have both a streamline and'
         write(6,*) 'polynomial hole.  pick one (istream,  ipoly)'
         write(6,*) ' '
         go to 666
      end if
      go to 777
 666  stop
 777  continue

      write(12,*) 'rmax,rmaxd,rmaxi,rmind,zmin,thetmu0'
      write(12,*) rmax,rmaxd,rmaxi,rmind,zmin,thetmu0
      
      pi=4.0d0*datan(1.0d0)
      r2p=2.0d0*pi
      write(6,*) 'pi,r2p',pi,r2p
      thete=r2p*thete/360.d0
      phie=r2p*phie/360.d0
      sinte=dsin(thete)
      coste=dcos(thete)
      sinpe=dsin(phie)
      cospe=dcos(phie)
c      g2=g**2
c      if(g.eq.0.) isot=1
      hit=0.d0
      htot=0.d0

c     limb darkening = 0 for constant, 1 for eddington.
c     make new table if you want something else
      if (limb.eq.1) call table

c     calculate density based on input mass
c     convert r to solar units.

      z1=zmin
      if(rmaxd.gt.0) then
         zmax=z1*rmaxd**b
      else
         zmax=0.d0
      end if


	return
	end
