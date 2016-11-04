      program tts

c     radiative transfer in a disk+envelope around a star.
c     height of disk, z, goes as z1*r**b.
c     density goes as c*r**-a.
c     opacity is dust with either rayleigh or H-G phase function.
c     output is stokes parameters as a function of angle
c
c history:
c 00/03/19 (mjw):  update version stamp for changes on 03/19.  Also,
c			change version name to VGER (from mctherm)
c

      integer iwav,MXWAV,nwav
	parameter (MXWAV=15)
      character cwav(MXWAV)*10,cfldust(MXWAV)*50,cvers*10,
     1   ctherm(MXWAV)*8

	data cvers/'20000814a'/

      open(unit=12,file='disk.dat',status='unknown')
	write(12,'(a,a)')cvers,' - VGER version'

c     input and setup
      call setup(cwav,cfldust,ctherm,nwav,MXWAV)

	do iwav=1,nwav

	   call initarr
      
	   call setup_wave(iwav,cfldust,ctherm,MXWAV)

c        radiative transfer
         call disk

c        convert stokes arrays to flux arrays
         call diskflux

         call output(iwav,cwav,ctherm,MXWAV)

	end do

      close(12)


      end
