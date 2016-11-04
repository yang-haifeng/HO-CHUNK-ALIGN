c     *********************************************************

      subroutine output(iwav,cwav,ctherm,MXWAV)

      implicit none 

      include 'tts.txt'
      include 'out.txt'
      include 'stokes.txt'
      include 'opacin.txt'
      include 'vger.txt'
      include 'grid.txt'

      integer iwav,MXWAV
      character cwav(MXWAV)*(*),ctherm(MXWAV)*(*)

      real dflux,scatave,cpusec,cpuhrs,cpumin,sec,time(2)
     1   ,fq,etime,fqp,fup,fvp,fi,fu,fv,rnorm,fp,totvg,dq,du,dv,di

      integer i,ia,imin,ihrs,ix,iy,it,ip,ir

      integer MXNAM
      parameter (MXNAM=20)
      character filnam(MXNAM)*25
	      
      call namer(1,cwav(iwav),filnam,MXNAM)

      if (iwriteim.eq.1) then
         do it=1,nmu
	      call namer(it,cwav(iwav),filnam,MXNAM)
               !make 7-element array for filnam for i,q,u,v
               open(unit=21,file=filnam(1),status='unknown')
               do iy=1,nx
                  write(21,*) (sngl(image(ix,iy,it)),ix=1,nx)
               end do
               close(21)
               open(unit=21,file=filnam(2),status='unknown')
               open(unit=22,file=filnam(3),status='unknown')
               open(unit=23,file=filnam(4),status='unknown')
               open(unit=24,file=filnam(5),status='unknown')
               do iy=1,nx
                  write(21,*) (sngl(imagei(ix,iy,it)),ix=1,nx)
                  write(22,*) (sngl(imageq(ix,iy,it)),ix=1,nx)
                  write(23,*) (sngl(imageu(ix,iy,it)),ix=1,nx)
                  write(24,*) (sngl(imagev(ix,iy,it)),ix=1,nx)
               end do
               close(21)
               close(22)
               close(23)
               close(24)
	   ! and imagei,q,u,v over nx
	   end do
	endif

      open(unit=21,file=filnam(8),status='unknown')
      do it=1,nmu
c         write(21,*) sngl(star(it))
         if (star(it).gt.0.d0) then
            dq=sngl(starq(it)/star(it))
            du=sngl(staru(it)/star(it))
            dv=sngl(starv(it)/star(it))
         else
            dq=0.
            du=0.
            dv=0.
         endif
         write(21,*) sngl(star(it)),dq,du,dv
      end do
      close(21)

      !name these for each wavelength     
      open(unit=11,file=filnam(6),status='unknown')
      open(unit=13,file=filnam(7),status='unknown')
      write(11,*) 
     1'cost,          i,          q/i,          u/i,          v/i'
      write(13,*) 
     1'cost,          i,          q/i,          u/i,          v/i'
      rnorm=float(np)/float(nmu-1)
      do ia=1,3
      do i=1,nmu
         if(spi(i,ia).eq.0.0d+00) then
            fqp=0.
            fup=0.
            fvp=0.
         else
            fqp=sngl(spq(i,ia)/spi(i,ia))
            fup=sngl(spu(i,ia)/spi(i,ia))
            fvp=sngl(spv(i,ia)/spi(i,ia))
         end if
         if(si(i,ia).eq.0.d0) then
            fq=0.
            fu=0.
            fv=0.
         else
            fq=sngl(sq(i,ia)/si(i,ia))
            fu=sngl(su(i,ia)/si(i,ia))
            fv=sngl(sv(i,ia)/si(i,ia))
         end if
         if(i.eq.1.or.i.eq.nmu) then
            fi=sngl(spi(i,ia))/rnorm*2
         else
            fi=sngl(spi(i,ia))/rnorm
         end if
         write(11,902) sngl(u(i)),sngl(si(i,ia)),fq,fu,fv
         write(13,902) sngl(u(i)),fi,fqp,fup,fvp
      end do
      end do
	   close(11)
	   close(13)
902   format(F10.4,4(2X,1pe12.5))

      open(unit=11,file=filnam(9),status='unknown')
	write(11,*) 'standard deviations of i,q,u,v'
	write(11,*) 'normalize to i for q,u,v to compare directly with'
	write(11,*) '    flux file (assumes i has no error) '
      write(11,*) 
     1'  cost,       sig_i,         sig_q/i,      sig_u/i,      sig_v/i'
      do ia=1,3
      do i=1,nmu
         if(nums(i,ia).le.1) then
            fi=0.
            fq=0.
            fu=0.
            fv=0.
            fp=0.
         else
            fi=sngl(si2(i,ia)*si(i,ia))
            fq=sngl(sq2(i,ia))
            fu=sngl(su2(i,ia))
            fv=sngl(sv2(i,ia))
            fp=1./sqrt(float(nums(i,ia)))*sngl(si(i,ia))
	   endif
         write(11,902) sngl(u(i)),fi,fq,fu,fv
      end do
	enddo
	close(11)

c     flux
      flux=flux+aflux
c     flux which hits disk
      dflux=sngl(tot/flux)
c     scattered flux
      if(nscat.gt.0) then
         scatave=nscat/sflux
      else
         scatave=0.
      end if
      sflux=sflux/flux
c     absorbed flux
      aflux=aflux/flux
c     cpu time
      cpusec=etime(time)
      cpuhrs=cpusec/3600.
      ihrs=int(cpuhrs)
      cpumin=(cpuhrs-ihrs)*60.
      imin=int(cpumin)
      sec=(cpumin-imin)*60.

      write(12,*) 'rmax,rmaxd,zmax',rmax,rmaxd,zmax
      write(12,*) 'rmaxi,rmind',rmaxi,rmind
      write(12,*) 'photons   ',np
      write(12,*) 'taur ',taur
      write(12,*) 'total flux', flux
      write(12,*) 'flux which hits disk+envelope',dflux,dflux*flux
      write(12,*) 'scattered flux',sflux,sflux*flux
      write(12,*) 'ave number of scatters in this flux',scatave
      write(12,*) 'flux which gets absorbed',aflux,aflux*flux
      write(12,*) 'cputime: ',ihrs,' hrs ',imin,' min ',sec,' sec'

      it=int(thete*180.001/pi+1)
      call namer(it,cwav(iwav),filnam,MXNAM)

      if (ipeel.eq.1) then
         open(unit=13,file='e'//filnam(6),status='unknown')
	   write(13,*) 'flux and polarization in the 3 apertures '
         write(13,*)
     1'  I              Q/I              U/I               V/I'
         write(13,*) 'total flux'
         do i=1,3
            if(ti(i).eq.0.d0) then
               fq=0.
               fu=0.
               fv=0.
            else
               fq=sngl(tq(i)/ti(i))
               fu=sngl(tu(i)/ti(i))
               fv=sngl(tv(i)/ti(i))
            end if
            write(13,*) sngl(ti(i)),fq,fu,fv
         enddo
         write(13,*) 'scattered flux only'
         do i=1,3
            if(tis(i).eq.0.d0) then
               fq=0.
               fu=0.
               fv=0.
            else
               fq=sngl(tqs(i)/tis(i))
               fu=sngl(tus(i)/tis(i))
               fv=sngl(tvs(i)/tis(i))
            end if
            write(13,*) sngl(tis(i)),fq,fu,fv
         enddo
         write(13,*) 'total - scatter flux ( = direct)'
          do i=1,3
             di=ti(i)-tis(i)
             dq=tq(i)*ti(i)-tqs(i)*tis(i)
             du=tu(i)*ti(i)-tus(i)*tis(i)
             dv=tv(i)*ti(i)-tvs(i)*tis(i)
             if (di.gt.0.) then
                dq=dq/di
                du=du/di
                dv=dv/di
             else
                dq=0.
                du=0.
                dv=0.
             endif
             write(13,*) di,dq,du,dv
          enddo
          write(13,*) 'direct flux from star only'
          if (tid.gt.0) then
             tqd=tqd/tid
             tud=tud/tid
             tvd=tvd/tid
          else
             tqd=0
             tud=0
             tvd=0
          endif
          write(13,*) tid,tqd,tud,tvd
          write(13,*) ' standard deviations in the 3 apertures '
         write(13,*)
     1' sigI          sig(Q/I)     sig(U/I)     sig(V/I)     Poisson(I)'
         do ia=1,3
            if(numt(ia).le.1) then
               fi=0.
               fq=0.
               fu=0.
               fv=0.
               fp=0.
            else
               fi=sngl(ti2(ia)*ti(ia))
               fq=sngl(tq2(ia))
               fu=sngl(tu2(ia))
               fv=sngl(tv2(ia))
               fp=1./sqrt(float(numt(ia)))*sngl(ti(ia))
	      endif
            write(13,*) fi,fq,fu,fv,fp
         end do
         close(13)

         open(unit=21,file='e'//filnam(2),status='unknown')
         open(unit=22,file='e'//filnam(3),status='unknown')
         open(unit=23,file='e'//filnam(4),status='unknown')
         open(unit=24,file='e'//filnam(5),status='unknown')
            do iy=1,nxhst
               write(21,*) (sngl(tihst(ix,iy)),ix=1,nxhst)
               write(22,*) (sngl(tqhst(ix,iy)),ix=1,nxhst)
               write(23,*) (sngl(tuhst(ix,iy)),ix=1,nxhst)
               write(24,*) (sngl(tvhst(ix,iy)),ix=1,nxhst)
            end do
         close(21)
         close(22)
         close(23)
         close(24)
	endif

c     write out average stokes intensities in grid 
      open(unit=15,file='aveI.unf',status='unknown',form='unformatted')
      write(15) nrg-1,ntg-1,npg-1
      write(15) (ravearr(ir)/autors,ir=1,nrg-1)
      write(15) (thetavearr(it),it=1,ntg-1)
      write(15) (phiavearr(ip),ip=1,npg-1)
      write(15) (((pathI(ir,it,ip),ir=1,nrg-1),it=1,ntg-1)
     1     ,ip=1,npg-1)
      close(15)

c     write out average stokes intensities in grid
      open(unit=15,file='aveQ.unf',status='unknown',form='unformatted')
      write(15) nrg-1,ntg-1,npg-1
      write(15) (ravearr(ir)/autors,ir=1,nrg-1)
      write(15) (thetavearr(it),it=1,ntg-1)
      write(15) (phiavearr(ip),ip=1,npg-1)
      write(15) (((pathQ(ir,it,ip),ir=1,nrg-1),it=1,ntg-1)
     1     ,ip=1,npg-1)
      close(15)

c     write out average stokes intensities in grid
      open(unit=15,file='aveU.unf',status='unknown',form='unformatted')
      write(15) nrg-1,ntg-1,npg-1
      write(15) (ravearr(ir)/autors,ir=1,nrg-1)
      write(15) (thetavearr(it),it=1,ntg-1)
      write(15) (phiavearr(ip),ip=1,npg-1)
      write(15) (((pathU(ir,it,ip),ir=1,nrg-1),it=1,ntg-1)
     1     ,ip=1,npg-1)
      close(15)

c     write out average stokes intensities in grid
      open(unit=15,file='aveV.unf',status='unknown',form='unformatted')
      write(15) nrg-1,ntg-1,npg-1
      write(15) (ravearr(ir)/autors,ir=1,nrg-1)
      write(15) (thetavearr(it),it=1,ntg-1)
      write(15) (phiavearr(ip),ip=1,npg-1)
      write(15) (((pathV(ir,it,ip),ir=1,nrg-1),it=1,ntg-1)
     1     ,ip=1,npg-1)
      close(15)

      if (iveeg.eq.1) then
         open(unit=21,file=filnam(10),status='unknown')
	   write(21,*) '      F         Q/F         U/F           V/F'
	   totvg=ivgflux
	   write(21,*) (totvg),sngl(qvgflux/totvg),
     1   sngl(uvgflux/totvg),sngl(vvgflux/totvg)
	   write(21,*) '    errF       errQ/F      err(U/F)      err(V/F)'
	   write(21,*) sngl(ivgerr)*(totvg),sngl(qvgerr),
     1   sngl(uvgerr),sngl(vvgerr)
	   close(21)
c        if you really want vger images, have to assign them names
c        in namer.  for now, just assume doing one wavelength.
         open(unit=21,file=filnam(11),status='unknown')
         open(unit=22,file=filnam(12),status='unknown')
         open(unit=23,file=filnam(13),status='unknown')
         open(unit=24,file=filnam(14),status='unknown')
	   do it=1,ncvg
	     write(21,*) (ivg(it,ip),ip=1,npvg)
	     write(22,*) (qvg(it,ip),ip=1,npvg)
	     write(23,*) (uvg(it,ip),ip=1,npvg)
	     write(24,*) (vvg(it,ip),ip=1,npvg)
	   end do
	   close(21)
	   close(22)
	   close(23)
	   close(24)
	endif
            
      return
      end
