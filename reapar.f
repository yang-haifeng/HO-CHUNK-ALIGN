      SUBROUTINE REAPAR(cflpar,iopar
     &     )

C***********************************************************************
C Subroutine REAPAR handles the reading of input parameters from the
C ".par" file, as well as elementary processing with those input
C parameters to generate arrays.
c
C Original version created by P.J.Flatau, Colorado State Univ.
C Modified by B.T.Draine, Princeton Univ. Obs.
C ported from program DDSCAT by m. wolff
c
C History:
C 95/01/17 (mjw): stripped out of DDSCAT and ported to BLOB
C 99/12/02 (mjw): stripped out of BLOB and ported to MCTHERM
c                 parameter OCCULT defined in calling routine (SETUP)
c 00/03/19 (mjw): added some Vger stuff (Vger position in units of rmax)
c 00/06/30 (baw): redo some spot input
c 00/07/27 (mjw): remove extraneous commas that drove lf95 crazy
c                   (occurrences of "write(*,*),something")
C
C***********************************************************************
C

      implicit none

C     .. Scalar Arguments ..
      character cflpar*(*)
      integer iopar


C     .. Local Scalars ..
      character cline*70,cmsgnm*70
      character*3 cbub,chole,climb,cpeel,cpoly,cstream
     1   ,cwrite,cspot,cveeg,cnspot
C     ..
C     .. Local Arrays ..
C     ..
C     .. External Subroutines ..
      external errmsg, wrimsg
C     ..
C     .. Intrinsic Functions ..
C     ..
C     .. Data statements ..
C     ..
c     .. Include statements ..

      include 'opacin.txt'
      include 'out.txt'
      include 'random.txt'
      include 'tabl.txt'
      include 'tts.txt'
      include 'spot.txt'
      include 'vger.txt'


C***********************************************************************

      open(unit=iopar,file=cflpar,status='old')

      READ(IOPAR,FMT=*,ERR=40)CLINE
      CALL WRIMSG(' ',CLINE)

C***********************************************************************
C ... Preliminaries
c
      READ(IOPAR,FMT=*,ERR=40)CLINE
      CALL WRIMSG(' ',CLINE)
      READ(IOPAR,FMT=*,ERR=40)np
      WRITE(CMSGNM,FMT='(i11,a)')np
     &     ,' = NP = Number of photons from central star'
      CALL WRIMSG('REAPAR',CMSGNM)
      READ(IOPAR,FMT=*,ERR=40)npout
      WRITE(CMSGNM,FMT='(i11,a)')npout
     &     ,' = NPOUT = Number of photons from outside star'
      CALL WRIMSG('REAPAR',CMSGNM)
      READ(IOPAR,FMT=*,ERR=40)i1
      WRITE(CMSGNM,FMT='(i10,a)')i1
     &     ,' = I1 = random number seed'
      CALL WRIMSG('REAPAR',CMSGNM)
      READ(IOPAR,FMT=*,ERR=40)iwrite
      WRITE(CMSGNM,FMT='(i9,a)')iwrite
     &     ,' = IWRITE = no. of photons between timing calls'
      CALL WRIMSG('REAPAR',CMSGNM)
c      READ(IOPAR,FMT=*,ERR=40)ctherm
c      if ( (ctherm.eq.'YES').OR.(ctherm.eq.'yes') ) then
c         itherm = 1
c         WRITE(CMSGNM,FMT='(a)')
c     &        ,'ITHERM=1, Calculating thermal emission + scattering'
c      else
c         itherm = 0
c         WRITE(CMSGNM,FMT='(a)')
c     &        ,'ITHERM=0, Calculating ONLY scattering'
c      endif
c      CALL WRIMSG('REAPAR',CMSGNM)
      READ(IOPAR,FMT=*,ERR=40)cpeel
      if ( (cpeel.eq.'YES').OR.(cpeel.eq.'yes') ) then
         ipeel = 1
         WRITE(CMSGNM,FMT='(a)')
     &        'IPEEL=1, Using peeling-off procedure'
      else
         ipeel = 0
         WRITE(CMSGNM,FMT='(a)')
     &        'IPEEL=0, Using standard MC (no peeling-off)'
      endif
      CALL WRIMSG('REAPAR',CMSGNM)
      READ(IOPAR,FMT=*,ERR=40)cveeg
      if ( (cveeg.eq.'YES').OR.(cveeg.eq.'yes') ) then
         iveeg = 1
         WRITE(CMSGNM,FMT='(a)')
     &        'IVEEG=1, calculating radiation field at point inside'
      else
         iveeg = 0
         WRITE(CMSGNM,FMT='(a)')
     &        'IVEEG=0, not calculating rad. field at point inside'
      endif
      CALL WRIMSG('REAPAR',CMSGNM)
      READ(IOPAR,FMT=*,ERR=40)cwrite
      if ( (cwrite.eq.'YES').OR.(cwrite.eq.'yes') ) then
         iwriteim = 1
         WRITE(CMSGNM,FMT='(a)')
     &        'IWRITEIM=1, writing out all images'
      else
         iwriteim = 0
         WRITE(CMSGNM,FMT='(a)')
     &        'IWRITEIM=0, only writing peeling-off images'
      endif
      CALL WRIMSG('REAPAR',CMSGNM)
      READ(IOPAR,FMT=*,ERR=40)climb
      if ( (climb.eq.'YES').OR.(climb.eq.'yes') ) then
         limb = 1
         WRITE(CMSGNM,FMT='(a)')
     &        'LIMB=1, Limb-darkening for central source'
      else
         limb = 0
         WRITE(CMSGNM,FMT='(a)')
     &        'LIMB=0, NO limb-darkening for central source'
      endif
      CALL WRIMSG('REAPAR',CMSGNM)

C***********************************************************************
C ... Central source properties
c
      READ(IOPAR,FMT=*,ERR=40)CLINE
      CALL WRIMSG(' ',CLINE)
      READ(IOPAR,FMT=*,ERR=40)rstar
      WRITE(CMSGNM,FMT='(1pe13.6,a)')rstar
     &     ,' = RSTAR = Stellar radius in Solar radii'
      CALL WRIMSG('REAPAR',CMSGNM)
      READ(IOPAR,FMT=*,ERR=40)tstar
      WRITE(CMSGNM,FMT='(1pe13.6,a)')tstar
     &     ,' = TSTAR = Blackbody temperature of central star (K)'
      CALL WRIMSG('REAPAR',CMSGNM)

      READ(IOPAR,FMT=*,ERR=40)massc
      WRITE(CMSGNM,FMT='(1pe13.6,a)')massc
     &     ,' = MASSC = Mass of central star (for TSC properties)'
      CALL WRIMSG('REAPAR',CMSGNM)

      READ(IOPAR,FMT=*,ERR=40) cspot
      if ( (cspot.eq.'YES').OR.(cspot.eq.'yes') ) then
         ispot = 1
         WRITE(CMSGNM,FMT='(a)')
     &        'ISPOT=1, emission from star+spot'
	   spotflag = 1   !require to run spotset first time
      else
         ispot = 0
         WRITE(CMSGNM,FMT='(a)')
     &        'ISPOT=0, uniform emission from star'
      endif
      CALL WRIMSG('REAPAR',CMSGNM)

      READ(IOPAR,FMT=*,ERR=40)tspot
      WRITE(CMSGNM,FMT='(1pe13.6,a)')tspot
     &     ,' = TSTAR = Blackbody temperature of spot (K)'
      CALL WRIMSG('REAPAR',CMSGNM)

      READ(IOPAR,FMT=*,ERR=40)thspot
      WRITE(CMSGNM,FMT='(1pe13.6,a)')thspot
     &     ,' = THSPOT = radius of spot in degrees'
      CALL WRIMSG('REAPAR',CMSGNM)

      READ(IOPAR,FMT=*,ERR=40)thspotin
      WRITE(CMSGNM,FMT='(1pe13.6,a)')thspotin
     &     ,' = THSPOTIN = inner radius of spot (for ring)'
      CALL WRIMSG('REAPAR',CMSGNM)

      READ(IOPAR,FMT=*,ERR=40)spotlat
      WRITE(CMSGNM,FMT='(1pe13.6,a)')spotlat
     &     ,' = SPOTLAT = latitude of spot (degrees)'
      CALL WRIMSG('REAPAR',CMSGNM)

      READ(IOPAR,FMT=*,ERR=40)spotlon
      WRITE(CMSGNM,FMT='(1pe13.6,a)')spotlon
     &     ,' = SPOTLON = longitude of spot (degrees)'
      CALL WRIMSG('REAPAR',CMSGNM)

      READ(IOPAR,FMT=*,ERR=40) cnspot
      if ( (cnspot.eq.'ONE').OR.(cnspot.eq.'one') ) then
         nspot = 1
         WRITE(CMSGNM,FMT='(a)')
     &        'NSPOT=1, one spot'
      else
         nspot = 2
         WRITE(CMSGNM,FMT='(a)')
     &        'NSPOT=2, two spots'
      endif
      CALL WRIMSG('REAPAR',CMSGNM)

C***********************************************************************
C ... Disk properties
c
      READ(IOPAR,FMT=*,ERR=40)CLINE
      CALL WRIMSG(' ',CLINE)
      READ(IOPAR,FMT=*,ERR=40)massd
      WRITE(CMSGNM,FMT='(1pe13.6,a)')massd
     &     ,' = MASSD = Disk mass in solar masses'
      CALL WRIMSG('REAPAR',CMSGNM)
      READ(IOPAR,FMT=*,ERR=40)rmaxd
      WRITE(CMSGNM,FMT='(1pe13.6,a)')rmaxd
     &     ,' = RMAXD = Maximum disk radius'
      CALL WRIMSG('REAPAR',CMSGNM)
      READ(IOPAR,FMT=*,ERR=40)rmind
      WRITE(CMSGNM,FMT='(1pe13.6,a)')rmind
     &     ,' = RMIND = Maximum disk radius'
      CALL WRIMSG('REAPAR',CMSGNM)
      READ(IOPAR,FMT=*,ERR=40)rddust
      WRITE(CMSGNM,FMT='(1pe13.6,a)')rddust
     &     ,' = RDDUST = Minimum disk opacity radius in Rstar'
      CALL WRIMSG('REAPAR',CMSGNM)
      READ(IOPAR,FMT=*,ERR=40)zmin
      WRITE(CMSGNM,FMT='(1pe13.6,a)')zmin
     &     ,' = ZMIN = Scale height of disk at Rstar'
      CALL WRIMSG('REAPAR',CMSGNM)
      READ(IOPAR,FMT=*,ERR=40)a
      WRITE(CMSGNM,FMT='(1pe13.6,a)')a
     &     ,' = A = Disk density exponent (~ r^(-a))'
      CALL WRIMSG('REAPAR',CMSGNM)
      READ(IOPAR,FMT=*,ERR=40)b
      WRITE(CMSGNM,FMT='(1pe13.6,a)')b
     &     ,' = B = Disk scale height exponent (~ r^(b))'
      CALL WRIMSG('REAPAR',CMSGNM)

C***********************************************************************
C ... Envelope properties
c
      READ(IOPAR,FMT=*,ERR=40)CLINE
      CALL WRIMSG(' ',CLINE)

      READ(IOPAR,FMT=*,ERR=40)rmax
      WRITE(CMSGNM,FMT='(1pe13.6,a)')rmax
     &     ,' = RMAX = Maximum envelope radius in AU'
      CALL WRIMSG('REAPAR',CMSGNM)
      READ(IOPAR,FMT=*,ERR=40)rmine
      WRITE(CMSGNM,FMT='(1pe13.6,a)')rmine
     &     ,' = RMINE = Minimum envelope radius in AU'
      CALL WRIMSG('REAPAR',CMSGNM)
      READ(IOPAR,FMT=*,ERR=40)rate
      WRITE(CMSGNM,FMT='(1pe13.6,a)')rate
     &     ,' = RATE = Mass infall rate for TSC env.'//
     &     '(M_O/yr)'
      CALL WRIMSG('REAPAR',CMSGNM)
      READ(IOPAR,FMT=*,ERR=40)rc
      WRITE(CMSGNM,FMT='(1pe13.6,a)')rc
     &     ,' = RC = centrifugal radius for TSC env.'
      CALL WRIMSG('REAPAR',CMSGNM)
      READ(IOPAR,FMT=*,ERR=40)chole
      READ(IOPAR,FMT=*,ERR=40)thet1
      READ(IOPAR,FMT=*,ERR=40)cpoly
      READ(IOPAR,FMT=*,ERR=40)thet2
      READ(IOPAR,FMT=*,ERR=40)ex1
      READ(IOPAR,FMT=*,ERR=40)ex2
      READ(IOPAR,FMT=*,ERR=40)z01
      READ(IOPAR,FMT=*,ERR=40)z02
      READ(IOPAR,FMT=*,ERR=40)exf
      READ(IOPAR,FMT=*,ERR=40)rhoconst1
      READ(IOPAR,FMT=*,ERR=40)rhoconst2
      READ(IOPAR,FMT=*,ERR=40)rhoamb
      READ(IOPAR,FMT=*,ERR=40)zflowmin
      READ(IOPAR,FMT=*,ERR=40)cstream
      READ(IOPAR,FMT=*,ERR=40)rchole      
      READ(IOPAR,FMT=*,ERR=40)thetmu0

      if ( (chole.eq.'YES').OR.(chole.eq.'yes') ) then
         ihole = 1
         WRITE(CMSGNM,FMT='(a)')
     &        'IHOLE=1, Carving cavity out'
         CALL WRIMSG('REAPAR',CMSGNM)
         WRITE(CMSGNM,FMT='(f6.1,a)')thet1
     &        ,' = THET1 = Opening angle of inner cavity wall'
         CALL WRIMSG('REAPAR',CMSGNM)
         if ( (cpoly.eq.'YES').OR.(cpoly.eq.'yes') ) then
            ipoly = 1
            WRITE(CMSGNM,FMT='(a)')
     &           'IPOLY=1, Polynomial-shaped cavity'
            CALL WRIMSG('REAPAR',CMSGNM)
            WRITE(CMSGNM,FMT='(f6.1,a)')thet2
     &           ,' = THET2 = Opening angle of outer cavity wall'
         else
            ipoly = 0
            WRITE(CMSGNM,FMT='(a)')
     &           'IPOLY=0, NON-polynomial-shaped cavity'
         endif
         CALL WRIMSG('REAPAR',CMSGNM)
         WRITE(CMSGNM,FMT='(f7.2,a)')ex1
     &        ,' = EX1 = Inner Cavity wall thickness'
         CALL WRIMSG('REAPAR',CMSGNM)
         WRITE(CMSGNM,FMT='(f7.2,a)')ex2
     &        ,' = EX2 = Outer Cavity wall thickness'
         CALL WRIMSG('REAPAR',CMSGNM)
         WRITE(CMSGNM,FMT='(1pe11.4,a)')z01
     &        ,' = Z01 = Height of inner wall at w=0'
         CALL WRIMSG('REAPAR',CMSGNM)
         WRITE(CMSGNM,FMT='(1pe11.4,a)')z02
     &        ,' = Z02 = Height of outer wall at w=0'
         CALL WRIMSG('REAPAR',CMSGNM)
         WRITE(CMSGNM,FMT='(f7.2,a)')exf
     &        ,' = EXF = exponent for cavity density power-law'
         CALL WRIMSG('REAPAR',CMSGNM)
         WRITE(CMSGNM,FMT='(1pe11.4,a)')rhoconst1
     &        ,' = RHOCONST1 = Coefficient for inner cavity density'
         CALL WRIMSG('REAPAR',CMSGNM)
         WRITE(CMSGNM,FMT='(1pe11.4,a)')rhoconst2
     &        ,' = RHOCONST2 = Coefficient for outer cavity density'
         CALL WRIMSG('REAPAR',CMSGNM)
         WRITE(CMSGNM,FMT='(1pe11.4,a)')rhoamb
     &        ,' = ambient cloud density (gm/cm^3)'
         CALL WRIMSG('REAPAR',CMSGNM)
         WRITE(CMSGNM,FMT='(1pe11.4,a)')zflowmin
     &        ,' = minimum radius of outflow (AU)'
         CALL WRIMSG('REAPAR',CMSGNM)
         if ( (cstream.eq.'YES').OR.(cstream.eq.'yes') ) then
            istream = 1
            WRITE(CMSGNM,FMT='(a)')
     &           'ISTREAM=1, streamline hole'
            CALL WRIMSG('REAPAR',CMSGNM)
            WRITE(CMSGNM,FMT='(1pe11.4,a)')rchole
     &           ,' = RCHOLE = Streamline hole size in AU'
            CALL WRIMSG('REAPAR',CMSGNM)
            WRITE(CMSGNM,FMT='(f6.2,a)')thetmu0
     &           ,' = THETMU0 = Opening angle of streamline hole (deg)'
         
         else
            istream = 0
            WRITE(CMSGNM,FMT='(a)')
     &           'ISTREAM=0'
         endif
      else
         ihole = 0
         WRITE(CMSGNM,FMT='(a)')
     &        'IHOLE=0, NO cavity'
      endif
      CALL WRIMSG('REAPAR',CMSGNM)
      READ(IOPAR,FMT=*,ERR=40)cbub
      READ(IOPAR,FMT=*,ERR=40)nbub
      READ(IOPAR,FMT=*,ERR=40)zbub1
      READ(IOPAR,FMT=*,ERR=40)zbub2
      READ(IOPAR,FMT=*,ERR=40)buboa
      if ( (cbub.eq.'YES').OR.(cbub.eq.'yes') ) then
         ibub = 1
         WRITE(CMSGNM,FMT='(a)')
     &        'IBUB=1, including bubble'
         CALL WRIMSG('REAPAR',CMSGNM)
         WRITE(CMSGNM,FMT='(f7.2,a)')nbub
     &        ,' = NBUB = Exponent describing bubble shape'
         CALL WRIMSG('REAPAR',CMSGNM)
         WRITE(CMSGNM,FMT='(1pe11.4,a)')zbub1
     &        ,' = ZBUB1 = height of bubble inner wall'
         CALL WRIMSG('REAPAR',CMSGNM)
         WRITE(CMSGNM,FMT='(1pe11.4,a)')zbub2
     &        ,' = ZBUB2 = height of bubble outer wall'
         CALL WRIMSG('REAPAR',CMSGNM)
         WRITE(CMSGNM,FMT='(1pe11.4,a)')buboa
     &        ,' = BUBOA = Opening angle of hole in bubble'
      else
         ibub = 0
         WRITE(CMSGNM,FMT='(a)')
     &        'IBUB=0, NO bubble'
      endif
      CALL WRIMSG('REAPAR',CMSGNM)

C***********************************************************************
C ... Vger stuff
c
c XVG0,YVG0,ZVG0 = location of VGER in units of RMAX
c XSRC0,YSRC0,ZSRC0 = location of outer illum. (RMAX)
c ISRC0,QSRC0,USRC0,VSRC0 = normalized stokes for incident outside illum.
c
      READ(IOPAR,FMT=*,ERR=40)CLINE
      CALL WRIMSG(' ',CLINE)
      READ(IOPAR,FMT=*,ERR=40)xvg0,yvg0,zvg0
      READ(IOPAR,FMT=*,ERR=40)xsrc0,ysrc0,zsrc0
      READ(IOPAR,FMT=*,ERR=40)isrc0,qsrc0,usrc0,vsrc0
      if (CVEEG.eq.'YES') then
         write(CMSGNM,fmt='(3(1pe11.4,1x),a)') xvg0,yvg0,zvg0,
     &        ' - VGER location in units of RMAX'
         CALL WRIMSG('REAPAR',CMSGNM)
         write(CMSGNM,fmt='(3(1pe11.4,1x),a)') xsrc0,ysrc0,zsrc0,
     &        ' - src illum. in units of RMAX'
         CALL WRIMSG('REAPAR',CMSGNM)
         write(CMSGNM,fmt='(4(1pe11.4,1x),a)')isrc0,qsrc0,usrc0,
     &        vsrc0,
     &       ' - incident stokes'
         CALL WRIMSG('REAPAR',CMSGNM)
      endif
      
      
C***********************************************************************
C ... output data
c
      READ(IOPAR,FMT=*,ERR=40)CLINE
      CALL WRIMSG(' ',CLINE)
      READ(IOPAR,FMT=*,ERR=40)rmaxi
      WRITE(CMSGNM,FMT='(1pe13.6,a)')rmaxi
     &     ,' = RMAXI = Image half-size in AU'
      CALL WRIMSG('REAPAR',CMSGNM)
      READ(IOPAR,FMT=*,ERR=40)aperture2
      WRITE(CMSGNM,FMT='(3(f8.1,1x),a)')aperture2
     &     ,' = APERTURE(1-3) = radius of aper 1-3 in AU'
      CALL WRIMSG('REAPAR',CMSGNM)
      READ(IOPAR,FMT=*,ERR=40)thete
      WRITE(CMSGNM,FMT='(f7.2,a)')thete
     &     ,' = THETE = theta angle (deg) of high S/N image'
      CALL WRIMSG('REAPAR',CMSGNM)
      READ(IOPAR,FMT=*,ERR=40)phie
      WRITE(CMSGNM,FMT='(f7.2,a)')phie
     &     ,' = PHIE = phi angle (deg) of high S/N image'
      CALL WRIMSG('REAPAR',CMSGNM)
      READ(IOPAR,FMT=*,ERR=40)xymaxdens
      WRITE(CMSGNM,FMT='(f7.2,a)')xymaxdens
     &     ,' =XYMAXDENS = xmax for xydens.dat in units of rstar'
      CALL WRIMSG('REAPAR',CMSGNM)

      CLOSE (IOPAR)
      RETURN
   40 CONTINUE
      CALL ERRMSG('FATAL','REAPAR',' Error reading .par file')
      RETURN
      END




