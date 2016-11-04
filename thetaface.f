      subroutine thetaface(iface,ifound,t,ux,uy,uz,x1,y1,z1,tan_th2_1,
     &     tan_th2_2)
c
c Computes the shortest distance (t) to constant theta (polar angle)
c surface set of (tan(theta1),tan(theta2)) from the point (x1,y1,z1)
c along the line defined by the directional cosines (ux,uy,uz).
c Intersections with negative "t" are NOT selected (i.e., are
c OPPOSITE of desired direction of propagation).
c
c NOTE:  constant theta surfaces are assumed to be centered at origin
c        with z-axis defining theta=0.  surfaces are cylindrical cones.
c            (x/a)**2 + (y/a)**2 = (z/c)**2, tan(theta) = a/c
c     ->  x**2 + y**2 - (z*tan(theta))**2 = 0.
c
c
c Output values:
c     iface = which constant radius surface is closest
c           = 0 for theta1, 1 for theta2, and -1 for NEITHER
c     ifound = flag for non-negative 
c            = .T. for positive "t" found, .F. NOT.
c     t = smallest distance for surface intersection
c         (negative value indicated desired interesction does NOT exist...
c          but user should test for this case using "ifound".
c
c Input values:
c     ux,uy,uz = directional cosines along x,y,z directions
c     x1,y1,z1 = point from which distance is calculated
c     tan_th2_1,_2 = (tan(theta))**2 for two constant surfaces
c        NOTE: ***** it is necessary to flag tan_th2 array using a -1.
c		to indicate theta=90.
c	 Program ASSUMES that theta=90 can occur in ONLY ONE of the
c	 tan_th2 values.
c
c m.j. wolff/b.a. whitney, 2000/02/04
c history:
c 00/02/09 (mjw):  add case to trap theta=90 as indicated by a negative
c                     value of tan_th2.  phiface is called to find the
c		      distance to the theta=90 (xy-plane) surface.  THIS
c	 	      WILL SLOW DOWN ALL theta<>90 calls as well.
c 00/03/16 (mjw):  Using double precision (real*8).
c                  Define array IND as integer (was accidentally real).
c                  Identify only positive roots (t.gt.0)
c                
c


      implicit none

c ... scalar arguments
      integer iface
      real*8 ux,uy,uz,x1,y1,z1,t,tan_th2_1,tan_th2_2
      logical ifound

c ... local scalars
      integer i,npos
      real*8 aa,bb,cc,descr,tmp1,tmp2,tmp3,tan_th2a,tan_th2b
      
c ... local arrays
      integer ind(4)
      real*8 root(4),posroot(4)

      tmp1 = x1*ux + y1*uy
      tmp2 = ux*ux + uy*uy
      tmp3 = x1*x1 + y1*y1
      
c
c put minimum value of tan_th2 in 'a' variable, and switch the order
c of the "face" offsets if necessary.  This makes it easier to test
c for possible occurance of theta=90 (as flagged by -1. value)
c
      if (tan_th2_2.lt.tan_th2_1) then
         tan_th2a = tan_th2_2
         tan_th2b = tan_th2_1
c ... IND array indcates which "face" is associated with which roots
         ind(1) = 1
         ind(2) = 1
         ind(3) = 0
         ind(4) = 0
      else
         tan_th2a = tan_th2_1
         tan_th2b = tan_th2_2
         ind(1) = 0
         ind(2) = 0
         ind(3) = 1
         ind(4) = 1
      endif

c if tan_th2a lt 0, then theta=90 surface is now xy-plane.
      if(tan_th2a.lt.0.d0) then
c          print*,'trapping theta=90'
          call phiface(iface,ifound,t,0.d0,0.d0,1.d0,0.d0,
     &         0.d0,0.d0,1.d0,0.d0,
     &         ux,uy,uz,x1,y1,z1)
          root(1) = t
          root(2) = -999.d0
      else 
         bb = 2.d0*(tmp1 - tan_th2a*uz*z1)
         aa = tmp2 - tan_th2a*uz*uz
         cc = tmp3 - z1*z1*tan_th2a
         descr = bb*bb - 4.d0*aa*cc
         if (descr.lt.0.d0) then
            root(1) = -999.d0
            root(2) = -999.d0
         else
            descr = dsqrt(descr)
            root(1) = (-bb + descr) / (2.d0*aa)
            root(2) = (-bb - descr) / (2.d0*aa)
         endif
      endif

c ASSUME tan_th2b CAN NEVER be .lt.0
      bb = 2.d0*(tmp1 - tan_th2b*uz*z1)
      aa = tmp2 - tan_th2b*uz*uz
      cc = tmp3 - z1*z1*tan_th2b
      descr = bb*bb - 4.d0*aa*cc
      if (descr.lt.0.d0) then
         root(3) = -999.d0
         root(4) = -999.d0
      else
         descr = dsqrt(descr)
         root(3) = (-bb + descr) / (2.d0*aa)
         root(4) = (-bb - descr) / (2.d0*aa)
      endif

c ... find only positive roots
      npos = 0
      do i=1,4
         if (root(i).gt.0.d0) then
            npos = npos + 1
            posroot(npos) = root(i)
         endif
      end do

c      print*,'root: ',root

c all roots .lt. 0
      if (npos.eq.0) then
         ifound = .false.
         t = -1.d0
         iface = -1
c at least one positive
      else
         ifound = .true.
         if (npos.eq.1) then
            t = posroot(1)
         elseif (npos.eq.2) then
            t = min(posroot(1),posroot(2))
         elseif (npos.eq.3) then
            t = min(posroot(1),posroot(2),posroot(3))
         else
            t = min(posroot(1),posroot(2),posroot(3),posroot(4))
         endif
c ... find which index "t" originated
         do i=1,4
            if (root(i).eq.t) iface=ind(i)
         end do
      endif

      return
      end
