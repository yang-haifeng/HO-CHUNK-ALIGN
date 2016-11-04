      subroutine phiface(iface,ifound,t,A1,B1,C1,D1,A2,B2,C2,D2,
     &     ux,uy,uz,x1,y1,z1)

c
c Computes the shortest distance (t) to constant phi (azimuthal angle)
c surface set of (phi1,phi2) from the point (x1,y1,z1)
c along the line defined by the directional cosines (ux,uy,uz).
c Intersections with negative "t" are NOT selected (i.e., are
c OPPOSITE of desired direction of propagation).
c
c NOTE:  constant phi surfaces are planes that contain the z-axis.
c            Ax + By +Cz + D = 0, where the coefficients are
c            pre-computed.
c
c Output values:
c     iface = which constant radius surface is closest
c           = 0 for phi1, 1 for phi2, and -1 for NEITHER
c     ifound = flag for non-negative 
c            = .T. for positive "t" found, .F. NOT.
c     t = smallest distance for surface intersection
c         (negative value indicated desired interesction does NOT exist...
c          but user should test for this case using "ifound".
c
c Input values:
c     ux,uy,uz = directional cosines along x,y,z directions
c     x1,y1,z1 = point from which distance is calculated
c     A1...D1,A2...D2 = plane equation coefficients for two constant surfaces
c
c m.j. wolff/b.a. whitney, 2000/02/04
c history:
c
c 00/03/16 (mjw):  Using double precision (real*8).
c                  Redefine array IND as integer (was real).
c                  Identify only roots greater than 0.
c
c
c

      implicit none

c ... scalar arguments
      integer iface
      real*8 A1,A2,B1,B2,C1,C2,D1,D2,ux,uy,uz,x1,y1,z1,t
      logical ifound

c ... local scalars
      integer i
      real*8 denom1,denom2,rootmax,rootmin

c ... local arrays
      integer ind(2)
      real*8 root(2)

c ... IND array indcates which "face" is associated with which roots
      ind(1) = 0
      ind(2) = 1

      denom1 = (A1*ux + B1*uy + C1*uz)
      if (denom1.eq.0.d0) then
         root(1) = -2.d0
      else
         root(1) = -(A1*x1 + B1*y1 + C1*z1 + D1) / denom1
      endif

      denom2 = (A2*ux + B2*uy + C2*uz)
      if (denom2.eq.0.d0) then
         root(2) = -2.d0
      else
         root(2) = -(A2*x1 + B2*y1 + C2*z1 + D2) / denom2
      endif

c      print*,denom1,denom2
c       print*,'phiface: ',root1,root2

      rootmax = max(root(1),root(2))
      rootmin = min(root(1),root(2))

c both roots .le. 0
      if (rootmax.le.0.d0) then
         ifound = .false.
         t = -1.
c one root .lt. 0
      elseif (rootmin.lt.0.d0) then
         ifound = .true.
         t = rootmax
c both roots ge 0
      else
         ifound = .true.
         t = rootmin
      endif

c ... find which index "t" originated
      iface = -1
      do i=1,2
         if (root(i).eq.t) iface=ind(i)
      end do

      return
      end




