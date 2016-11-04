      subroutine radface(iface,ifound,t,ux,uy,uz,x1,y1,z1,R2_1,R2_2)
c
c Computes the shortest distance (t) to constant radius surface from set of
c radii=(R1,R2) from the point (x1,y1,z1) along the line defined by the
c directional cosines (ux,uy,uz).  Intersections with negative "t" are
c NOT selected (i.e., are OPPOSITE of desired direction of propagation).
c
c NOTE:  spherical surfaces are assumed to be centered at origin -- 
c            x**2 + y**2 + z**2 = R**2
c
c Output values:
c     iface = which constant radius surface is closest
c           = 0 for R1, 1 for R2, and -1 for NEITHER
c     ifound = flag for non-negative 
c            = .T. for positive "t" found, .F. NOT.
c     t = smallest distance for surface intersection
c         (negative value indicated desired interesction does NOT exist...
c          but user should test for this case using "ifound".
c
c Input values:
c     ux,uy,uz = directional cosines along x,y,z directions
c     x1,y1,z1 = point from which distance is calculated
c     R1,R2 = radii of two constant surfaces
c
c m.j. wolff/b.a. whitney, 2000/02/04
c history:
c 00/02/28 (mjw):  disallow "0" as a valid distance
c 00/03/16 (mjw):  now passing R**2 instead of R
c                  using double precision (real*8)
c

      implicit none

c ... scalar arguments
      integer iface
      real*8 ux,uy,uz,x1,y1,z1,t,R2_1,R2_2
      logical ifound

c ... local arrays
      real*8 ind(4),root(4),posroot(4)

c ... local scalars
      integer i,npos
      real*8 bb,cc,descr

c ... IND array indcates which "face" is associated with which roots
      ind(1) = 0
      ind(2) = 0
      ind(3) = 1
      ind(4) = 1

      bb = 2.d0*(x1*ux + y1*uy + uz*z1)
c  simplify since aa = 1.
      cc = x1*x1 + y1*y1 + z1*z1 - R2_1
      descr = bb*bb - 4.d0*cc

c      print*,'descr = ',descr,bb,cc
      if (descr.lt.0.d0) then
         root(1) = -999.d0
         root(2) = -999.d0
      else
         descr = dsqrt(descr)
         root(1) = (-bb + descr) / 2.d0
         root(2) = (-bb - descr) / 2.d0
      endif

      cc = x1*x1 + y1*y1 + z1*z1 - R2_2
      descr = bb*bb - 4.d0*cc
c      print*,'descr = ',descr,bb,cc
      if (descr.lt.0.) then
         root(3) = -999.d0
         root(4) = -999.d0
      else
         descr = dsqrt(descr)
         root(3) = (-bb + descr) / 2.d0
         root(4) = (-bb - descr) / 2.d0
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
c      print*,'posroot: ',posroot

c all roots .lt. 0
      if (npos.eq.0) then
         ifound = .false.
         t = -1.
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




