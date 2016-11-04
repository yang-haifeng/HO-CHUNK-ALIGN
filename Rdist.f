      subroutine Rdist(ifound,t,ux,uy,uz,x1,y1,z1,R)
c
c Computes the shortest distance (t) to constant radius surface R
c from the point (x1,y1,z1) along the line defined by the
c directional cosines (ux,uy,uz).  Intersections with negative "t" are
c NOT selected (i.e., are OPPOSITE of desired direction of propagation).
c
c NOTE:  spherical surfaces are assumed to be centered at origin -- 
c            x**2 + y**2 + z**2 = R**2
c
c Output values:
c     ifound = flag for non-negative 
c            = .T. for positive "t" found, .F. NOT.
c     t = smallest distance for surface intersection
c         (negative value indicated desired interesction does NOT exist...
c          but user should test for this case using "ifound".
c
c Input values:
c     ux,uy,uz = directional cosines along x,y,z directions
c     x1,y1,z1 = point from which distance is calculated
c     R = radius of constant surface
c
c m.j. wolff/b.a. whitney, 2000/02/04
c history:
c
c

      implicit none

c ... scalar arguments
      real*8 ux,uy,uz,x1,y1,z1,t,R
      logical ifound

c ... local arrays
      real root(2)

c ... local scalars
      real bb,cc,descr,rootmax,rootmin


      bb = 2.*(x1*ux + y1*uy + uz*z1)
c  simplify since aa = 1.
      cc = x1*x1 + y1*y1 + z1*z1 - R*R
      descr = bb*bb - 4.*cc

c      print*,'descr = ',descr,bb,cc
      if (descr.lt.0.) then
         root(1) = -999.
         root(2) = -999.
c         print*,'shit...no intersection'
      else
         descr = sqrt(descr)
         root(1) = (-bb + descr) / 2.
         root(2) = (-bb - descr) / 2.
      endif

      rootmax = max(root(1),root(2))
      rootmin = min(root(1),root(2))

c both roots .lt. 0
      if (rootmax.lt.0.) then
         ifound = .false.
         t = -1.
c one root .lt. 0
      elseif (rootmin.lt.0.) then
         ifound = .true.
         t = rootmax
c both roots ge 0
      else
         ifound = .true.
         t = rootmin
      endif

      return
      end




