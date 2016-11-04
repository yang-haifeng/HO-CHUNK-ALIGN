      subroutine find_wall_2D(dt,ioffset,r2i,r2ip1,
     &     tan_th2j,tan_th2jp1,ux,uy,uz,x1,y1,z1)
c
c   This routine calculates the distance to all faces of a 2D grid cell
c   (in a spherical -- r,theta) from a point (x1,y1,z1) along the
c   directional cosine vector (ux,uy,yz).
c     NOTE: ***** it is necessary to flag tan_th2 array using a -1.
c               to indicate theta=90.
c        Program ASSUMES that theta=90 can occur in ONLY ONE of the
c        tan_th2 values.
c
c  Output values:
c     dt = distance from (x1,y1,z1) along (ux,uy,uz) to closest face(s),
c          along for possibility of passing through cell vertex (i.e.,
c     ioffset = integer offsets for new cell position.  
c          intersecting two constant surfaces).
c
c  Input values:
c     r2i,r2ip1 = squared radii of two constant surfaces r(i),r(i+1)
c     tan_th2j,jp1 = (tan(theta))**2 for two constant surfaces
c        theta(j),theta(j+1).  see NOTE above.
c     ux,uy,uz = directional cosines along x,y,z directions
c     x1,y1,z1 = point from which distance is calculated
c
c m.j. wolff/b.a. whitney, 2000/03/16 (created from 3D version)
c history:
c


      implicit none
c ... 
      integer ioffset(3)
      real*8 dt,r2i,r2ip1,tan_th2j,tan_th2jp1,ux,uy,uz,x1,y1,z1

c ..
      integer i,iface
      real*8 t
      real*8 dist(3)
      integer idist(3)
      logical ifound

c find shortest non-negative distance to two constant radius surfaces
c where negative ifound = .false. indicates no intersection    
      call radface(iface,ifound,t,ux,uy,uz,x1,y1,z1,r2i,r2ip1)
      dist(1) = t
      idist(1) = iface
c     print*,'RAD - iface,ifound,t:',iface,ifound,t

c find shortest positive distance to two constant theta surfaces
c where negative ifound = .false. indicates no intersection    
      call thetaface(iface,ifound,t,ux,uy,uz,x1,y1,z1,tan_th2j,
     &     tan_th2jp1)
c     print*,'THETA - iface,ifound,t:',iface,ifound,t
      dist(2) = t
      idist(2) = iface

c find minimum value of t (distance) which is still ge 0
c (of course, assume that at least one non-zero value in DIST)
c t=0 is NOT allowed to be "shortest distance".  offsets
c are either -1 or 1.  if "ri" is closest face, then idist=0.
c if "rip1" is closest, then idist=1
      t = max(dist(1),dist(2))
      do i=1,2
         if ( (dist(i).gt.0.d0).AND.(dist(i).lt.t) ) then
              t = dist(i)
          endif
      enddo
      do i=1,2
         ioffset(i) = 0
         if (dist(i).eq.t) then
c "shortest distance" offsets are either -1 or 1
            ioffset(i) = -1 + 2*idist(i)
         endif
      enddo
c set "winning distance" to dt to pass back
      dt = t

      ioffset(3)=0

      return
      end
