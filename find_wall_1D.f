      subroutine find_wall_1D(dt,ioffset,r2i,r2ip1,
     &     ux,uy,uz,x1,y1,z1)
c
c   This routine calculates the distance to all faces of a 1D grid cell
c   (in a spherical -- r) from a point (x1,y1,z1) along the
c   directional cosine vector (ux,uy,yz).
c
c  Output values:
c     dt = distance from (x1,y1,z1) along (ux,uy,uz) to closest face(s),
c          along for possibility of passing through cell vertex (i.e.,
c     ioffset = integer offsets for new cell position.  
c          intersecting two constant surfaces).
c
c  Input values:
c     r2i,r2ip1 = squared radii of two constant surfaces r(i),r(i+1)
c     ux,uy,uz = directional cosines along x,y,z directions
c     x1,y1,z1 = point from which distance is calculated
c
c m.j. wolff/b.a. whitney, 2000/03/16 (created from 3D version)
c history:
c     2000/12/27 (baw) created from 2D version
c


      implicit none
c ... 
      integer ioffset(3)
      real*8 dt,r2i,r2ip1,ux,uy,uz,x1,y1,z1

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
      ioffset(1) = -1 + 2*idist(1)
      dt = t

c      do i=1,1
c         ioffset(i) = 0
c         if (dist(i).eq.t) then
c "shortest distance" offsets are either -1 or 1
c            ioffset(i) = -1 + 2*idist(i)
c         endif
c      enddo
c set "winning distance" to dt to pass back
c      dt = t

      ioffset(3)=0
      ioffset(2)=0

      return
      end
