      subroutine find_wall(dt,ioffset,A1,A2,B1,B2,C1,C2,D1,D2,r2i,r2ip1,
     &     tan_th2j,tan_th2jp1,ux,uy,uz,x1,y1,z1)
c
c   This routine calculates the distance to all faces of a grid cell
c   (in a spherical -- r,theta,phi) from a point (x1,y1,z1) along the
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
c     A1...D1,A2...D2 = plane equation coefficients for two constant surfaces
c     r2i,r2ip1 = squared radii of two constant surfaces r(i),r(i+1)
c     tan_th2j,jp1 = (tan(theta))**2 for two constant surfaces
c        theta(j),theta(j+1).  see NOTE above.
c     ux,uy,uz = directional cosines along x,y,z directions
c     x1,y1,z1 = point from which distance is calculated
c
c m.j. wolff/b.a. whitney, 2000/02/04
c history:
c 00/02/09 (mjw):  add correct algorithm for finding minimum non-negative
c                      distance.
c                  add documentation (such as it is)
c 00/02/18 (mjw):  trap case of point on wall surface (i.e. t=0 ) and 
c                     force ioffset=0.
c 00/02/28 (mjw):  removed t=0 as possible choice (RADFACE will not pass
c                     back as such).  NOW: only offset values are
c                     -1 or 1.
c 00/03/16 (mjw):  pass r**2 instead r
c                  double precision
c


      implicit none
c ... 
      integer ioffset(3)
      real*8 A1,A2,B1,B2,C1,C2,D1,D2,dt
      real*8 r2i,r2ip1,tan_th2j,tan_th2jp1,ux,uy,uz,x1,y1,z1

c ..
      integer i,iface
      real*8 t
      real*8 dist(3)
      integer idist(3)
      logical ifound

c     print*, 'in find_wall '
c     print*, 'inputs'
c     print*,'a1,a2,b1,b2 ',a1,a2,b1,b2
c     print*,'c1,c2,d1,d2 ',c1,c2,d1,d2
c     print*,'r2i,r2ip1 ',r2i,r2ip1 
c     print*,'tan_th2j,tan_th2jp1 ', tan_th2j,tan_th2jp1
c     print*,'ux,uy,uz ',ux,uy,uz
c     print*,'y1,z1 ',y1,z1
c find shortest non-negative distance to two constant radius surfaces
c where negative ifound = .false. indicates no intersection    
      call radface(iface,ifound,t,ux,uy,uz,x1,y1,z1,r2i,r2ip1)
      dist(1) = t
      idist(1) = iface
c     if (t.lt.0.d0) then
c     print*,'RAD - iface,ifound,t:',iface,ifound,t
c     endif

c find shortest positive distance to two constant theta surfaces
c where negative ifound = .false. indicates no intersection    
      call thetaface(iface,ifound,t,ux,uy,uz,x1,y1,z1,tan_th2j,
     &     tan_th2jp1)
c     if (t.lt.0.d0) then
c     print*,'THETA - iface,ifound,t:',iface,ifound,t
c     endif
      dist(2) = t
      idist(2) = iface

c find shortest non-negative distance to two constant phi surfaces
c where negative ifound = .false. indicates no intersection    
      call phiface(iface,ifound,t,A1,B1,C1,D1,A2,B2,
     &     C2,D2,ux,uy,uz,x1,y1,z1)
c     if (t.lt.0.d0) then
c     print*,'PHI - iface,ifound,t:',iface,ifound,t
c     endif
      dist(3) = t
      idist(3) = iface

c find minimum value of t (distance) which is still ge 0
c (of course, assume that at least one non-zero value in DIST)
c t=0 is NOT allowed to be "shortest distance".  offsets
c are either -1 or 1.  if "ri" is closest face, then idist=0.
c if "rip1" is closest, then idist=1
      t = max(dist(1),dist(2),dist(3))
      do i=1,3
         if ( (dist(i).gt.0.d0).AND.(dist(i).lt.t) ) then
              t = dist(i)
          endif
      enddo
      do i=1,3
         ioffset(i) = 0
         if (dist(i).eq.t) then
c "shortest distance" offsets are either -1 or 1
            ioffset(i) = -1 + 2*idist(i)
         endif
      enddo
c set "winning distance" to dt to pass back
      dt = t

      return
      end


