      subroutine rotate(ir,it,ip,pi,r2p,sqp,sup,sqpb,supb
     1     ,cost,sint,phi,cosp,sinp
     1     ,costb,sintb,phib,cospb,sinpb)

c     BAW, 02/14/02 rotate Stokes vectors into frame with B-field
c     along z.

c     input:  ir,it,ip,pi,r2p,sqp,sup,cost,sint,phi,cosp,sinp
c     output:  sqpb,supb,costb,sintb,phib,cospb,sinpb

      implicit none

      real*8 pi,r2p,cost,sint,phi,cosp,sinp
     1     ,costbz,sintbz,phibz,cospbz,sinpbz,sinpb,cospb
     1     ,costb,sintb,phib,cos2,sin2,ir1,sup,sqp,sqpb,supb
      integer ir,it,ip

      include 'grid.txt'

c     find angle of B-field w.r.t. z-axis

      if (ir.ge.nrg.or.ir.lt.1) print*,'ERROR ir off in rotate'
      if (it.lt.0.or.it.gt.ntg) print*,'ERROR it off in rotate'
      if (ip.lt.0.or.ip.gt.npg) print*,'ERROR ip off in rotate'

      costbz=costbzarr(ir,it,ip)
      sintbz=sintbzarr(ir,it,ip)
      phibz=phibzarr(ir,it,ip)
      cospbz=cospbzarr(ir,it,ip)
      sinpbz=sinpbzarr(ir,it,ip)

c     print*,costbz,sintbz,phibz,cospbz,sinpbz
c      print*,ir,it,ip,costbz
c      print*,''

c     now find angle between photon and field

      call findangle5(pi,r2p,costbz,sintbz,phibz,cospbz,sinpbz
     1     ,cost,sint,phi,cosp,sinp
     1     ,costb,sintb,phib,ir1)
      
      cospb=dcos(phib)
      sinpb=dsin(phib)

      if (ir1.gt.0.d0) then
         cos2=dcos(2.d0*ir1)
         sin2=dsin(2.d0*ir1)
         sqpb=cos2*sqp+sin2*sup
         supb=-sin2*sqp+cos2*sup
      else
         sqpb=sqp
         supb=sup
      endif

      return
      end


