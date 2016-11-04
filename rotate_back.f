      subroutine rotate_back(ir,it,ip,pi,r2p,sqp,sup,sqpb,supb
     1     ,cost,sint,phi,cosp,sinp
     1     ,costb,sintb,phib,cospb,sinpb)

c     BAW, 02/14/02 rotate Stokes vectors back from B-field frame
c     to observer frame.

c     input: costb,sintb,phib,cospb,sinpb,sqpb,supb
c     output: cost,sint,phi,cosp,sinp,sqp,sup

c     calculated b-field angles in subroutine rotate.  should not
c     have changed.  could redo here, since it's just an array call.

      implicit none

      real*8 pi,r2p,cost,sint,phi,cosp,sinp
     1     ,costbz,sintbz,phibz,cospbz,sinpbz,sinpb,cospb
     1     ,costb,sintb,phib,cos2,sin2,ir1,sip,sup,sqp,sqpb,supb
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

c      print*,costbz,sintbz,phibz,cospbz,sinpbz
c      print*,''

      call findangle4(pi,r2p,costbz,sintbz,phibz,cospbz,sinpbz
     1     ,costb,sintb,phib,cospb,sinpb
     1     ,cost,sint,phi,ir1)

      cosp=dcos(phi)
      sinp=dsin(phi)
      
      if (ir1.gt.0) then
         cos2=dcos(2.d0*ir1)
         sin2=dsin(2.d0*ir1)
         sqp=cos2*sqpb+sin2*supb
         sup=-sin2*sqpb+cos2*supb
      else
         sqp=sqpb
         sup=supb
      endif

      return
      end
