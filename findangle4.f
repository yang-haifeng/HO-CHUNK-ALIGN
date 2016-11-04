      subroutine findangle4(pi,r2p,cost1,sint1,phi1,cosp1,sinp1
     1 ,cost2,sint2,phi2,cosp2,sinp2,costnew,sintnew,phinew,ri1)

c     HISTORY
c     2/14/02 BAW, modified findangle3 to rotate stokes back
c     from magnetic field direction to observer frame.

c     9/7/99  BAW, add i1 angle for stokes scattering diagram.

c	input:
c	cost1,sint1,phi1    --B-field direction
c	cost2,sint2,phi2    --photon dir. w.r.t. B-field
c						  
c	output:
c	costnew,sintnew,phinew,ri1   ---tnew is scattering angle in scattering frame
c					 ---phinew is azimuthal angle in frame of scattering
c	see MARS notebook p. 4-7 for drawings of coordinate systems
c     new drawings on pieces of paper...(2/14/02)

c     note:  when calling for pathfinder and MGS coords, the positions are
c     theta1 and phi1, and the photon angles are theta2 and phi2.

      implicit none

      real*8 pi,r2p,cost1,sint1,cost2,sint2,phi1,phi2,costnew,
     1     sintnew,phinew,diff,cosx,cosp1,cosp2,sinp1,sinp2,ri1,x
     1     ,cosdel,del
c     integer flag
      character*20 routine
      real*8 check
      
      routine='findangle4'

c     from drawings, let thetabz=theta1, thetab=theta2, etc.

c     in all cases, cos(y)=-cospb=-cosp2 (see drawings)
      costnew=cost1*cost2-sint1*sint2*cosp2
      costnew=check(costnew,routine)
      sintnew=dsqrt(1.d0-costnew**2)

      if (sint1.gt.1.e-8) then
         if (sintnew.gt.1.e-8) then
            cosdel=(cost2-(cost1*costnew))/(sint1*sintnew)
            cosdel=check(cosdel,routine)
            del=dacos(cosdel)
            if (phi2.lt.pi) then
               phinew=del+phi1
               if (phinew.gt.r2p) phinew=phinew-r2p
            else
               phinew=phi1-del
               if (phinew.lt.0.) phinew=phinew+r2p
            endif
         else
            phinew=0.d0     !at the pole
         endif
      else
         phinew=phi2
      endif
         
      if (sint2.gt.1.e-8) then
         if (sintnew.gt.1.e-8) then
            cosx=(cost1-cost2*costnew)/sint2/sintnew
            cosx=check(cosx,routine)
            x=dacos(cosx)
            if (phinew.gt.phi1) then
               diff=phinew-phi1
               if (diff.lt.pi) then
                  ri1=x
               else
                  ri1=r2p-x
               endif
            else
               diff=phi1-phinew
               if (diff.lt.pi) then
                  ri1=r2p-x
               else
                  ri1=x
               endif
            endif
         else
            ri1=0.d0            !lined up with B-z axis I think
         endif
      else
         ri1=0.d0            !B along z, no stokes rotation 
      endif

      return
      end



