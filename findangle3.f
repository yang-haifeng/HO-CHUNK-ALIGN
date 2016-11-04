      subroutine findangle3(pi,r2p,cost1,sint1,phi1,cosp1,sinp1
     1 ,cost2,sint2,phi2,cosp2,sinp2,costnew,sintnew,phinew,ri1)

c     HISTORY
c     9/7/99  BAW, add i1 angle for stokes scattering diagram.

c	input:
c	cost1,sint1,phi1    --incident photon direction (t1,phi1)
c						  (or normal to surface)
c	cost2,sint2,phi2    --scattered photon direction (t2,phi2)
c						  
c	output:
c	costnew,sintnew,phinew   ---tnew is scattering angle in scattering frame
c					 ---phinew is azimuthal angle in frame of scattering
c	see MARS notebook p. 4-7 for drawings of coordinate systems

c     note:  when calling for pathfinder and MGS coords, the positions are
c     theta1 and phi1, and the photon angles are theta2 and phi2.

      implicit none
      real*8 pi,r2p,cost1,sint1,cost2,sint2,phi1,phi2,costnew,
     1  sintnew,phinew,diff,cosx,cosp1,cosp2,sinp1,sinp2,ri1,x
c      integer flag
      character*20 routine
      real*8 check

      routine='findangle'

      costnew=cost1*cost2+(sint1*sint2*(cosp1*cosp2+sinp1*sinp2))
      costnew=check(costnew,routine)
      sintnew=dsqrt(1.d0-costnew**2)

c      sint1 and sintnew should both be greater than 0.  but they could be 0.
      if (sint1.gt.1.d-8) then
        if (sintnew.gt.1.d-8) then
          cosx=(cost2-(cost1*costnew))/(sint1*sintnew)
          cosx=check(cosx,routine)
          x=dacos(cosx)
          if (phi2.gt.phi1) then
            diff=phi2-phi1
            if (diff.gt.pi) then
              diff=r2p-diff
              ri1=x
              phinew=pi+x
            else
              phinew=pi-x
              ri1=r2p-x
            endif
          else
            diff=phi1-phi2
            if (diff.gt.pi) then
              diff=r2p-diff
              phinew=pi-x
              ri1=r2p-x
            else
              phinew=pi+x
              ri1=x
            endif
          endif
        else
          phinew=0.d0        !should be random, doesn't matter, means at pole 
        endif
      else
        phinew=phi2
      endif

      return
      end