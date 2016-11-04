      subroutine findangle5(pi,r2p,cost1,sint1,phi1,cosp1,sinp1
     1 ,cost2,sint2,phi2,cosp2,sinp2,costnew,sintnew,phinew,ri1)

c     HISTORY
c     2/16/02 BAW, modified findangle3 for B-field case. 
c     define stokes rotation w.r.t. theta2, not theta1

c     9/7/99  BAW, add i1 angle for stokes scattering diagram.

c	input:
c	cost1,sint1,phi1    -- B-field direction (t1,phi1)
c	cost2,sint2,phi2    -- photon direction (t2,phi2)
c						  
c	output:
c	costnew,sintnew,phinew,ri1   --- photon direction w.r.t. B-field

c	see MARS notebook p. 4-7 for drawings of coordinate systems


      implicit none
      real*8 pi,r2p,cost1,sint1,cost2,sint2,phi1,phi2,costnew
     1     ,sintnew,phinew,diff,cosx,cosp1,cosp2,sinp1,sinp2,ri1,x
     1     ,y,cosy
c      integer flag
      character*20 routine
      real*8 check

      routine='findangle5'

      costnew=cost1*cost2+(sint1*sint2*(cosp1*cosp2+sinp1*sinp2))
      costnew=check(costnew,routine)
      sintnew=dsqrt(1.d0-costnew**2)

c      sint1 and sintnew should both be greater than 0.  but they could be 0.
      if (sint1.gt.1.d-8) then
        if (sintnew.gt.1.d-8) then
          cosy=(cost2-(cost1*costnew))/(sint1*sintnew)
          cosy=check(cosy,routine)
          y=dacos(cosy)
          if (phi2.gt.phi1) then
            diff=phi2-phi1
            if (diff.gt.pi) then
              diff=r2p-diff
              phinew=pi+y
            else
              phinew=pi-y
            endif
          else
            diff=phi1-phi2
            if (diff.gt.pi) then
              diff=r2p-diff
              phinew=pi-y
            else
              phinew=pi+y
            endif
          endif
        else
           phinew=0.d0       !should be random, doesn't matter, means at pole 
        endif
      else
        phinew=phi2
      endif

c     position angle
      if (sint2.gt.1.e-8) then
         if (sintnew.gt.1.e-8) then
            cosx=(cost1-cost2*costnew)/sint2/sintnew
            cosx=check(cosx,routine)
            x=dacos(cosx)
            if (phi2.gt.phi1) then
               diff=phi2-phi1
               if (diff.lt.pi) then
                  ri1=r2p-x
               else
                  ri1=x
               endif
            else
               diff=phi1-phi2
               if (diff.lt.pi) then
                  ri1=x
               else
                  ri1=r2p-x
               endif
            endif
         else
c     sintnew=sin(thetab)=0, you are lined up with field.
c     not sure what this means...x has no meaning
            ri1=0.d0
         endif
      else
c     sint2=0 means theta=0; I think pa changes by 180 but that's 0
         ri1=0.d0
      endif
               
      return
      end
