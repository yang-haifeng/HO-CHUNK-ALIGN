c     ******************************************************

      subroutine testgrid(xp,yp,zp,rsq,rtot,ii,routine,iphot,r2p)

c     test to see if photon position corresponds to grid
c     position, ii(1),ii(2),ii(3)

c     input:
c        xp,yp,zp,rsq  photon position
c        ii(3), grid indices
c        routine, name of subroutine called from
c        iphot, photon number

c     output:
c       ii(3), updated if necessary.

c     ******************************************************

      implicit none

      include 'grid.txt'

      real*8 xp,yp,zp,rsq,beta,lp1,r2p,cosb,rtot
      integer ii(3),ir,it,ip,iphot
      character routine*20

      ir=ii(1)
      it=ii(2)
      ip=ii(3)
c     print*,'rsq,iphot,r2arr(ir),r2arr(ir+1)'
c     print*,rsq,iphot
c     print*,r2arr(ir)
c     print*,r2arr(ir+1)
c     check
c     rsq should be between r2arr(ir) and r2arr(ir+1)
      if (rsq.lt.r2arr(ir)) then
         if (rsq.lt.r2arr(ir-1)) then
            print*, 'ERROR, rsq way off! ',routine
            print*,'rtot,rarr(ir),ir',rtot,rarr(ir),ir
            call locate(r2arr,nrg,nrg,rsq,ir)
            if (ir.lt.1) then
               print*,'ERROR testgrid, ir out of limits'
               print*,'rsq,ir',rsq,ir
               ir=1
            endif
            if (ir.gt.nrg) then
               print*,'ERROR testgrid, ir out of limits'
               print*,'rsq,ir',rsq,ir
               ir=nrg
            endif
            print*,'new ir',ir
            ii(1)=ir
         else
            print*,'WARNING, setting ir=ir-1 ',routine
            print*,'iphot,rsq, rsq(ir),rsq(ir-1),ir ',
     $           iphot,rsq,r2arr(ir),r2arr(ir-1),ir
            ir=ir-1
            ii(1)=ir
         endif
      else if (rsq.gt.r2arr(ir+1)) then
         if (rsq.gt.r2arr(ir+2)) then
            print*,'ERROR, rsq way off! too high',routine
            print*,'iphot,rtot,rarr(ir),ir',
     $           iphot,rtot,rarr(ir),ir
            call locate(r2arr,nrg,nrg,rsq,ir)
            if (ir.lt.1) then
               print*,'ERROR testgrid, ir out of limits'
               print*,'rsq,ir',rsq,ir
               ir=1
            endif
            if (ir.gt.nrg) then
               print*,'ERROR testgrid, ir out of limits'
               print*,'rsq,ir',rsq,ir
               ir=nrg
            endif
            print*,'new ir',ir
            ii(1)=ir
         else
            print*,'WARNING, setting ir=ir+1',routine
            print*,'iphot,rsq, rsq(ir),rsq(ir+1),ir ',
     $           iphot,rsq,r2arr(ir),r2arr(ir+1),ir
            ir=ir+1
            ii(1)=ir
         endif
      endif

      if (ntg.gt.1) then
c     beta should be between thetarr(it) and thetarr(it+1)
c     cos(beta) should be between costarr(it) and costarr(it+1)
c     beta=dacos(zp/dsqrt(rsq))
         cosb=zp/rtot
         if (cosb.gt.costarr(it)) then
            if (cosb.gt.costarr(it-1)) then
               beta=dacos(cosb)
               print*, 'ERROR, beta way off! ',routine
               print*,'iphot',iphot
               print*,'zp,rtot,cosb,beta',zp,rtot,cosb,beta
               print*,'it,beta,thetarr(it),thetarr(it-1)',
     1              it,beta,thetarr(it),thetarr(it-1) 
               call locate(thetarr,ntg,ntg,beta,it)
               if (it.lt.1) then
                  print*,'ERROR testgrid, it out of limits'
                  print*,'beta,it',beta,it
                  it=1
               endif
               if (it.gt.ntg) then
                  print*,'ERROR testgrid, it out of limits'
                  print*,'beta,it',beta,it
                  it=ntg
               endif
               print*,'new it,thet(it) ',it,thetarr(it)
               ii(2)=it
            else
               print*,'WARNING, setting it=it-1 ',routine
               print*,'iphot,cosb,costarr(it),costarr(it-1),it'
               print*,iphot,cosb,costarr(it),costarr(it-1),it
               print*,'r,x,y,z',rtot,xp,yp,zp
               it=it-1
               ii(2)=it
            endif
         else if (cosb.lt.costarr(it+1)) then
            if (cosb.lt.costarr(it+2)) then
               beta=dacos(cosb)
               print*,'ERROR, beta way off! ',routine
               print*,'iphot',iphot
               print*,'zp,rtot,cosb,beta',zp,rtot,cosb,beta
               print*,'it,beta,thetarr(it),thetarr(it+1)',
     1              it,beta,thetarr(it),thetarr(it+1) 
               call locate(thetarr,ntg,ntg,beta,it)
               if (it.lt.1) then
                  print*,'ERROR testgrid, it out of limits'
                  print*,'beta,it',beta,it
                  it=1
               endif
               if (it.gt.ntg) then
                  print*,'ERROR testgrid, it out of limits'
                  print*,'beta,it',beta,it
                  it=ntg
               endif
               print*,'new it,thet(it) ',it,thetarr(it)
               ii(2)=it
            else
               print*,'WARNING, setting it=it+1 ',routine
               print*,'iphot,cosb,costarr(it),costarr(it+1),it'
               print*,iphot,cosb,costarr(it),costarr(it+1),it
               print*,'r,x,y,z',rtot,xp,yp,zp
               it=it+1
               ii(2)=it
            endif
         endif
      endif

c     do the same for phi
c     l should be between phiarr(ip) and phiarr(ip+1)
      if (npg.gt.1) then
         lp1=datan2(yp,xp)
         if (lp1.le.0.d0) lp1=lp1+r2p
         if (lp1.lt.phiarr(ip)) then
            if (lp1.lt.phiarr(ip-1)) then
               print*, 'tauint, ERROR, longitude way off! '
               call locate(phiarr,npg,npg,lp1,ip)
               if (ip.lt.1) then
                  print*,'ERROR testgrid, ip  out of limits'
                  print*,'lp1,ip',lp1,ip
                  ip=1
               endif
               if (ip.gt.npg) then
                  print*,'ERROR testgrid, ip  out of limits'
                  print*,'lp1,ip',lp1,ip
                  ip=npg
               endif
               print*,'new ip',ip
               ii(3)=ip
            else
               ip=ip-1
               ii(3)=ip
               print*,'WARNING, setting ip=ip-1 '
               print*,'long, phiarr(ip-1),ip-1 ',lp1,
     $              phiarr(ip),ip
            endif
         else if (lp1.gt.phiarr(ip+1)) then
            if (lp1.gt.phiarr(ip+2)) then
               print*,'tauint, ERROR, longitude way off! '
               call locate(phiarr,npg,npg,lp1,ip)
               if (ip.lt.1) then
                  print*,'ERROR testgrid, ip  out of limits'
                  print*,'lp1,ip',lp1,ip
                  ip=1
               endif
               if (ip.gt.npg) then
                  print*,'ERROR testgrid, ip  out of limits'
                  print*,'lp1,ip',lp1,ip
                  ip=npg
               endif
               print*,'new ip',ip
               ii(3)=ip
            else
               ip=ip+1
               ii(3)=ip
               print*,'WARNING, setting ip=ip+1'
               print*,'long, phiarr(ip+1),ip+1 ',
     $              lp1,phiarr(ip),ip
            endif
         endif
      endif

      return
      end


