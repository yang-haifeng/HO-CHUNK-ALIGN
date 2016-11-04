      function check(x, routine)

c     x is cos(theta) or sin(theta)
c     check to see if it's greater than 1

      implicit none
c      real*8 x
      real*8 x
      real*8 check
      character*20 routine,cmsger*50
c      intrinsic sngl
      
      if (x .gt. 1.d0) then
        if (x .gt. 1.01d0) then
            write(cmsger,'(a,f14.11)')'cos,sin > 1',x
            call ERRMSG('WARNING',routine,cmsger)
            print*,'hi'
        end if
        x = 1.d0
      else if (x .lt. -1.d0) then
        if (x .lt. -1.01d0) then
            print*,'hi'
            write(cmsger,'(a,f14.11)')'cos,sin < -1',x
            call ERRMSG('WARNING',routine,cmsger)
        end if
        x = -1.d0
      end if

c      check = sngl(x)
      check=x

      return
      end



