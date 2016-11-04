      double precision function dcosd(angle)

      implicit none

c ... arguments
      real*8 angle

c ... local
      real*8 ang,pi

      pi = 3.1415926535898d0

      ang = angle * pi / 180.d0
      dcosd = dcos(ang)

      return
      end

      double precision function dsind(angle)

      implicit none

c ... arguments
      real*8 angle
c ... local
      real*8 ang,pi

      pi = 3.1415926535898d0

      ang = angle * pi / 180.d0
      dsind = dsin(ang)

      return
      end

      double precision function dtand(angle)

      implicit none

c ... arguments
      real*8 angle
c ... local
      real*8 ang,pi

      pi = 3.1415926535898d0

      ang = angle * pi / 180.d0
      dtand = dtan(ang)

      return
      end
