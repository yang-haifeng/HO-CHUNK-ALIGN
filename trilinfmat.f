      real*8 function trilinfmat(il,im,ix,iy,iz,xval,yval,zval)
c
c     trilinearly interpolate Fmat table.  assumes you already 
c     called a 'locate' subroutine to find position (ix,iy,iz) in array.
c  
c     b.whitney, space science institute
c     dec 12, 2000

c     history
c     00/04/12  (baw) use include file since I couldn't figure out how
c     to read the 4-D array into a 2-D array without creating
c     a big fat 4-D array.
c     00/12/12 (baw) convert bilinear routine into trilinear

      implicit none

      include 'tab.txt'

      integer il,im,ix,iy,iz
      real*8 xx,yy,zz,xval,yval,zval

      xx=(xval-cosbarr(ix))/(cosbarr(ix+1)-cosbarr(ix))
      yy=(yval-cosarr(iy))/(cosarr(iy+1)-cosarr(iy))
      zz=(zval-dphi(iz))/(dphi(iz+1)-dphi(iz))

      trilinfmat=fmat(il,im,ix,iy,iz)*(1.d0-xx)*(1.d0-yy)*(1.d0-zz)
     1     +fmat(il,im,ix+1,iy,iz)*xx*(1.d0-yy)*(1.d0-zz)
     1     +fmat(il,im,ix,iy+1,iz)*(1.d0-xx)*yy*(1.d0-zz)
     1     +fmat(il,im,ix,iy,iz+1)*(1.d0-xx)*(1.d0-yy)*zz
     1     +fmat(il,im,ix+1,iy+1,iz)*xx*yy*(1.d0-zz)
     1     +fmat(il,im,ix+1,iy,iz+1)*xx*(1.d0-yy)*zz
     1     +fmat(il,im,ix,iy+1,iz+1)*(1.d0-xx)*yy*zz
     1     +fmat(il,im,ix+1,iy+1,iz+1)*xx*yy*zz

      return
      end

