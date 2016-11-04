c     ***********************************************************
      subroutine densenv(rad,thetin,costin,sina,pi,pihalf,phi
     1     ,dens)
      
c     calculates density in envelope. called during grid setup
c     so doesn't need to be optimized.

c  history:
c  00/09/06 (baw): write subroutine, modified from opacenv.f
c  01/02/17 (baw): combine opacitybub subroutine into here.

      implicit none
      include 'tts.txt'
      include 'opacin.txt' 

      real*8 rad,dens,cosa,rado,radi,phi,thet,thetin,pihalf,y
     $     ,pi,r2,costin,m,n,fact,sina,xmu,rx,rp,zup,xmu0,factor
     $     ,rx2,xmu0new,zlo,zp

      integer iflag

c ... begin g77
      real*8 dcosd,dsind,dtand
      external dcosd,dsind,dtand
c     ... end g77

      dens=2.0833e-17

      return 
      end

