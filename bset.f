      subroutine bset(rad,rc,cost,sint,phi,r2p,pi
     1     ,costbz,sintbz,phibz,cospbz,sinpbz,gs,a,b,c,pwr)

c JPS added the possibility of having the B field parallel to the outflow 


c     BAW 2/21/02
c     calculate B-field direction in grid
c     note, rad,rc are normalized to rmax

      implicit none

      real*8 rad,cost,sint,costbz,sintbz,phibz,cospbz,sinpbz,a,b
     1     ,Br,Bth,Btot,Bz,x2,x,phi,r2p,pi,rc,Bphi,z,w,sinpb,phib
     1     ,rcyl,Brt,delphi,pwr,sinp,cosp,Bx,By,Bw,c,gs

      real*8 check

      real*8 a_cav,r_cav,r_cyl_k,z_cav,r_env

      sinp=sin(phi)
      cosp=cos(phi)

      go to 2       !Galli & Shu + toroidal component
c      go to 5       !Galli & Shu + toroidal non-continuous
c      go to 10      !radial field
c      go to 15      !axial field (along z)
c      go to 20      !field along x
c      go to 25      !straight field, arbitrary direction
c      go to 30      !straight field, 30 degree direction

c     Galli & Shu field, zeroth order
      x=gs*rad           !x=.5 at r=rmax
      x2=x**2
      Br=(1.+x2)**2/(4.*x2)*cost
      Bth=-dsqrt((1.+x2)**3/(8.*x))*sint
      Bphi=0.d0
      go to 500

c      Btot=dsqrt(Br**2+Bth**2)
c      Bz=Br*cost-Bth*sint
c      costbz=Bz/Btot
c      sintbz=dsqrt(1.d0-costbz**2)
c      if (cost.gt.0.d0) then
c         phibz=phi
c      else
c         phibz=phi+pi
c         if (phibz.gt.r2p) phibz=phibz-r2p
c      endif
c      cospbz=cos(phibz)
c      sinpbz=sin(phibz)
c      go to 999

 2    continue
c     Galli & Shu + phi component
      x=rad*gs
      x2=x**2
c      a=2.
c      b=0.03
c      c=1.d0
c      pwr=1.5
      rcyl=rad*sint
      z=abs(rad*cost)
c      signz=z/abs(z)
      Br=(1.+x2)**2/(4.*x2)*cost
      Bth=-dsqrt((1.+x2)**3/(8.*x))*sint
      Brt=dsqrt(Br**2+Bth**2)
c      Bphi=c*(a*z/rc)**(-pwr)*(b*rcyl/rc)**(-pwr)*Brt
      Bphi=c*exp(-a*z/rc)*(b*rcyl/rc)**(-pwr)*Brt
      go to 500

c      Bx=Br*sint*cosp+Bth*cost*cosp-Bphi*sinp
c      By=Br*sint*sinp+Bth*cost*sinp+Bphi*cosp
c      Bz=Br*cost-Bth*sint
c      Bw=dsqrt(Bx**2+By**2)
c      Btot=dsqrt(Br**2+Bth**2+Bphi**2)
c      if (Bw.gt.1.d-6) then
c         phibz=datan2(By,Bx)
c      else
c         phibz=phi
c      endif
c      costbz=Bz/Btot
c      sintbz=dsqrt(1.d0-costbz**2)
c      cospbz=cos(phibz)
c      sinpbz=sin(phibz)
      
c      delphi=asin(Bphi/Btot)
c      if (cost.gt.0.d0) then
c         phibz=phi+delphi
c         if (phibz.gt.r2p) phibz=phibz-r2p
c      else
c         phibz=phi-delphi
c         if (phibz.gt.r2p) phibz=phibz-r2p
c      endif

 5    continue
c     Galli & Shu field + toroid, zeroth order
c     make a non-continuous change from G-S to toroidal field at
c     specified radius and z
      b=5.d0*rc
      a=0.4*b
      x=rad/2.           !x=.5 at r=rmaxe
      x2=x**2
      Br=(1.+x2)**2/(4.*x2)*cost
      Bth=-dsqrt((1.+x2)**3/(8.*x))*sint
      z=(rad*cost)
      w=rad*sint
c     Bphi=(z/rc)**(-1)*(2.*rad/rc)**(-1)*10.d0*dsqrt(Br**2+Bth**2)
      if (rad.lt.b.and.abs(z).lt.a) then 
         costbz=0.d0
         sintbz=1.d0
         if (z.gt.0.d0) then
            phibz=phi+pi/2.d0
            if (phibz.gt.r2p) phibz=phibz-r2p
         else
            phibz=phi-pi/2.d0
            if (phibz.lt.0.d0) phibz=phibz+r2p
         endif
      else
         Bphi=0.d0
         Btot=dsqrt(Br**2+Bth**2+Bphi**2)
         Bz=Br*cost-Bth*sint
         costbz=Bz/Btot
         sintbz=dsqrt(1.d0-costbz**2)
c     sinpb=Bphi/Btot
c     phib=asin(sinpb)
         if (cost.gt.0.d0) then
            phibz=phi
         else
            phibz=phi+pi
            if (phibz.gt.r2p) phibz=phibz-r2p
         endif
      endif
      cospbz=cos(phibz)
      sinpbz=sin(phibz)
      go to 999

 10   continue
c     radial B-field

c      costbz=cost
c      sintbz=sint
c      phibz=phi
c      cospbz=cos(phi)
c      sinpbz=sin(phi)
c      go to 999

c     testing, gives same results
      Br=1.
      Bth=0.
      Bphi=0.
      go to 500

 15   continue
c     B-field along z
      costbz=1.d0
      sintbz=0.d0
      phibz=phi
      cospbz=cos(phi)
      sinpbz=sin(phi)
      go to 999

 20   continue
c     B-field along x
c      costbz=0.d0
c      sintbz=1.d0
c      phibz=0.d0
c      cospbz=1.d0
c      sinpbz=0.d0
c     go to 999

c     testing, gives same results, 
      Br=sint*cosp
      Bth=cost*cosp
      Bphi=-sinp

      go to 500

 25   continue
c     arbitrary aligned B-field
      costbz=cos(30.d0*pi/180.d0)
      sintbz=sin(30.d0*pi/180.d0)
      phibz=138.d0*pi/180.d0
      cospbz=cos(phibz)
      sinpbz=sin(phibz)
      go to 999

 30   continue
c     30 degree B-field
      costbz=0.866025404d0
      sintbz=0.5d0
      phibz=0.d0
      cospbz=1.d0
      sinpbz=0.d0
      go to 999

 500   continue

      Bx=Br*sint*cosp+Bth*cost*cosp-Bphi*sinp
      By=Br*sint*sinp+Bth*cost*sinp+Bphi*cosp
      Bz=Br*cost-Bth*sint
      Bw=dsqrt(Bx**2+By**2)

c JPS changed the B field to parallel to the outflow in the outflow:
c Comments from IDL test code:
c Force the cavity to have B in the Z direction only.
c Formula for POLYN models is z=a*r_cyl^1.5 where a=R_env/(R_env*tan(theta))^1.5
c Stark et al. (2006)
c The AFGL 2591 models all had R_env=10^5 AU and theta~10-25.
c theta=13.2, so tan(theta)=0.23455
c theta=11.0 so tan(theta)=0.19438
c theta=7.5 has tan(theta)=0.13165
c We want to look at only a small fraction of the total cavity, so
c since the max of r_cyl=1,

c Calling sequence is quite different!

c      subroutine bset(rad,rc,cost,sint,phi,r2p,pi
c     1     ,costbz,sintbz,phibz,cospbz,sinpbz,gs,a,b,c,pwr)

      r_env=1.
      a_cav=r_env/(r_env*0.13165)**1.5

      r_cav=0.
      z_cav=z
      if (z_cav .gt. 0.) r_cav=(z_cav/a_cav)**0.6666667

      r_cyl_k=rcyl
c  arr=where(r_cyl_k lt r_cavity and r_cyl_k gt 0.,count)
c  print,'k,z_cavity,r_cavity,number of elements in cavity',k,' ',z_cavity,r_cavity,' ',count
c  list=f_where_list(r_cyl_k,arr)

c      for l=0,count-1 do begin
c      bx[list[0,l],list[1,l],k]=0.
c      by[list[0,l],list[1,l],k]=0.
c      endfor

      if (r_cyl_k .lt. r_cav) then
      Bx=0.
      By=0.
      Bw=0.
      endif

c      endfor

c end of JPS section

      Btot=dsqrt(Br**2+Bth**2+Bphi**2)
      if (Bw.gt.1.d-6) then
         phibz=datan2(By,Bx)
         if (phibz.lt.0.d0) phibz=phibz+2.d0*pi
      else
         phibz=phi
      endif
      costbz=Bz/Btot
      sintbz=dsqrt(1.d0-costbz**2)
      cospbz=cos(phibz)
      sinpbz=sin(phibz)

 999  continue




      return
      end


