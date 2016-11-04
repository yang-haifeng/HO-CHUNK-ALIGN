      SUBROUTINE locate(xx,nn,n,x,j)
c     from numerical recipes
c     searches an ordered table, using bisection

c     baw 8/1/2000, allow dimensioned array size to be bigger than size used.
c     dimensioned array size is nn; usable is n.

      INTEGER j,n,nn
      REAL*8 x,xx(nn)
      INTEGER jl,jm,ju
      jl=0
      ju=n+1
10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      goto 10
      endif
      j=jl
      return
      END
