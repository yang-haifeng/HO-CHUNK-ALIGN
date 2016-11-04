      subroutine table
      
c	makes a table of integrated probability distribution function.
c	for limb darkening
	
      implicit none

      include 'tabl.txt'
      integer i
	
      delmu=1.0/(float(ntab))
      xmu(1)=0.0
	
      do i=2,ntab
        xmu(i)=delmu*float(i)
        prob(i)=0.5*(xmu(i)**3+xmu(i)**2)
      end do
	
      write(6,*) 'in table, prob(n), should equal 1, ',prob(ntab)
      print*,'ntab',ntab
	
      return
      end
	
	
	
c	*********************************************************

c	*********************************************************

	subroutine darkening(random,mu)

c	interpolates from table (subroutine table) to sample
c	mu from probability distribution function.
c	mu ranges from 0 to 1 

      implicit none

	include 'tabl.txt'
	
	integer i,n
	real*8 random,mu,check
	
c	do i=1,ntab
c	  check=random-prob(i)
c	  if (check.lt.0.0) then
c	    n=i
c	    go to 10
c	  end if
c	end do

c     locate index in prob array
	call locate(prob,ntab,ntab,random,n)
	
c	write(6,*) 'table for mu must be wrong...'
	
c10	mu=xmu(n)-(prob(n)-random)/(prob(n)-prob(n-1))*delmu
10	mu=xmu(n+1)-(prob(n+1)-random)/(prob(n+1)-prob(n))*delmu

      if (mu.gt.1.0.or.mu.lt.0.) then
	print*,'dammit'
	continue
	endif

	return
	end

