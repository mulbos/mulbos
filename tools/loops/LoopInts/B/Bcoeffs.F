* Bcoeffs.F

#include "defs.h"

*#define DEBUG_Bcoeffs

***********************************************************************
	subroutine CYPBcoeff(para, B, ldpara)
***********************************************************************
* Possible approximation here to avoid numerical cancellation for 
* m1 = m2 << p. See CYPSetM0(). This is a very good approximation.
* Author: Le Duc Ninh (2011)
***********************************************************************
	implicit none
	integer ldpara
	CVAR para(ldpara,Pbb)
	CVAR paraC(Pbb)
	CVAR B(Nbb)
	integer i
	logical dump

#include "lt.h"

***[DEBUG
	dump = ibits(debugkey, DebugB, 1) .ne. 0
	if( dump ) call CYPDumpPara(2, para, ldpara, "CYPBcoeff")
***]DEBUG

	do i = 1, Nbb
	 B(i) = 0
	enddo
	
	do i = 1, Pbb
	 paraC(i) = para(1, i)
	enddo

***[set small masses to 0 when it's safe
	call CYPSetM0(paraC, Pbb)
***]set small masses to 0 when it's safe
	
	call CYPBscalar(paraC, B(bb0))	

	call CYPBtensor(paraC, B)

***[DEBUG
	if( dump ) call CYPDumpCoeff(2, B)
***]DEBUG
	end
***********************************************************************

	subroutine CYPBscalar(para, fx)
	implicit none
	CVAR para(Pbb), paraC(Pbb), fx

#include "lt.h"

#ifdef DEBUG2
	print *, "CYPBscalar in Bcoeffs.F:"
	print *, "    p1   =", para(1)
	print *, "    m1 =", para(2)
	print *, "    m2 =", para(3)
#endif

	iers = 2

	if ( abs(para(3)) .le. 0 ) then
	paraC(1) = para(1)
	paraC(2) = para(3)
	paraC(3) = para(2)
	else
	paraC(1) = para(1)
	paraC(2) = para(2)
	paraC(3) = para(3)
	endif

        call CYPBnx(0, paraC, fx)

***[check for errors
	if ( mod(iers,4) .eq. 0 ) then 
	print *, "CYPBscalar in Bcoeffs.F:"
	print *, "    p1   =", para(1)
	print *, "    m1 =", para(2)
	print *, "    m2 =", para(3)
	print *, "    B0 =", fx
	endif
***]check for errors
	iers = 1

	end

***************************************************************
	subroutine CYPBtensor(para, B)
**********************************************************
* Method: If abs(p) is large enough we use Passarino-Veltman 
* method, otherwise the coefficients B(bb11..1) are calculated 
* directly like the B0 function.   
**********************************************************
	implicit none
	integer i
	CVAR para(Pbb)
	CVAR p, m1, m2, f
	CVAR B(Nbb)

	CVAR CYPA0s, CYPA00s
	external CYPA0s, CYPA00s
	
#include "lt.h"

	do i = 2, Nbb
	 B(i) = 0
	enddo

	if ( lambda .eq. -2 ) return

	p = para(1)
	m1 = para(2)
	m2 = para(3)

	f = p - m2 + m1
	mmax = max(abs(m1), abs(m2))

*** derivatives
* can go up to rank 5 or more
        call CYPdBnx(0, para, B(dbb0))
        call CYPdBnx(1, para, B(dbb1))
        call CYPdBnx(2, para, B(dbb11))

***[Rank 1
	if ( abs(p) .le. stab2*mmax ) then
        call CYPBnx(1, para, B(bb1))
	else 
	B(bb1) = (CYPA0s(m1) - CYPA0s(m2) - f*B(bb0))/p/2
	endif

*** derivatives (p cannot be zero here)
        if ( abs(B(dbb1) - nan) .le. 0 )
     &	B(dbb1) = (-f*B(dbb0) - B(bb0))/p/2 - B(bb1)/p
***]Rank 1

***[Rank 2
	B(bb00) = (CYPA0s(m2) + 2*m1*B(bb0) + f*B(bb1) 
     &	 + (m1 + m2 - p/3)*xrat)/6

*** derivatives
	B(dbb00) = (2*m1*B(dbb0) + f*B(dbb1) + B(bb1)  
     &	 - 1/3D0*xrat)/6

	if ( abs(p) .le. stab2*mmax ) then
        call CYPBnx(2, para, B(bb11))
	else
	B(bb11) = (2*CYPA0s(m2) - 2*m1*B(bb0) - 4*f*B(bb1) 
     &	 + (p/3 - m1 - m2)*xrat)/p/6
        endif 

*** derivatives (p cannot be zero here)
        if ( abs(B(dbb11) - nan) .le. 0 )
     &	B(dbb11) = (- 2*m1*B(dbb0) - 4*f*B(dbb1) - 4*B(bb1)
     &	 + 1D0/3*xrat)/p/6 - B(bb11)/p
***]Rank 2

***[Rank 3
	  B(bb001) = .125D0*( 2*m1*B(bb1) - CYPA0s(m2) +
     &      f*(B(bb11) + 1/6D0*xrat) - .5D0*(m1 + m2)*xrat)

	if ( abs(p) .le. stab2*mmax ) then
        call CYPBnx(3, para, B(bb111))
	else
	  B(bb111) = -.25D0*(CYPA0s(m2) + 2*m1*B(bb1) +
     &      f*(3*B(bb11) + 1/6D0*xrat) - .5D0*(m1 + m2)*xrat)/p
	endif
***]Rank 3

***[Rank 4
      B(bb0011)=CYPA0s(m2)+2*m1*B(bb11)+f*B(bb111)
     & -(3*p-5*m1-15*m2)/30D0*xrat
      B(bb0011)=B(bb0011)/10D0
      
      B(bb0000)=CYPA00s(m2)+2*m1*B(bb00)+f*B(bb001)
     > +(p**2-5*p*(m1+m2)+10*(m1**2+m1*m2+m2**2))/60D0*xrat
      B(bb0000)=B(bb0000)/10D0

	if ( abs(p) .le. stab2*mmax ) then
        call CYPBnx(4, para, B(bb1111))
	else
      B(bb1111)=(CYPA0s(m2)-f*B(bb111)-6*B(bb0011))/p/2D0
	endif
***]Rank 4

***[Rank 5
      B(bb00001)=-CYPA00s(m2)+2*m1*B(bb001)+f*B(bb0011)
     > -(p*p-4*p*m1-6*p*m2+5*m1**2+10*m1*m2+15*m2**2)/120D0*xrat
      B(bb00001)=B(bb00001)/12D0
      
      B(bb00111)=-CYPA0s(m2)+2*m1*B(bb111)+f*B(bb1111)
     > +(2*p-3*m1-12*m2)/30D0*xrat
      B(bb00111)=B(bb00111)/12D0

	if ( abs(p) .le. stab2*mmax ) then
        call CYPBnx(5, para, B(bb11111))
	else
      B(bb11111)=8*B(bb00111) + B(bb1111)*f + CYPA0s(m2)
      B(bb11111)=-B(bb11111)/p/2D0
	endif
***]Rank 5

	end

************************************************************************
        subroutine CYPBnx(n, para, fx)
************************************************************************ 
* Def. CYPBnx = \int_0^1 (-x)^n [delta + ln(mu) - ln(M2)] dx
* with M2 = -p*x*(1-x) + m1*(1-x) + m2*x - i*eps, n >= 0.   
* Ref. Denner and Dittmaier hep-ph/0509141
************************************************************************ 
	implicit none
	integer i, n
	CVAR para(Pbb)
	CVAR p, m1, m2, fx
	CVAR x(2), cx(2), ctemp
	RVAR ix(2), itemp
        CVAR CYPFnx
	external CYPFnx

#include "lt.h"

c      print*,'from CYPFnx',n, x
      
	p = para(1)
	m1 = para(2)
	m2 = para(3)      
      
      fx = 0
      
***[singular cases
* B(bb11..1) is soft finite. It can be collinear divergence 
* if p = m1 = m2 = 0.
* B0 is UV divergent.

	if( abs(p) + abs(m1) + abs(m2) .le. 0 ) then
        if( lambda .ge. 0 ) then
* UV
	 fx = delta*(-1)**n/(1+n)

        else 
* collinear
	 if( lambda .eq. -1 ) fx = -1D0*(-1)**n/(1+n)
        endif 

	return
	endif
***]singular cases

***[regular cases
	if ( lambda .lt. 0 ) return

* UV part
        fx = delta

***[special cases

***[two zero masses
	if( abs(m1) + abs(m2) .le. 0 ) then
	x(1) = 1
        fx = fx - CYPLOG(-p/mudim,-1D0) + 1D0/(1 + n) - CYPFnx(n,x(1),1)
	fx = fx*(-1)**n/(1+n)
	return
	endif
***]two zero masses

***[one zero mass
* m1 eq 0
	if( abs(m1) .le. 0 ) then
	 if ( abs(p - m2) .le. 0 ) then
	  fx = fx + 2D0/(1+n) 
	  fx = fx - CYPLOG(p/mudim, -1D0)
	  fx = fx*(-1)**n/(1+n)
	  return
	 endif

	 if ( abs(p) .le. 0 ) then
	  fx = fx + 1D0/(1+n)
	  fx = fx - CYPLOG(m2/mudim, -1D0)
	  fx = fx*(-1)**n/(1+n)
	  return
	 endif

***[m2 ne 0, p ne 0 and m2 ne p
	  fx = fx + 1D0/(1+n) - CYPLOG(p/mudim, -1D0)
	  x(1) = (p - m2)/p
	  ix(1) = 1
	  if ( DBLE(p) .lt. 0 ) ix(1) = -ix(1)
	  fx = fx - CYPLOG(-x(1), -ix(1)) - CYPFnx(n,x(1),ix(1))
	  fx = fx*(-1)**n/(1+n)
	  return
***]m2 ne 0, p ne 0 and m2 ne p

	endif

* m2 eq 0
	if( abs(m2) .le. 0 ) then

	 if ( abs(p) .le. 0 ) then
	  x(1) = 1
	  fx = fx - CYPLOG(m1/mudim, -1D0) - CYPFnx(n,x(1),1)
	  fx = fx*(-1)**n/(1+n)
	  return
	 endif

	  x(1) = 1
	  fx = fx - CYPLOG(-p/mudim, -1D0) - CYPFnx(n,x(1),1)
	  x(1) = m1/p
	  ix(1) = -1
	  if ( DBLE(p) .lt. 0 ) ix(1) = -ix(1)
	  fx = fx - CYPLOG(-x(1), -ix(1)) - CYPFnx(n,x(1),ix(1))
	  fx = fx*(-1)**n/(1+n)
	  return
	endif
***]one zero mass

***[zero momentum, mi ne 0
	if ( abs(p) .le. 0 ) then
	 if ( abs(m2 - m1) .le. 0 ) then
	  fx = fx - CYPLOG(m2/mudim, -1D0)
	  fx = fx*(-1)**n/(1+n)
	  return
	 endif
	
	x(1) = m1/(m1 - m2)
	ix(1) = -1
	if ( DBLE(m1 - m2) .lt. 0 ) ix(1) = -ix(1)
	fx = fx - CYPLOG((m2-m1)/mudim, -1D0) - CYPLOG(-x(1), -ix(1))
     &   - CYPFnx(n,x(1),ix(1))
        fx = fx*(-1)**n/(1+n)
	return 
	endif
***]zero momentum, mi ne 0

***]special cases

***[general case: p ne 0 and mi ne 0
	ctemp = 1D0
	itemp = -1
	if ( DBLE(p) .lt. 0 ) itemp = 1 
	call CYPRoot2(ctemp, (m2-m1-p)/p, m1/p, itemp, x, ix)

	fx = fx - CYPLOG(p/mudim, -1D0) 

#ifdef DEBUG_Bcoeffs
	print*,'Bcoeffs.F: roots x1, x2 =', x
#endif

	do i = 1,2
	cx(i) = x(i)
	if ( DIMAG(x(i)) .eq. 0 ) cx(i) = x(i) + cI*ix(i)
	enddo
	ctemp = m2/p
	if ( DIMAG(ctemp) .eq. 0 ) ctemp = ctemp - cI/p

        fx = fx - CYPetax(-cx(1), -cx(2), ctemp)
	do i = 1, 2
	fx = fx - CYPLOG(-x(i), -ix(i)) - CYPFnx(n,x(i),ix(i)) 
	enddo
	fx = fx*(-1)**n/(1+n)
***]general case: p ne 0 and mi ne 0

***]regular cases
	end

************************************************************************
      CVAR function CYPFnx(n, x, is)
************************************************************************ 
* Def. CYPFnx = (n+1)\int_0^1 t^n ln(1 - t/x) dt, n >= 0.
* Ref. Denner and Dittmaier hep-ph/0509141
************************************************************************ 
      CVAR x, ctemp
      RVAR is
      integer n, l
      
#include "lt.h"

c      print*,'from CYPFnx',n, prec
      
      CYPFnx = 0
      
      if (abs(x) .le. 0) then 
      print*,'Error in CYPFnx(): x vanishes'
      stop
      endif
      
      if (abs(x) .lt. 10D0) then 
      
      if (abs(x - 1) .le. 0) then 
       CYPFnx = 0
      else
       CYPFnx = (1-x**(n+1))*CYPLOG(1-1/x,is)
      endif
      
      do l = 0,n
      CYPFnx = CYPFnx - x**(n-l)/(l+1) 
      enddo
      
      else
      CYPFnx = CYPLOG(1-1/x,is)
      
      ctemp = 1D0 
      l = n
      do while ( abs(ctemp) .gt. prec )
      l = l +1 

      ctemp = 1/x**(l-n)/(l+1) 
      CYPFnx = CYPFnx + ctemp
      enddo      
      
      endif

      end

************************************************************************
        subroutine CYPdBnx(n, para, fx)
************************************************************************ 
* Def. CYPdBnx = (-1)^n \int_0^1 x^(n+1)*(1-x)/M2 dx
* with M2 = -p*x*(1-x) + m1*(1-x) + m2*x - i*eps, n >= 0.   
************************************************************************ 
	implicit none
	integer i, n
	CVAR para(Pbb)
	CVAR p, m1, m2, fx
	CVAR x(2)
	RVAR ix(2), itemp
        CVAR CYPHnx, CYPLnx
	external CYPHnx, CYPLnx

#include "lt.h"

      	p = para(1)
	m1 = para(2)
	m2 = para(3)      
      
        fx = 0

***[vanishing scales
	if( abs(p) + abs(m1) + abs(m2) .le. 0 ) return 
***]vanishing scales
      
***[singular cases
***[soft
	if( n .eq. 0 ) then
	if( abs(p - m2) + abs(m1) .le. 0 .or. 
     &	 abs(p - m1) + abs(m2) .le. 0 ) then
*** dimensional regularization 
	 if( lambda .lt. 0 ) then
	  if( lambda .eq. -1 ) fx = -1/(2*p)
          return
	 endif  

*** the finite part
	 fx = -(-log(p) + 2)/2/p

*** dimensional regularization 
	 if( lambda .eq. 0 ) then
          fx = fx - log(mudim)/2/p

*** mass regularization
	 else
          fx = fx - log(lambda)/2/p
	 endif 
         
	return
	endif
	endif
***]soft
***]singular cases      

***[regular cases
	if ( lambda .lt. 0 ) return

***[special cases

***[two zero masses
	if( abs(m1) + abs(m2) .le. 0 ) then
	x(1) = 1
        fx = CYPHnx(n - 1,x(1),1D0)*(-1)**n/p
	return
	endif
***]two zero masses

***[one zero mass
* m1 eq 0
	if( abs(m1) .le. 0 ) then
	 if ( abs(p - m2) .le. 0 ) then
	  fx = CYPHnx(n - 1,czero,1D0)*(-1)**n/p
	  return
	 endif

	 if ( abs(p) .le. 0 ) then
	  fx = CYPHnx(n,czero,1D0)*(-1)**n/m2
	  return
	 endif

***[m2 ne 0, p ne 0 and m2 ne p
	  x(1) = (p - m2)/p
	  ix(1) = 1
	  if ( DBLE(p) .lt. 0 ) ix(1) = -ix(1)
	  fx = CYPHnx(n - 1,x(1),ix(1))*(-1)**n/p
	  return
***]m2 ne 0, p ne 0 and m2 ne p

	endif

* m2 eq 0
	if( abs(m2) .le. 0 ) then

	 if ( abs(p) .le. 0 ) then
	  x(1) = 1
	  fx = -CYPHnx(n,x(1),1D0)*(-1)**n/m1
	  return
	 endif

	  if ( abs(p - m1) .le. 0 ) then 

* For n > 0 the results are finite. 
* We will recalculate them by using tensor reduction formulae later.
* fx = nan here is temporary.
	   if( n .eq. 0 )
     &	   print*,'Warning in CYPdBnx(): x1=x2=1 (pinch at endpoint)'
	   fx = nan

	   return
	  else 

	  x(1) = 1
	  fx = CYPHnx(n,x(1),1D0)
	  x(1) = m1/p
	  ix(1) = -1
	  if ( DBLE(p) .lt. 0 ) ix(1) = -ix(1)
	  fx = fx - CYPHnx(n,x(1),ix(1))

	  fx = fx*(-1)**n/(p - m1)
	  endif
	  return
	endif
***]one zero mass

***[zero momentum, mi ne 0
	if ( abs(p) .le. 0 ) then
	 if ( abs(m2 - m1) .le. 0 ) then
	  fx = (-1)**n/m1/(6+5*n+n**2)
	  return
	 endif
	
	x(1) = m1/(m1 - m2)
	ix(1) = -1
	if ( DBLE(m1 - m2) .lt. 0 ) ix(1) = -ix(1)
	fx = CYPHnx(n,x(1),ix(1))*(-1)**n/(m2 - m1)
	return 
	endif
***]zero momentum, mi ne 0

***]special cases

***[general case: p ne 0 and mi ne 0
	itemp = -1
	if ( DBLE(p) .lt. 0 ) itemp = 1 
	call CYPRoot2(cone, (m2-m1-p)/p, m1/p, itemp, x, ix)

* to deal with numerical cancellation
        if ( abs(sqrt(p)-sqrt(m1)-sqrt(m2)) .le. 0 .or. 
     &   abs(sqrt(p) - abs(sqrt(m1)-sqrt(m2))) .le. 0 ) then
         x(1) = -(m2-m1-p)/p/2
	 x(2) = x(1)
	endif 

c	print*,'Bcoeffs.F x1, x2 =', x(1), x(2)


	if ( abs(x(1) - x(2)) .le. 0 ) then 
         if ( DIMAG(x(1)) .ne. 0 .or. DBLE(x(1)) .lt. 0 
     &	  .or. DBLE(x(1)) .gt. 1 ) then
	   fx = (-1)**n*CYPLnx(n,x(1),ix(1))/p 
         else   

* For n > 0 the results are finite. 
* We will recalculate them by using tensor reduction formulae later.
* fx = nan here is temporary.
	   if( n .eq. 0 )
     &	   print*,'Warning in CYPdBnx(): Real 0<x1=x2<1 (pinch)'
	   fx = nan

	   return
         endif 

	else 
	 do i = 1, 2
	 fx = fx - (-1)**i*CYPHnx(n,x(i),ix(i)) 
	 enddo
	 fx = fx*(-1)**n/p/(x(1) - x(2))
	endif

***]general case: p ne 0 and mi ne 0

***]regular cases
	end

************************************************************************
      CVAR function CYPHnx(n, x, is)
************************************************************************ 
* Def. CYPHnx = \int_0^1 t^(n+1)*(1-t)/(t - x) dt, n >= -1.
************************************************************************ 
      CVAR x, CYPGnx
      RVAR is
      integer n
      external CYPGnx
      
#include "lt.h"

c      print*,'from CYPHnx',n
      
      CYPHnx = 0
      
      if (n .lt. -1) then 
       print*,'Error in CYPHnx(): invalid input (n < -1)'
       stop
      endif            
      
      if (abs(x) .le. 0 .and. n .gt. -1) then 
       CYPHnx = 1D0/(2 + 3*n + n**2)
       return
      endif

      if (abs(1 - x) .le. 0) then 
       CYPHnx = -1D0/(2 + n)
       return
      endif

      CYPHnx = CYPGnx(n + 1, x, is) - CYPGnx(n + 2, x, is)

      end

************************************************************************
      CVAR function CYPLnx(n, x, is)
************************************************************************ 
* Def. CYPLnx = \int_0^1 t^(n+1)*(1-t)/(t - x)^2 dt, n >= -1.
************************************************************************ 
      CVAR x, CYPKnx, CYPGnx
      RVAR is
      integer n
      external CYPKnx, CYPGnx
      
#include "lt.h"

c      print*,'from CYPLnx',n
      
      CYPLnx = 0

      if (n .lt. -1) then 
       print*,'Error in CYPLnx(): invalid input (n < -1)'
       stop
      endif            
      
      if (abs(x) .le. 0 .and. n .gt. 0) then 
       CYPLnx = 1D0/(n + n**2)
       return
      endif

      if (abs(1 - x) .le. 0) then 
       CYPLnx = -CYPGnx(n + 1, x, is)
       return
      endif

      CYPLnx = CYPKnx(n + 1, x, is) - CYPKnx(n + 2, x, is)

      end

************************************************************************
      CVAR function CYPGnx(n, x, is)
************************************************************************ 
* Def. CYPGnx = \int_0^1 t^n/(t - x) dt, n >= 0.
************************************************************************ 
      CVAR x, ctemp
      RVAR is, k1, nk1
      integer n, k, j
      
#include "lt.h"

c      print*,'from CYPGnx',n
      
      CYPGnx = 0
      
      if (n .lt. 0) then 
       print*,'Error in CYPGnx(): invalid input (n < 0)'
       stop
      endif            
      
      if (abs(x) .le. 0 .and. n .gt. 0) then 
       CYPGnx = 1D0/n
       return
      endif      
      
      if (abs(x*(1 - x)) .le. 0) then 
      print*,'Error in CYPGnx(): invalid input [x(1-x)=0]'
      stop
      endif

      if ( abs(x) .gt. 10D0 ) then
      ctemp = 1D0
      k = 0
      do while ( abs(ctemp) .gt. prec )
      k = k + 1 
      ctemp = 1D0/(k + n)/x**k
      CYPGnx = CYPGnx - ctemp
      enddo

      return
      endif

      CYPGnx = CYPLOG(1 - 1/x, is)*x**n 

      if (n .gt. 0) then 
      do k = 1, n
       k1 = 1
       nk1 = 1
       do j = 1, k
        k1 = k1*j
	nk1 = nk1*(n - k + j)
       enddo

      CYPGnx = CYPGnx + nk1/k1*x**(n-k)*((1 - x)**k - (-x)**k)/k
      enddo
      endif 

      end

************************************************************************
      CVAR function CYPKnx(n, x, is)
************************************************************************ 
* Def. CYPKnx = \int_0^1 t^n/(t - x)^2 dt, n >= 0.
************************************************************************ 
      CVAR x, ctemp
      RVAR is, k1, nk1
      integer n, k, j
      
#include "lt.h"

c      print*,'from CYPKnx',n
      
      CYPKnx = 0
      
      if (n .lt. 0) then 
       print*,'Error in CYPKnx(): invalid input (n < 0)'
       stop
      endif            
      
      if (abs(x) .le. 0 .and. n .gt. 1) then 
       CYPKnx = 1D0/(n-1)
       return
      endif      
      
      if (abs(x*(1 - x)) .le. 0) then 
      print*,'Error in CYPKnx(): invalid input [x(1-x)=0]'
      stop
      endif

      if ( abs(x) .gt. 10D0 ) then
      ctemp = 1D0
      k = 0
      do while ( abs(ctemp) .gt. prec )
      k = k + 1 
      ctemp = k/x**(k + 1)/(k + n)
      CYPKnx = CYPKnx + ctemp
      enddo

      return
      endif

      CYPKnx = -x**n/x/(1-x) + n*CYPLOG(1 - 1/x, is)*x**(n-1) 

      if (n .gt. 1) then 
      do k = 2, n
       k1 = 1
       nk1 = 1
       do j = 1, k
        k1 = k1*j
	nk1 = nk1*(n - k + j)
       enddo

      CYPKnx = CYPKnx + 
     & nk1/k1*x**(n-k)*((1 - x)**(k-1) - (-x)**(k-1))/(k-1)
      enddo
      endif

      end
