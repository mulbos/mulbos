* D0div.F
#include "defs.h"

***********************************************************************
* Calculate the infrared (soft and collinear) scalar 4-point integrals.
* Ref0. Denner and Dittmaier arXiv: 1005.2076
* Ref1. Beenakker and Denner Nucl. Phys. B338, 349
* Ref2. Ellis and Zanderighi arXiv: 0712.1851v3 
* Convention: same as Ref0 for mass regularization (MR), 
* same as Ref2 for dimensional regularization (DR). 
* Translation rules from Ref0 to Ref2 for DR: 
* Ref2 = z0 + z1/eps + z2/eps^2,
* Ref0 = d0 + d1/eps + d2/eps^2,
* then: z2 = d2, z1 = d1, z0 = d0 + d2*pi^2/6. 
* Input: q1,q2,q3,q4,q1q2,q2q3,n1,n2,n3,n4 and lambda. 
* Output: D0, fail

* todo: many singular cases to be done.
* Author: Le Duc Ninh

***********************************************************************
	subroutine CYPD0ir_4m0_4p0a(fail,p1p2,p2p3,D0)
***********************************************************************
* Ref0. eq. (4.71)
* Author: Le Duc Ninh
***********************************************************************
	implicit none
	CVAR D0
        CVAR p1p2, p2p3
        CVAR rmu, rla
	integer fail
#include "lt.h"

	D0 = 0
	fail = 0

	rmu = sqrt(mudim)
	if ( lambda .gt. 0 ) then
	 rla = sqrt(lambda)
	else
	 rla = lambda
	endif
***[4-soft 4-coll
	 if( lambda .lt. 0 ) then
	  if( lambda .eq. -2 ) D0 = 4/(p1p2*p2p3)
	  if( lambda .eq. -1 ) D0 = 2/(p1p2*p2p3)*
     &	  (CYPLOG(-mudim/p1p2,1D0) + CYPLOG(-mudim/p2p3,1D0))
	  return
	 endif	

	 D0 = (2*CYPLOG(-mudim/p2p3,1D0)*CYPLOG(-mudim/p1p2,1D0) -
     &   pi**2)/(p1p2*p2p3)
***]4-soft 4-coll
	end
***********************************************************************
	subroutine CYPD0ir_3m0_2p0a(fail,p3,p1p2,p2p3,D0)
***********************************************************************
* Ref0. eq. (4.45)
* Author: Le Duc Ninh
***********************************************************************
	implicit none
	CVAR D0
        CVAR p3, p1p2, p2p3, m4, rm4
        CVAR rmu, rla
	integer fail
#include "lt.h"

	D0 = 0
	fail = 0

	rmu = sqrt(mudim)
	if ( lambda .gt. 0 ) then
	 rla = sqrt(lambda)
	else
	 rla = lambda
	endif

***[3-soft 2-coll
	 m4 = p3
	 rm4 = sqrt(m4)

	 if( lambda .lt. 0 ) then
	  if( lambda .eq. -2 ) D0 = 2/(p1p2*(p2p3 - m4))
	  if( lambda .eq. -1 ) D0 = 1/(p1p2*(p2p3 - m4))*
     &	  (CYPLOG(-mudim/p1p2,1D0) + 2*CYPLOG(rmu*rm4/(m4-p2p3),1D0))
	  return
	 endif	
	
	 D0 = (2*CYPLOG(-mudim/p1p2,1D0)*CYPLOG(rmu*rm4/(m4-p2p3),1D0) -
     &   pi**2/2D0)/(p1p2*(p2p3 - m4))
***]3-soft 2-coll
	end
***********************************************************************
	subroutine CYPD0ir_4m0_3p0a(fail,p4,p1p2,p2p3,D0)
***********************************************************************
* Ref0. eq. (4.60)
* Author: Le Duc Ninh
***********************************************************************
	implicit none
	CVAR D0
        CVAR p4, p1p2, p2p3
        CVAR rmu, rla
* These are temporary variables. 
	CVAR x(2)
	RVAR ix(2)

	integer fail
#include "lt.h"

	D0 = 0
	fail = 0

	rmu = sqrt(mudim)
	if ( lambda .gt. 0 ) then
	 rla = sqrt(lambda)
	else
	 rla = lambda
	endif

***[2-soft 3-coll
c	print*,'Ref0. (4.60)'

	x(1) = p4/p2p3
	ix(1) = DBLE(p2p3 - p4)
	if ( ix(1) .eq. 0 ) ix(1) = 1D0

	x(2) = p4/p1p2
	ix(2) = DBLE(p1p2 - p4)
	if ( ix(2) .eq. 0 ) ix(2) = 1D0

	if( lambda .lt. 0 ) then
	 if( lambda .eq. -2 ) D0 = 2/(p1p2*p2p3)
	 if( lambda .eq. -1 ) D0 = 
     &    2*(CYPLOG(-mudim/p1p2,1D0)+CYPLOG(x(1),ix(1)))
     &	  /(p1p2*p2p3)
	 return
	endif

	D0 = 2*CYPLOG(-mudim/p2p3,1D0)*CYPLOG(-mudim/p1p2,1D0) - 
     &  CYPLOG(-mudim/p4,1D0)**2 - 2*CYPSp(1-x(1),-ix(1)) - 
     &  2*CYPSp(1-x(2),-ix(2)) - 2*pi**2/3
	D0 = (D0 + 2*pi**2/6)/(p1p2*p2p3)
***]2-soft 3-coll
	end
***********************************************************************
	subroutine CYPD0ir_3m0_2p0b(fail,p3,p4,p1p2,p2p3,D0)
***********************************************************************
* Ref0. eq. (4.36)
* Author: Le Duc Ninh
***********************************************************************
	implicit none
	CVAR D0
        CVAR p3, p4, p1p2, p2p3, m4, rm4
        CVAR rmu, rla
* These are temporary variables. 
	CVAR x(1)
	RVAR ix(1)

	integer fail
#include "lt.h"

	D0 = 0
	fail = 0

	rmu = sqrt(mudim)
	if ( lambda .gt. 0 ) then
	 rla = sqrt(lambda)
	else
	 rla = lambda
	endif

***[2-soft 2-coll
	m4 = p4
	rm4 = sqrt(m4)

* only DR possible

c	print*,'Ref0. eq. (4.36)'

	 if( lambda .lt. 0 ) then
	  if( lambda .eq. -2 ) D0 = 3/2D0/p1p2/(p2p3 - m4)
	  if( lambda .eq. -1 ) D0 = 1/p1p2/(p2p3 - m4)*(
     &     CYPLOG(-mudim/p1p2,1D0) - CYPLOG(rmu*rm4/(m4 - p3),1D0) 
     &     + 2*CYPLOG(rmu*rm4/(m4 - p2p3),1D0))
	  return
	 endif

	ix(1) = DBLE(p3 - p2p3)
	if ( ix(1) .eq. 0 ) ix(1) = 1D0

	D0 = -CYPLOG(rmu*rm4/(m4 - p3),1D0)**2 
     &  + 2*CYPLOG(-mudim/p1p2,1D0)*CYPLOG(rmu*rm4/(m4 - p2p3),1D0) 
     &  - 2*CYPSp((p3 - p2p3)/(m4 - p2p3), ix(1)) - 2*pi**2/3D0 
     &  + pi**2/4D0
        D0 = D0/p1p2/(p2p3 - m4)
***]2-soft 2-coll
	end
***********************************************************************
	subroutine CYPD0ir_2m0_1p0a(fail,p2,p3,p4,p1p2,p2p3,D0)
***********************************************************************
* Ref0. eq. (4.25)
* Author: Le Duc Ninh
***********************************************************************
	implicit none
	CVAR D0
        CVAR p2, p3, p4, p1p2, p2p3, m3, rm3, m4, rm4
        CVAR rmu, rla
	CVAR x34
	RVAR ix34, ix34s

	integer fail
#include "lt.h"

	D0 = 0
	fail = 0

	rmu = sqrt(mudim)
	if ( lambda .gt. 0 ) then
	 rla = sqrt(lambda)
	else
	 rla = lambda
	endif

***[2-soft 1-coll
	m3 = p2
	m4 = p4
	rm3 = sqrt(m3)
	rm4 = sqrt(m4)

	 if( lambda .lt. 0 ) then
	  if( lambda .eq. -2 ) D0 = 1/((p1p2 - m3)*(p2p3 - m4))
	  if( lambda .eq. -1 ) D0 = 1/((p1p2 - m3)*(p2p3 - m4))*(
     &     CYPLOG(rmu*rm3/(m3 - p1p2),1D0) + 
     &    CYPLOG(rmu*rm4/(m4 - p2p3),1D0))
	  return
	 endif

	call CYPKam(p3, rm3, rm4, x34, ix34, ix34s)

	D0 = 2*CYPLOG(rmu*rm3/(m3 - p1p2),1D0)*CYPLOG(rmu*rm4/(m4 - p2p3),1D0) 
     &  - CYPLOG(x34,ix34)**2 - 2*pi**2/3

        D0 = (D0 + pi**2/6)/((p1p2 - m3)*(p2p3 - m4))
***]2-soft 1-coll
	end
***********************************************************************
	subroutine CYPD0ir_2m0_0p0a(fail,p1,p3,p1p2,p2p3,D0)
***********************************************************************
* Ref0. eq. (4.9) and (4.11) 
* Author: Le Duc Ninh
***********************************************************************
	implicit none
	CVAR D0
        CVAR p1, p3, p1p2, p2p3, m2, rm2, m4, rm4
        CVAR rmu, rla
	CVAR x24
	RVAR ix24, ix24s
	integer fail
#include "lt.h"

	D0 = 0
	fail = 0

	rmu = sqrt(mudim)
	if ( lambda .gt. 0 ) then
	 rla = sqrt(lambda)
	else
	 rla = lambda
	endif

***[2-soft 0-coll
	m2 = p1
	m4 = p3
	rm2 = sqrt(m2)
	rm4 = sqrt(m4)

	call CYPKam(p2p3, rm2, rm4, x24, ix24, ix24s)

	if ( abs(x24 - 1) .le. 0 ) then
	 if( lambda .lt. 0 ) then
	  if( lambda .eq. -1 ) D0 = 1/(rm2*rm4*p1p2)
	  return
	 endif

	 D0 = -CYPLOG(-p1p2, -1D0)/(rm2*rm4*p1p2)
	 if( lambda .eq. 0 ) then
          D0 = D0 + log(mudim)/(rm2*rm4*p1p2)
	 else
          D0 = D0 + log(lambda)/(rm2*rm4*p1p2)
	 endif

	else
	 if( lambda .lt. 0 ) then
	  if( lambda .eq. -1 ) D0 = -2*x24*CYPLOG(x24,ix24)/
     &	    (rm2*rm4*p1p2*(1 - x24**2))	   
	  return
	 endif

	 D0 = CYPLOG(-p1p2, -1D0)*2*x24*CYPLOG(x24,ix24)/
     &	  (rm2*rm4*p1p2*(1 - x24**2))
	 if( lambda .eq. 0 ) then
          D0 = D0 - log(mudim)*2*x24*CYPLOG(x24,ix24)/
     &	  (rm2*rm4*p1p2*(1 - x24**2))
	 else
          D0 = D0 - log(lambda)*2*x24*CYPLOG(x24,ix24)/
     &	  (rm2*rm4*p1p2*(1 - x24**2))
	 endif
	
	endif
***]2-soft 0-coll
	end
***********************************************************************
	subroutine CYPD0ir_3m0_2p0c(fail,p3,p4,p1p2,p2p3,m4,D0)
***********************************************************************
* Ref0. eq. (4.29)
* Author: Le Duc Ninh
***********************************************************************
	implicit none
	CVAR D0
        CVAR p3, p4, p2p3, p1p2, m4, rm4
        CVAR rmu, rla
* These are temporary variables. 
	CVAR x(3)
	RVAR ix(3)
	integer fail
#include "lt.h"

	D0 = 0
	fail = 0

	rmu = sqrt(mudim)
	if ( lambda .gt. 0 ) then
	 rla = sqrt(lambda)
	else
	 rla = lambda
	endif

***[1-soft 2-coll
	rm4 = sqrt(m4)

	 ix(1) = DBLE(p2p3 - p3)
	 if ( ix(1) .eq. 0 ) ix(1) = 1D0

	 ix(2) = DBLE(p2p3 - p4)
	 if ( ix(2) .eq. 0 ) ix(2) = 1D0

	 ix(3) = DBLE(p4 - p3)
	 if ( ix(3) .eq. 0 ) ix(3) = 1D0

	 if( lambda .lt. 0 ) then
	  if( lambda .eq. -2 ) D0 = 1/(p1p2*(p2p3 - m4))
	  if( lambda .eq. -1 ) D0 = 1/(p1p2*(p2p3 - m4))*
     &	  (CYPLOG(-mudim/p1p2,1D0) + CYPLOG((m4-p3)/(m4-p2p3),ix(1)) + 
     &    CYPLOG((m4-p4)/(m4-p2p3),ix(2)))
	  return
	 endif

	 D0 = CYPLOG(-mudim/p1p2,1D0)**2/2 + CYPLOG(-mudim/p1p2,1D0)*
     &   (CYPLOG((m4-p3)/(m4-p2p3),ix(1)) + 
     &    CYPLOG((m4-p4)/(m4-p2p3),ix(2))) - 
     &    CYPLOG((m4-p3)/(m4-p4),ix(3))**2/2 
 
         D0 = D0 - 
     &   2*CYPSp((p3 - p2p3)/(m4 - p2p3), -ix(1)) - 
     &   2*CYPSp((p4 - p2p3)/(m4 - p2p3), -ix(2))	 

	x(1) = p1p2/(p3 - m4)
	ix(1) = DBLE(p3 - m4 - p1p2)
	if ( ix(1) .eq. 0 ) ix(1) = 1D0

	x(2) = m4/(m4 - p4)
	ix(2) = 1D0
	D0 = D0 + CYPLi2n(x, ix, 2) - pi**2/6
	D0 = D0/p1p2/(p2p3-m4)
***]1-soft 2-coll
	end
***********************************************************************
	subroutine CYPD0ir_2m0_1p0b(fail,p2,p3,p4,p1p2,p2p3,m3,D0)
***********************************************************************
* Ref0. eq. (4.19) 
* Author: Le Duc Ninh
***********************************************************************
	implicit none
	CVAR D0
        CVAR p2, p3, p4, p1p2, p2p3, m3, rm3, m4, rm4
        CVAR rmu, rla
	CVAR x34
	RVAR ix34, ix34s
* These are temporary variables. 
	CVAR x(2)
	RVAR ix(2)
	integer fail, l
#include "lt.h"

	D0 = 0
	fail = 0

	rmu = sqrt(mudim)
	if ( lambda .gt. 0 ) then
	 rla = sqrt(lambda)
	else
	 rla = lambda
	endif

***[1-soft 1-coll
	rm3 = sqrt(m3)
	m4 = p4
	rm4 = sqrt(m4)

	call CYPKam(p3, rm3, rm4, x34, ix34, ix34s)
	ix(1) = DBLE(p1p2 - p2)
	if ( ix(1) .eq. 0 ) ix(1) = 1D0

	 if( lambda .lt. 0 ) then
	  if( lambda .eq. -2 ) D0 = 1/(p2p3 - m4)/(p1p2-m3)/2
	  if( lambda .eq. -1 ) D0 = (CYPLOG(rm4*rmu/(m4 - p2p3),1D0) + 
     &    CYPLOG((m3 - p2)/(m3 - p1p2),ix(1)))
     &	  /(p2p3 - m4)/(p1p2-m3)
	  return
	 endif

	D0 = CYPLOG(rm4*rmu/(m4 - p2p3),1D0)*
     &	2*CYPLOG((m3 - p2)/(m3 - p1p2),ix(1)) + 
     &  CYPLOG(rm4*rmu/(m4 - p2p3),1D0)**2 - pi**2/12

	D0 = D0 - 2*CYPSp((p2 - p1p2)/(m3 - p1p2), -ix(1))

        do l = 1, -1, -2
	x(1) = (m4 - p2p3)/(m3 - p2)
	ix(1) = DBLE(m4 - p2p3) - DBLE(m3 - p2) 
	if ( ix(1) .eq. 0 ) ix(1) = 1D0
	if ( abs(m3) .le. 0 ) then
	 x(2) = (m4 - p3)/m4
	 x(2) = -x(2)*(l-1)/2
	 ix(2) = -1D0
	else
	 x(2) = x34**l*rm3/rm4
	ix(2) = ix34*l 
	endif
	D0 = D0 + CYPLi2n(x, ix, 2)
	enddo

	D0 = D0/(p1p2 - m3)/(p2p3 - m4)
***]1-soft 1-coll
	end
***********************************************************************
	subroutine CYPD0ir_1m0_0p0a(fail,p1,p2,p3,p4,p1p2,p2p3,m3,D0)
***********************************************************************
* Ref0. eqs. (4.4) and (4.5)
* Author: Le Duc Ninh
***********************************************************************
	implicit none
	CVAR D0
        CVAR p1, p2, p3, p4, p1p2, p2p3, m2, rm2, m3, rm3, m4, rm4
        CVAR rmu, rla
* These are not temporary variables. 
	CVAR x23, x24, x34
	RVAR ix23, ix24, ix34, ix23s, ix24s, ix34s
* These are temporary variables. 
	CVAR x(3)
	RVAR ix(3)
	integer fail, l, k
#include "lt.h"

	D0 = 0
	fail = 0

	rmu = sqrt(mudim)
	if ( lambda .gt. 0 ) then
	 rla = sqrt(lambda)
	else
	 rla = lambda
	endif

***[1-soft 0-coll
	m2 = p1
	m4 = p4
	rm2 = sqrt(m2)
	rm3 = sqrt(m3)
	rm4 = sqrt(m4)

	call CYPKam(p2p3, rm2, rm4, x24, ix24, ix24s)
	call CYPKam(p2, rm2, rm3, x23, ix23, ix23s)
	call CYPKam(p3, rm3, rm4, x34, ix34, ix34s)
       
*** dimensional regularization 
	 if( lambda .lt. 0 ) then
	  if( lambda .eq. -1 ) then
	   if ( abs(x24 - 1) .le. 0 ) then
	    D0 = 1/(rm2*rm4*(p1p2 - m3)*2)
	   else
	    D0 = -x24*CYPLOG(x24,ix24)/
     &	    (rm2*rm4*(p1p2 - m3)*(1 - x24**2))	   
           endif
	  endif
          return
	 endif 

***[the finite part
	if ( abs(m3) .le. 0 ) then
	x(2) = rm2*(m4 - p3)/rm4/(m2 - p2)
	ix(2) = -DBLE(m2 - p2) + DBLE(m4 - p3)
        if ( ix(2) .eq. 0 ) ix(2) = 1D0

***[case x24 = 1 (or p2p3 = (rm2 - rm4)**2 )
	if ( abs(x24 - 1) .le. 0 ) then
	D0 = log(rm2*rm4) - CYPLOG(m2 - p2, -1D0) 
     &	+ 2*CYPLOG(-p1p2, -1D0) - CYPLOG(m4 - p3, -1D0) 
     &  + 2 + CYPLOG(x(2), ix(2))*(1+x(2))/(1-x(2))
	D0 = -D0/(rm2*rm4*(p1p2 - m3)*2)

c	print*,'Check D0 02, x(2), x23, x24, x34 =', x(2), x23, x24, x34

*** dimensional regularization 
	 if( lambda .eq. 0 ) then
          D0 = D0 + log(mudim)/(rm2*rm4*(p1p2 - m3)*2)
*** mass regularization
	 else
          D0 = D0 + log(lambda)/(rm2*rm4*(p1p2 - m3)*2)
	 endif  

	return
	endif
***]case x24 = 1 (or p2p3 = (rm2 - rm4)**2 )

	D0 = CYPLOG(x24,ix24)*(-CYPLOG(x24,ix24)/2 + log(rm2*rm4) - 
     &  CYPLOG(m2 - p2, -1D0) + 2*CYPLOG(-p1p2, -1D0) - 
     &  CYPLOG(m4 - p3, -1D0)) + 0.5D0*CYPLOG(x(2), ix(2))**2

	x(1) = x24
	ix(1) = ix24
        do l = 1, -1, -2
	x(2) = x(2)**l
	ix(2) = ix(2)*l
	
c	x(3) = 1D0
c	ix(3) = 0

	D0 = D0 + CYPLi2n(x, ix, 2)
	enddo
	
	D0 = D0 + 2*CYPLOG(x24,ix24)*CYPLOG(1 - x24**2,-ix24s) + 
     &	CYPSp(x24**2,ix24s) - pi**2/6
	
	D0 = D0*x24/(rm2*rm4*(p1p2 - m3)*(1 - x24**2))	

	else
* abs(m3) > 0 here
***[case x24 = 1 (or p2p3 = (rm2 - rm4)**2 )
	if ( abs(x24 - 1) .le. 0 ) then
	D0 = 2*CYPLOG(rm3/(m3-p1p2),1D0) - 2
	x(2) = x23
	ix(2) = ix23
        do l = -1, 1, 2
	x(3) = x34**l
	ix(3) = ix34*l
	D0 = D0 - (CYPLOG(x(2), ix(2)) + CYPLOG(x(3), ix(3)))*
     &	(1 + x(2)*x(3))/(1 - x(2)*x(3)) 
	enddo
	D0 = D0/(rm2*rm4*(p1p2 - m3)*2)

c	print*,'Check D0 20, x(2), x23, x24, x34 =', x(2), x23, x24, x34

*** dimensional regularization 
	 if( lambda .eq. 0 ) then
          D0 = D0 + log(mudim)/(rm2*rm4*(p1p2 - m3)*2)
*** mass regularization
	 else
          D0 = D0 + log(lambda)/(rm2*rm4*(p1p2 - m3)*2)
	 endif  

	return
	endif
***]case x24 = 1 (or p2p3 = (rm2 - rm4)**2 )

	D0 = -2*CYPLOG(x24,ix24)*CYPLOG(rm3/(m3-p1p2),1D0) + 
     &	CYPLOG(x23,ix23)**2 + CYPLOG(x34,ix34)**2

	x(1) = x24
	ix(1) = ix24
        do k = 1, -1, -2
        do l = 1, -1, -2
	x(2) = x23**k
	x(3) = x34**l
	ix(2) = ix23*k
	ix(3) = ix34*l

	D0 = D0 + CYPLi2n(x, ix, 3)
	enddo
	enddo

	D0 = D0 + 2*CYPLOG(x24,ix24)*CYPLOG(1 - x24**2,-ix24s) + 
     &	CYPSp(x24**2,ix24s) - pi**2/6
	
	D0 = D0*x24/(rm2*rm4*(p1p2 - m3)*(1 - x24**2))

	endif
***]the finite part

*** dimensional regularization 
	 if( lambda .eq. 0 ) then
          D0 = D0 - log(mudim)*CYPLOG(x24,ix24)*x24/
     &	  (rm2*rm4*(p1p2 - m3)*(1 - x24**2))

*** mass regularization
	 else
          D0 = D0 - log(lambda)*CYPLOG(x24,ix24)*x24/
     &	  (rm2*rm4*(p1p2 - m3)*(1 - x24**2))
	 endif

	return
***]1-soft 0-coll
	end
***********************************************************************
	subroutine CYPD0ir_4m0_2p0a(fail,p2,p4,p1p2,p2p3,D0)
***********************************************************************
* Ref0. eqs. (4.52)
* Author: Le Duc Ninh
***********************************************************************
	implicit none
	CVAR D0
        CVAR p2, p4, p1p2, p2p3
        CVAR rmu, rla
* These are temporary variables. 
	CVAR x(2), cx(2)
	RVAR ix(2), icx(2)
	integer fail, l
#include "lt.h"

	D0 = 0
	fail = 0

	rmu = sqrt(mudim)
	if ( lambda .gt. 0 ) then
	 rla = sqrt(lambda)
	else
	 rla = lambda
	endif

***[2-coll (opposite)
c	print*,'Ref0. (4.52)'

	x(1) = p2/p1p2
	ix(1) = DBLE(p1p2 - p2)
	if ( ix(1) .eq. 0 ) ix(1) = 1D0

	x(2) = p4/p2p3
	ix(2) = DBLE(p2p3 - p4)
	if ( ix(2) .eq. 0 ) ix(2) = 1D0

* Landau singularity (double parton scattering sing.)
* return 0 here
	if ( abs(p1p2*p2p3 - p2*p4) .le. 0 ) return

	 if( lambda .lt. 0 ) then
	  if( lambda .eq. -1 ) D0 = 
     &    2*(CYPLOG(x(2),ix(2))+CYPLOG(x(1),ix(1)))
     &	  /(p1p2*p2p3 - p2*p4)
	  return
	 endif

	cx(1) = p1p2/p4
	icx(1) = DBLE(p4 - p1p2)
	if ( icx(1) .eq. 0 ) icx(1) = 1D0

	cx(2) = p2/p2p3
	icx(2) = DBLE(p2p3 - p2)
	if ( icx(2) .eq. 0 ) icx(2) = 1D0 

	D0 = (CYPLi2n(x, ix, 2) + 
     &  CYPLOG(-mudim/p1p2,1D0)*(CYPLOG(x(2),ix(2))+CYPLOG(x(1),ix(1)))
     &  + CYPSp(1-1/x(1),ix(1)) + CYPSp(1-cx(1),-icx(1)) - 
     &  CYPSp(1-cx(2),-icx(2)) - CYPSp(1-x(2),-ix(2)))*
     &  2/(p1p2*p2p3 - p2*p4)
***]2-coll (opposite)
	end
***********************************************************************
	subroutine CYPD0ir_2m0_1p0c(fail,p2,p3,p4,p1p2,p2p3,m3,m4,D0)
***********************************************************************
* Ref0. eqs. (4.13)
* Author: Le Duc Ninh
***********************************************************************
	implicit none
	CVAR D0
        CVAR p2, p3, p4, p1p2, p2p3, m3, m4, rm3, rm4
        CVAR rmu, rla
* These are not temporary variables. 
	CVAR x34
	RVAR ix34, ix34s
* These are temporary variables. 
	CVAR x(3), cx(4)
	RVAR ix(3), icx(4)
	integer fail, l, k
#include "lt.h"

	D0 = 0
	fail = 0

	rmu = sqrt(mudim)
	if ( lambda .gt. 0 ) then
	 rla = sqrt(lambda)
	else
	 rla = lambda
	endif

***[1-coll
	rm3 = sqrt(m3)
	rm4 = sqrt(m4)

	call CYPKam(p3, rm3, rm4, x34, ix34, ix34s)

* Ref0. (4.16)
	cx(1) = (m4 - p4)/(m3 - p1p2)
	icx(1) = DBLE(m4 - p4 - m3 + p1p2)
	if ( icx(1) .eq. 0 ) icx(1) = 1D0
	cx(2) = (m4 - p2p3)/(m3 - p2)
	icx(2) = DBLE(m4 - p2p3 - m3 + p2)
	if ( icx(2) .eq. 0 ) icx(2) = 1D0	

	if ( abs(m4) .gt. 0 ) then
	if ( abs(m3) .le. 0 ) then
	 cx(3) = (m4 - p3)/m4
	 icx(3) = -1D0
	 cx(4) = 0
	 icx(4) = 1D0
	else
	 cx(3) = (rm3/rm4)/x34
	 icx(3) = -ix34
	 cx(4) = (rm3/rm4)*x34
	 icx(4) = ix34
	endif
	
	else
	 cx(3) = nan
	 icx(3) = 0
	 cx(4) = m3/(m3 - p3)
	 icx(4) = 1D0
	endif

* Ref0. (4.13)
	 if( lambda .lt. 0 ) then
	  if( lambda .eq. -1 ) then
	   D0 = CYPLOG(cx(1),icx(1)) - CYPLOG(cx(2),icx(2))
	   D0 = D0/((p1p2 - m3)*(p2p3 - m4) - (p2 - m3)*(p4 - m4))
	  endif
          return
	 endif 

	   x(1) = (p1p2 - p2)/(m3 - p2)
	   ix(1) = DBLE(p1p2 - p2)
	   if ( ix(1) .eq. 0) ix(1) = 1D0
	   x(2) = (p4 - p2p3)/(m4 - p2p3)
	   ix(2) = DBLE(p4 - p2p3)
	   if ( ix(2) .eq. 0) ix(2) = 1D0
	   x(3) = p3/p4
	   ix(3) = DBLE(p4 - p3)
	   if ( ix(3) .eq. 0) ix(3) = 1D0
	   
           D0 = 2*(CYPSp(x(1),ix(1)) - CYPSp(x(2),ix(2)))
	   x(1) = cx(2)
	   ix(1) = icx(2)
	   x(2) = 1/cx(1)
	   ix(2) = -icx(1)
	   D0 = D0 - 2*CYPLi2n(x,ix,2)	   

	   if ( abs(m3) + abs(m4) .le. 0 ) then
* Ref0. (4.17)
	   D0 = D0 + 
     &	   (CYPLOG(cx(1),icx(1))**2 - CYPLOG(cx(2),icx(2))**2)/2 + 
     &     (CYPLOG(cx(1),icx(1)) - CYPLOG(cx(2),icx(2)))*
     &	   (CYPLOG(-mudim/p4,1D0) + CYPLOG(x(3),ix(3)))

	   else
* Ref0. (4.13)
	   D0 = D0 + 
     &	   (CYPLOG(cx(1),icx(1)) - CYPLOG(cx(2),icx(2)))*2*
     &	   CYPLOG(rmu*rm4/(m4 - p4),1D0)

	   do k = 1, 2
	    do l = 3, 4
	     x(1) = cx(k)
	     ix(1) = icx(k)
	     x(2) = cx(l)
	     ix(2) = icx(l)
	     D0 = D0 + (-1D0)**k*CYPLi2n(x,ix,2)
	    enddo
	   enddo
	   endif 
	   D0 = D0/((p1p2 - m3)*(p2p3 - m4) - (p2 - m3)*(p4 - m4))
***]1-coll
	end
***********************************************************************
	subroutine CYPD0ir_3m0q_2p0a(fail,p3,p1p2,p2p3,D0)
***********************************************************************
* Ref0. eqs. (4.72)
* Author: Le Duc Ninh
***********************************************************************
	implicit none
	CVAR D0
        CVAR p3, p1p2, p2p3, m4, rm4
        CVAR rmu, rla
	integer fail
#include "lt.h"

	D0 = 0
	fail = 0

	rmu = sqrt(mudim)
	if ( lambda .gt. 0 ) then
	 rla = sqrt(lambda)
	else
	 rla = lambda
	endif

***[3-soft 4-coll (2-quasi)
	m4 = p3
	rm4 = sqrt(m4)

	 if( lambda .lt. 0 ) then
	  if( lambda .eq. -2 ) D0 = 2/(p1p2*p2p3)
	  if( lambda .eq. -1 ) D0 = 1/(p1p2*p2p3)*
     &	  (CYPLOG(-mudim/p1p2,1D0) + 2*CYPLOG(-rmu*rm4/p2p3,1D0))
	  return
	 endif	

	 D0 = (2*CYPLOG(-rmu*rm4/p2p3,1D0)*CYPLOG(-mudim/p1p2,1D0) -
     &   pi**2/2D0)/(p1p2*p2p3)
***]3-soft 4-coll (2-quasi)
	end
***********************************************************************
	subroutine CYPD0ir_3m0q_2p0b(fail,p3,p4,p1p2,p2p3,D0)
***********************************************************************
* Ref0. eq. (4.61)
* Author: Le Duc Ninh
***********************************************************************
	implicit none
	CVAR D0
        CVAR p3, p4, p1p2, p2p3, m4, rm4
        CVAR rmu, rla
* These are temporary variables. 
	CVAR x(1)
	RVAR ix(1)

	integer fail
#include "lt.h"

	D0 = 0
	fail = 0

	rmu = sqrt(mudim)
	if ( lambda .gt. 0 ) then
	 rla = sqrt(lambda)
	else
	 rla = lambda
	endif

***[2-soft 3-coll (1-quasi)
	m4 = p4
	rm4 = sqrt(m4)

	if( lambda .lt. 0 ) then
	 if( lambda .eq. -2 ) D0 = 3/(2*p1p2*p2p3)
	 if( lambda .eq. -1 ) D0 = 
     &    (2*CYPLOG(-rmu*rm4/p2p3,1D0)-CYPLOG(-rmu*rm4/p3,1D0) + 
     &    CYPLOG(-mudim/p1p2,1D0))/(p1p2*p2p3)
	 return
	endif

	x(1) = p3/p2p3
	ix(1) = DBLE(p2p3 - p3)
	if ( ix(1) .eq. 0 ) ix(1) = 1D0

	D0 = 2*CYPLOG(-mudim/p1p2,1D0)*CYPLOG(-rmu*rm4/p2p3,1D0) - 
     &  CYPLOG(-rmu*rm4/p3,1D0)**2 - 2*CYPSp(1-x(1),-ix(1)) 
     &  - 2*pi**2/3
	D0 = (D0 + pi**2/4)/(p1p2*p2p3)
***]2-soft 3-coll (1-quasi)
	end
***********************************************************************
	subroutine CYPD0ir_2m0q_0p0a(fail,p1,p3,p1p2,p2p3,D0)
***********************************************************************
* Ref0. eqs. (4.47)
* Author: Le Duc Ninh
***********************************************************************
	implicit none
	CVAR D0
        CVAR p1, p3, p1p2, p2p3, m2, rm2, m4, rm4
        CVAR rmu, rla
	integer fail
#include "lt.h"

	D0 = 0
	fail = 0

	rmu = sqrt(mudim)
	if ( lambda .gt. 0 ) then
	 rla = sqrt(lambda)
	else
	 rla = lambda
	endif

***[2-soft 2-coll (2-quasi)
	m2 = p1
	m4 = p3
	rm2 = sqrt(m2)
	rm4 = sqrt(m4)

	 if( lambda .lt. 0 ) then
	  if( lambda .eq. -1 ) D0 = 2/(p1p2*(p2p3 - m4))*
     &	  CYPLOG(rm2*rm4/(m4-p2p3),1D0)
	  return
	 endif	

	if( lambda .eq. 0 ) then
	 D0 = 2*CYPLOG(-mudim/p1p2,1D0)*CYPLOG(rm2*rm4/(m4-p2p3),1D0)
     &	 /(p1p2*(p2p3 - m4))
	else
c	 print*,'D0ir check: new S2C2M1'
	 D0 = 2*CYPLOG(-lambda/p1p2,1D0)*CYPLOG(rm2*rm4/(m4-p2p3),1D0)
     &	 /(p1p2*(p2p3 - m4))
	endif
***]2-soft 2-coll (2-quasi)
	end
***********************************************************************
	subroutine CYPD0ir_1m0q_2p0a(fail,p1,p1p2,p2p3,D0)
***********************************************************************
* Ref0. eq. (4.75)
* Author: Le Duc Ninh
***********************************************************************
	implicit none
	CVAR D0
        CVAR p1, p1p2, p2p3, m2, rm2, m3, rm3, m4, rm4
        CVAR rmu, rla
	integer fail
#include "lt.h"

	D0 = 0
	fail = 0

	rmu = sqrt(mudim)
	if ( lambda .gt. 0 ) then
	 rla = sqrt(lambda)
	else
	 rla = lambda
	endif

***[1-soft 4-coll (4-quasi)
	m2 = p1
	rm2 = sqrt(m2)

	 if( lambda .lt. 0 ) then
	  if( lambda .eq. -1 ) D0 = 1/p1p2/p2p3*CYPLOG(-m2/p2p3,1D0)
	  return
	 endif

	if( lambda .eq. 0 ) then
	 D0 = 2*CYPLOG(-rmu*rm2/p1p2,1D0)
	else
	 D0 = 2*CYPLOG(-rla*rm2/p1p2,1D0)
	endif

	D0 = 1/p1p2/p2p3*(D0*CYPLOG(-m2/p2p3,1D0) 
     &   - pi**2/2D0)
***]1-soft 4-coll (4-quasi)
	end
***********************************************************************
	subroutine CYPD0ir_1m0q_0p0a(fail,p1,p2,p3,p4,p1p2,p2p3,m3,D0)
***********************************************************************
* Ref0. eq. (4.31)
* Author: Le Duc Ninh
***********************************************************************
	implicit none
	CVAR D0
        CVAR p1, p2, p3, p4, p1p2, p2p3, m2, rm2, m3, rm3, m4, rm4
        CVAR rmu, rla
* These are temporary variables. 
	CVAR x(3)
	RVAR ix(3)
	integer fail
#include "lt.h"

	D0 = 0
	fail = 0

	rmu = sqrt(mudim)
	if ( lambda .gt. 0 ) then
	 rla = sqrt(lambda)
	else
	 rla = lambda
	endif

***[1-soft 2-coll (2-quasi)
	m2 = p1
	m4 = p4
	rm2 = sqrt(m2)
	rm3 = sqrt(m3)
	rm4 = sqrt(m4)

c	print*,'Check D0 02, m2, m4, p2p3 =', m2, m4, p2p3

	 if( lambda .lt. 0 ) then
	  if( lambda .eq. -1 ) D0 = CYPLOG(-rm4*rm2/p2p3,1D0)
     &	  /p2p3/(p1p2-m3)
	  return
	 endif
	
	ix(1) = DBLE(p1p2 - p3)
	if ( ix(1) .eq. 0 ) ix(1) = 1D0

	ix(2) = DBLE(p1p2 - p2)
	if ( ix(2) .eq. 0 ) ix(2) = 1D0

	ix(3) = DBLE(p3 - p2)
	if ( ix(3) .eq. 0 ) ix(3) = 1D0
	
	D0 = -CYPLOG(-p2p3,-1D0)*CYPLOG(-rm4*rm2/p2p3,1D0)
     &	- (CYPLOG(-m4/p2p3,1D0)**2 + CYPLOG(-m2/p2p3,1D0)**2)/4 
     &  + CYPLOG(-m4/p2p3,1D0)*CYPLOG((m3-p3)/(m3-p1p2),ix(1)) 
     &  + CYPLOG(-m2/p2p3,1D0)*CYPLOG((m3-p2)/(m3-p1p2),ix(2)) 
     &  - CYPLOG((m3-p2)/(m3-p3),ix(3))**2/2

	x(1) = p2p3/(p2 - m3)
	ix(1) = DBLE(p2 - m3 - p2p3)
	if ( ix(1) .eq. 0 ) ix(1) = 1D0

	x(2) = m3/(m3 - p3)
	ix(2) = 1D0
	D0 = D0 + CYPLi2n(x, ix, 2) - pi**2/3 

	 if( lambda .eq. 0 ) then
          D0 = D0 + log(mudim)*CYPLOG(-rm4*rm2/p2p3,1D0)
	 else
c	  print*,'D0ir check: new S1C2M1'
          D0 = D0 + log(lambda)*CYPLOG(-rm4*rm2/p2p3,1D0)
	 endif
	 D0 = D0/p2p3/(p1p2-m3)
***]1-soft 2-coll (2-quasi)
	end
***********************************************************************
	subroutine CYPD0ir_2m0q_0p0b(fail,p1,p3,p4,p1p2,p2p3,D0)
***********************************************************************
* Ref0. eq. (4.39)
* Author: Le Duc Ninh
***********************************************************************
	implicit none
	CVAR D0
        CVAR p1, p3, p4, p1p2, p2p3, m2, rm2, m4, rm4
        CVAR rmu, rla
* These are temporary variables. 
	CVAR x(1)
	RVAR ix(1)
	integer fail
#include "lt.h"

	D0 = 0
	fail = 0

	rmu = sqrt(mudim)
	if ( lambda .gt. 0 ) then
	 rla = sqrt(lambda)
	else
	 rla = lambda
	endif

***[1-soft 2-coll (2-quasi)
	m2 = p1
	m4 = p4
	rm2 = sqrt(m2)
	rm4 = sqrt(m4)

* the soft epsilon pole can be translated to photon-mass log
	 if( lambda .lt. 0 ) then
	  if( lambda .eq. -1 ) D0 = 1/p1p2/(p2p3 - m4) 
     &    *CYPLOG(rm2*rm4/(m4 - p2p3),1D0)
	  return
	 endif

	 if( lambda .eq. 0 ) then
          D0 = log(mudim/m4)
	 else
          D0 = log(lambda/m4)
	 endif  

	ix(1) = DBLE(p1p2 - p3 + m4)
	if ( ix(1) .eq. 0 ) ix(1) = 1D0

	D0 = CYPLOG(rm2*rm4/(m4-p2p3),1D0)*
     &  (D0 + 2*CYPLOG((p3-m4)/p1p2,ix(1)))
     &  + CYPLOG(rm2*rm4/(m4 - p2p3),1D0)**2 + pi**2/6D0
	D0 = D0/p1p2/(p2p3 - m4)
***]1-soft 2-coll (2-quasi)
	end
***********************************************************************
	subroutine CYPD0ir_1m0q_0p0b(fail,p1,p2,p3,p4,p1p2,p2p3,m3,D0)
***********************************************************************
* Ref0. eq. (4.20)
* Author: Le Duc Ninh
***********************************************************************
	implicit none
	CVAR D0
        CVAR p1, p2, p3, p4, p1p2, p2p3, m2, rm2, m3, rm3, m4, rm4
        CVAR rmu, rla
* These are not temporary variables. 
	CVAR x34
	RVAR ix34, ix34s
* These are temporary variables. 
	CVAR x(2)
	RVAR ix(2)
	integer fail, l
#include "lt.h"

	D0 = 0
	fail = 0

	rmu = sqrt(mudim)
	if ( lambda .gt. 0 ) then
	 rla = sqrt(lambda)
	else
	 rla = lambda
	endif

***[1-soft 1-coll (1-quasi)
	m2 = p1
	m4 = p4
	rm2 = sqrt(m2)
	rm3 = sqrt(m3)
	rm4 = sqrt(m4)

	call CYPKam(p3, rm3, rm4, x34, ix34, ix34s)
	ix(1) = DBLE(p1p2 - p2)
	if ( ix(1) .eq. 0 ) ix(1) = 1D0

c	print*,'Check D0 01a, m2, m4, p2p3 =', m2, m4, p2p3

	 if( lambda .lt. 0 ) then
	  if( lambda .eq. -1 ) D0 = CYPLOG(rm4*rm2/(m4 - p2p3),1D0)
     &	  /(p2p3 - m4)/(p1p2-m3)
	  return
	 endif

	D0 = CYPLOG(rm4*rm2/(m4 - p2p3),1D0)*
     &	(2*CYPLOG((m3 - p2)/(m3 - p1p2),ix(1)) - log(m2)) + 
     &  CYPLOG(rm4*rm2/(m4 - p2p3),1D0)**2 - pi**2/6

c	print*,'Check D0 01b, m2, m4, x34 =', m2, m4, x34


        do l = 1, -1, -2
	x(1) = (m4 - p2p3)/(m3 - p2)
	ix(1) = DBLE(m4 - p2p3) - DBLE(m3 - p2) 
	if ( ix(1) .eq. 0 ) ix(1) = 1D0
	if ( abs(m3) .le. 0 ) then
	 x(2) = (m4 - p3)/m4
	 x(2) = -x(2)*(l-1)/2
	 ix(2) = -1D0
	else
	 x(2) = x34**l*rm3/rm4
	ix(2) = ix34*l 
	endif
	D0 = D0 + CYPLi2n(x, ix, 2)
	enddo     

	 if( lambda .eq. 0 ) then
          D0 = D0 + log(mudim)*CYPLOG(rm4*rm2/(m4 - p2p3),1D0)
	 else
c	  print*,'D0ir check: new S1C1M1'
          D0 = D0 + log(lambda)*CYPLOG(rm4*rm2/(m4 - p2p3),1D0)
	 endif  

	D0 = D0/(p1p2 - m3)/(p2p3 - m4)
***]1-soft 1-coll (1-quasi)
	end
***********************************************************************
	subroutine CYPD0ir_0m0q_4p0a(fail,p1p2,p2p3,m1,D0)
***********************************************************************
* Ref0. eq. (4.76)
* Author: Le Duc Ninh
***********************************************************************
	implicit none
	CVAR D0
        CVAR p1p2, p2p3, m1
        CVAR rmu, rla
	integer fail
#include "lt.h"

	D0 = 0
	fail = 0

	rmu = sqrt(mudim)
	if ( lambda .gt. 0 ) then
	 rla = sqrt(lambda)
	else
	 rla = lambda
	endif

***[0-soft 4-coll (4-quasi)
	 if( lambda .lt. 0 ) return

	D0 = 2/p1p2/p2p3*(CYPLOG(-m1/p1p2,1D0)*CYPLOG(-m1/p2p3,1D0) 
     &   - pi**2/2D0)
***]0-soft 4-coll (4-quasi)
	end
***********************************************************************
	subroutine CYPD0ir_1m0q_0p0c(fail,p1,p2,p3,p4,p1p2,p2p3,
     &  m3,m4,D0)
***********************************************************************
* Ref0. eq. (4.14)
* Author: Le Duc Ninh
***********************************************************************
	implicit none
	CVAR D0
        CVAR p1, p2, p3, p4, p1p2, p2p3, m2, rm2, m3, rm3, m4, rm4
        CVAR rmu, rla
* These are not temporary variables. 
	CVAR x34
	RVAR ix34, ix34s
* These are temporary variables. 
	CVAR x(3), cx(4)
	RVAR ix(3), icx(4)
	integer fail, l, k
#include "lt.h"

	D0 = 0
	fail = 0

	rmu = sqrt(mudim)
	if ( lambda .gt. 0 ) then
	 rla = sqrt(lambda)
	else
	 rla = lambda
	endif

***[0-soft 1-coll (1-quasi)
	m2 = p1
	rm2 = sqrt(m2)
	rm3 = sqrt(m3)
	rm4 = sqrt(m4)

c        print*,'D0ir.F: check 4.14'

	 if( lambda .lt. 0 ) return

	call CYPKam(p3, rm3, rm4, x34, ix34, ix34s)

* Ref0. (4.16)
	cx(1) = (m4 - p4)/(m3 - p1p2)
	icx(1) = DBLE(m4 - p4 - m3 + p1p2)
	if ( icx(1) .eq. 0 ) icx(1) = 1D0
	cx(2) = (m4 - p2p3)/(m3 - p2)
	icx(2) = DBLE(m4 - p2p3 - m3 + p2)
	if ( icx(2) .eq. 0 ) icx(2) = 1D0	

	if ( abs(m4) .gt. 0 ) then
	if ( abs(m3) .eq. 0 ) then
	 cx(3) = (m4 - p3)/m4
	 icx(3) = -1D0
	 cx(4) = 0
	 icx(4) = 1D0
	else
	 cx(3) = (rm3/rm4)/x34
	 icx(3) = -ix34
	 cx(4) = (rm3/rm4)*x34
	 icx(4) = ix34
	endif
	
	else
	 cx(3) = nan
	 icx(3) = 0
	 cx(4) = m3/(m3 - p3)
	 icx(4) = 1D0
	endif

* Ref0. (4.14)
	   if ( abs(m3) + abs(m4) .le. 0 ) then
* Ref0. (4.17)
c        print*,'D0ir.F: check 4.14 and 4.17'
	   x(1) = (m3 - p1p2)/(m3 - p2)
	   ix(1) = DBLE(-p1p2 + p2)
	   if ( ix(1) .eq. 0) ix(1) = 1D0
	   x(2) = (m4 - p4)/(m4 - p2p3)
	   ix(2) = DBLE(-p4 + p2p3)
	   if ( ix(2) .eq. 0) ix(2) = 1D0
	   x(3) = p3/p4
	   ix(3) = DBLE(p4 - p3)
	   if ( ix(3) .eq. 0) ix(3) = 1D0
	   D0 = - CYPLOG(x(1),ix(1))**2 + CYPLOG(x(2),ix(2))**2

	   x(1) = cx(2)
	   ix(1) = icx(2)
	   x(2) = 1/cx(1)
	   ix(2) = -icx(1)
	   D0 = D0 - 2*CYPLi2n(x,ix,2)

	   D0 = D0 + 
     &	   (CYPLOG(cx(1),icx(1))**2 - CYPLOG(cx(2),icx(2))**2)/2 + 
     &     (CYPLOG(cx(1),icx(1)) - CYPLOG(cx(2),icx(2)))*
     &	   (CYPLOG(-m2/p4,1D0) + CYPLOG(x(3),ix(3)))	   

	   elseif ( abs(m4) .le. 0 ) then
* cx(3) = nan if m4 = 0	   
* Ref0. (4.14b)
c        print*,'D0ir.F: check 4.14B'
	   D0 = (CYPLOG(cx(1),icx(1)) - CYPLOG(cx(2),icx(2)))*2*
     &	   CYPLOG(rm2*rm3/(m3 - p2),1D0)

	   x(1) = cx(2)
	   ix(1) = icx(2)
	   x(2) = 1/cx(1)
	   ix(2) = -icx(1)
	   D0 = D0 - 2*CYPLi2n(x,ix,2)

	   do k = 1, 2
	    do l = 3, 4
	     x(1) = 1/cx(k)
	     ix(1) = -icx(k)
	     if ( abs(cx(l) - nan) .le. prec) then
	     x(2) = 0
	     ix(2) = 0
	     else
	     x(2) = 1/cx(l)
	     ix(2) = -icx(l)
	     endif
	     D0 = D0 - (-1D0)**k*CYPLi2n(x,ix,2)
	    enddo
	   enddo 

	   else

* Ref0. (4.14a)
c        print*,'D0ir.F: check 4.14A'
	   D0 = (CYPLOG(cx(1),icx(1)) - CYPLOG(cx(2),icx(2)))*2*
     &	   CYPLOG(rm2*rm4/(m4 - p4),1D0)

	   x(1) = (m3 - p1p2)/(m3 - p2)
	   ix(1) = DBLE(-p1p2 + p2)
	   if ( ix(1) .eq. 0) ix(1) = 1D0
	   x(2) = (m4 - p4)/(m4 - p2p3)
	   ix(2) = DBLE(-p4 + p2p3)
	   if ( ix(2) .eq. 0) ix(2) = 1D0
	   D0 = D0 - CYPLOG(x(1),ix(1))**2 + CYPLOG(x(2),ix(2))**2

	   x(1) = cx(2)
	   ix(1) = icx(2)
	   x(2) = 1/cx(1)
	   ix(2) = -icx(1)
	   D0 = D0 - 2*CYPLi2n(x,ix,2)

	   do k = 1, 2
	    do l = 3, 4
	     x(1) = cx(k)
	     ix(1) = icx(k)
	     x(2) = cx(l)
	     ix(2) = icx(l)
	     D0 = D0 + (-1D0)**k*CYPLi2n(x,ix,2)
	    enddo
	   enddo 

c        print*,'D0ir.F: check 4.14 c', D0

	   endif
	   D0 = D0/((p1p2 - m3)*(p2p3 - m4) - (p2 - m3)*(p4 - m4))
c        print*,'D0ir.F: check 4.14', D0
***]0-soft 1-coll (1-quasi)
	end
***********************************************************************
      subroutine CYPD0per(mcase,per)
      implicit none
      integer mcase, per(10), 
     & permtable(10,0:15), i
      data permtable/
     & 1,2,3,4,5,6,7,8,9,10, ! 0, 4 masses zero (rot 1)
     & 1,2,3,4,5,6,7,8,9,10, ! 1, 3 masses zero (rot 1)
     & 1,6,3,5,4,2,7,8,10,9, ! 2, 3 masses zero (rot 2)
     & 1,2,3,4,5,6,7,8,9,10, ! 3, 2 masses zero (rot 1)
     & 5,3,6,1,4,2,7,9,10,8, ! 4, 3 masses zero (rot 4)
     & 5,2,6,4,1,3,7,9,8,10, ! 5, 2 masses zero (rot 3)
     & 4,6,2,5,1,3,7,10,8,9, ! 6, 2 masses zero (rot 5)
     & 1,2,3,4,5,6,7,8,9,10, ! 7, 1 masses zero (rot 1)
     & 2,3,4,1,6,5,8,9,10,7, ! 8, 3 masses zero (rot 10)
     & 2,5,4,6,1,3,8,9,7,10, ! 9, 2 masses zero (rot 9)
     & 6,4,5,2,1,3,8,10,7,9, ! 10, 2 masses zero (rot 11)
     & 1,5,3,6,2,4,8,7,9,10, ! 11, 1 masses zero (rot 7)
     & 3,4,1,2,5,6,9,10,7,8, ! 12, 2 masses zero (rot 17)
     & 5,1,6,3,2,4,9,7,8,10, ! 13, 1 masses zero (rot 13)
     & 4,1,2,3,6,5,10,7,8,9, ! 14, 1 masses zero (rot 19)
     & 1,2,3,4,5,6,7,8,9,10  ! 15, 0 masses zero (rot 1)
     &/

      do i = 1, 10
       per(i) = permtable(i, mcase)
      enddo
      end
