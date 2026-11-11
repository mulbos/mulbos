* func_basic.F
* basic functions

#include "defs.h"

	CVAR function CYPLOG(x, is)
************************************************************************
* Def. Log(x + i*is) where "is" is infinitesimal.
* Author: Le Duc Ninh (2011)
************************************************************************
	implicit none
	RVAR is
	CVAR x

#include "ff.h"

	CYPLOG = 0

***[singular cases
	if (abs(x) .le. 0) then
	 if ( runmode .ge. 3 ) then
	 print*,'Error in CYPLOG(): log(0) encountered, return 0'
	 call YPErrId(1)
	 endif
	 return
	endif
***]singular cases

	if (DIMAG(x) .eq. 0) then
	 if (DBLE(x) .gt. 0) then
	  CYPLOG = log(x)
	 else
	  if (abs(is) .eq. 0) then
	   if ( runmode .ge. 3 ) then
	   print*,'Error in CYPLOG(): x < 0 and is = 0, return 0'
	   call YPErrId(1)
	   endif
	   return
	  endif

	  CYPLOG = log(-x) + cipi*is/abs(is)
	 endif 	
	
	else
	CYPLOG = log(x)
	
	endif

	end
	
************************************************************************

	CVAR function CYPeta(ca,cb)
************************************************************************
* Def. eta = 2*i*pi*( thIm(-a)*thIm(-b)*thIm(a*b)
*				- thIm(a)*thIm(b)*thIm(-a*b) )
* with thIm(a) = theta(Im(a))
* WARNING: do not use this function for Real arguments with ieps. 
* See CYPetax() recommended.
* Author: Le Duc Ninh (2011)
************************************************************************
	implicit none
	CVAR ca,cb
	RVAR a,b,ab,rab

#include "ff.h"

	CYPeta = 0

	a = DIMAG(ca)
	b = DIMAG(cb)
	if ( a*b .lt. 0 ) return

	if ( a*b .eq. 0 ) then
	  if ( runmode .ge. 3 ) then
	  print*,'Error in CYPeta: Im(a, b) vanishes, return 0'
	  call YPErrId(1)
	  endif
	  return  
	endif

	rab = DBLE(ca)*DBLE(cb) - a*b
	ab = DBLE(ca)*b + a*DBLE(cb)

	if ( abs(ab) .le. 0 ) then
	  if ( runmode .ge. 3 ) then
	  print*,'Error in CYPeta: Im(ab) vanishes, return 0'
	  call YPErrId(1)
	  endif
	  return  
	endif

	if ( a .lt. 0 .and. b .lt. 0 .and. ab .gt. 0 ) then
	  CYPeta = 1
	elseif ( a .gt. 0 .and. b .gt. 0 .and. ab .lt. 0 ) then
	  CYPeta = -1
	else
	  CYPeta = 0
	endif
	
	CYPeta = 2*cipi*CYPeta

	end

************************************************************************
	CVAR function CYPetax(ca,cb,cc)
************************************************************************
* Def. etax = 2*i*pi*( thIm(-a)*thIm(-b)*thIm(c)
*				- thIm(a)*thIm(b)*thIm(-c) )
* with thIm(a) = theta(Im(a))
* Use this in any cases if possible.
* Author: Le Duc Ninh (2011)
************************************************************************

	implicit none
	CVAR ca,cb,cc
	RVAR a,b,ab

#include "ff.h"

	CYPetax = 0

	a = DIMAG(ca)
	b = DIMAG(cb)
	ab = DIMAG(cc)

	if ( a .gt. 0 .neqv. b .gt. 0 ) return

	if ( a*b*ab .eq. 0 ) then
	 if ( runmode .ge. 3 ) then
	 print*,'Error in CYPetax: Im(a, b, ab) vanishes, return 0'
	 call YPErrId(1)
	 endif
	 return  
	endif

	if ( a .lt. 0 .and. b .lt. 0 .and. ab .gt. 0 ) then
	 CYPetax = 1
	elseif ( a .gt. 0 .and. b .gt. 0 .and. ab .lt. 0 ) then
	 CYPetax = -1
	else
	 CYPetax = 0
	endif

	CYPetax = 2*cipi*CYPetax

	end

************************************************************************
	subroutine CYPRoot2(a, b, c, is, x, ix)
************************************************************************
* Ref0. Dao Thi Nhung and Le Duc Ninh arXiv: 0902.0325
* Find the 2 roots of the equation: a*x^2 + b*x + c + is*I = 0
* where "is" is infinitesimal, Ref0. Eq. (A.2).
* Input: a, b, c are complex; "is" is real.
* Output: x(2) and ix(2).
************************************************************************
	implicit none
	CVAR x(2)
	RVAR is, ix(2)
	CVAR a, b, c, d

#include "ff.h"

c	print*,'hello CYPRoot2'

	if ( abs(a) .le. 0 ) then
	 print*,'Error in CYPRoot2: a vanishes, STOP'
	 stop
	endif

	ix(1) = is

	if ( abs(c) .le. 0) then
	x(1) = 0
	x(2) = -b/a
	d = sqrt(b**2)
	if ( DBLE(d) .lt. 0 ) ix(1) = -ix(1) 

	else
	d = b**2 - 4*a*c
	d = sqrt(d)
	x(1) = (-b - d)/(2*a)
	x(2) = (-b + d)/(2*a)
	if ( DBLE(d) .lt. 0 ) ix(1) = -ix(1) 
	
	if ( abs(x(1)) .lt. abs(x(2)) ) then
	 x(1) = c/a/x(2)
	else
	 x(2) = c/a/x(1)
	endif 
	 
	endif

	ix(2) = -ix(1)

c	print*,'d =', d
c	print*,'x =', x

	end

************************************************************************
	CVAR function CYPSp(x, is)
************************************************************************
* Computes the dilogarithm (Li2 or Sp) for any (complex) x 	
* to a precision prec.
* Author: Le Duc Ninh (2011)
************************************************************************
	implicit none
	CVAR x, CYPSp1, CYPLOG
	RVAR is
	external CYPSp1, CYPLOG

#include "ff.h"

	CYPSp = 0

***[special cases
	if ( abs(x-1) .le. 0 ) then 
	 CYPSp = pi**2/6D0 
	 return
	endif  

	if ( abs(x-1) .le. prec ) then 
	 CYPSp = -(1-x) + pi**2/6D0 - CYPLOG(x, is)*CYPLOG(1-x, -is) 
	 return
	endif  
***]special cases

	if ( abs(x) .le. 1 .and. DBLE(x) .le. 0.5D0 ) then 
	 CYPSp = CYPSp1(x, is)

	elseif ( abs(x) .gt. 1 .and. DBLE(x) .le. 0.5D0 ) then
	 CYPSp = -CYPSp1(1/x, -is) - pi**2/6D0 - CYPLOG(-x, -is)**2/2
	
	elseif ( abs(x) .le. 1 .and. DBLE(x) .gt. 0.5D0 ) then
	 CYPSp = -CYPSp1(1-x, -is) + pi**2/6D0 
     &	  - CYPLOG(x, is)*CYPLOG(1-x, -is) 

	else
* here we have |x| > 1 AND Re(x) > 0.5	
	 if ( DBLE(1/x) .le. 0.5D0 ) then
	  CYPSp = -CYPSp1(1/x, -is) - pi**2/6D0 - CYPLOG(-x, -is)**2/2
	 else
	  CYPSp = -CYPSp1(1-x, -is) + pi**2/6D0 
     &	  - CYPLOG(x, is)*CYPLOG(1-x, -is) 
         endif
	 
	endif

	end

************************************************************************
	CVAR function CYPSp1(x, is)
************************************************************************
* Computes the dilogarithm (Li2 or Sp) for complex x with 
* |x| <= 1 and Re(x) <= 0.5 to a precision prec.
* Ref0. 't Hooft and Veltman, Nucl. Phys. B153 (1979) 365.
* Author: Le Duc Ninh (2011)
************************************************************************
	implicit none
	integer i, init, nn
	save init
	data init /0/
        parameter(nn = 60)
	CVAR x, z, CYPLOG, ctemp
	RVAR is, bern(nn), temp
	external CYPLOG

#include "ff.h"

	CYPSp1 = 0

c	print*,'From CYPSp1(): x, abs(x) =', x, abs(x)

	if ( abs(x) .gt. 1 .or. DBLE(x) .gt. 0.5D0 ) then 
	 print*,'Error in CYPSp11(): invalid input, STOP'
	 stop
	endif 

***[special cases
	if ( abs(x) .le. prec ) then 
	 CYPSp1 = x
	 return
	endif  

	if ( abs(x-0.5D0) .le. 0 ) then 
	 CYPSp1 = pi**2/12D0 - log(2D0)**2/2D0 
	 return
	endif
***]special cases

***[Bernouilli numbers
        if (init .eq. 0) then
        init = 1

        do i = 3, nn-1, 2
          bern(i) = 0
        enddo

        bern( 1) = -1D0/2D0
        bern( 2) =  1D0/6D0
        bern( 4) = -1D0/30D0
        bern( 6) =  1D0/42D0
        bern( 8) = -1D0/30D0
        bern(10) =  5D0/66D0
        bern(12) = -691D0/2730D0
        bern(14) =  7D0/6D0
        bern(16) = -3617D0/510D0
        bern(18) =  43867D0/798D0
        bern(20) = -174611D0/330D0
        bern(22) =  854513D0/138D0
        bern(24) = -236364091D0/2730D0
        bern(26) =  8553103D0/6D0
        bern(28) = -23749461029D0/870D0
        bern(30) =  8615841276005D0/14322D0
        bern(32) = -7709321041217D0/510D0
        bern(34) =  2577687858367D0/6D0
        bern(36) = -26315271553053477373D0/1919190D0
        bern(38) =  2929993913841559D0/6D0
        bern(40) = -261082718496449122051D0/13530D0

	bern(42) = 1520097643918070802691D0/1806D0
	bern(44) = -27833269579301024235023D0/690D0
	bern(46) = 596451111593912163277961D0/282D0
	bern(48) = -5609403368997817686249127547D0/46410D0
	bern(50) = 495057205241079648212477525D0/66D0
	bern(52) = -801165718135489957347924991853D0/1590D0
	bern(54) = 29149963634884862421418123812691D0/798D0
	bern(56) = -2479392929313226753685415739663229D0/870D0
	bern(58) = 84483613348880041862046775994036021D0/354D0
	bern(60) = -1215233140483755572040304994079820246041491D0/
     &   56786730D0
        endif
***]Bernouilli numbers

	z = -CYPLOG(1 - x, -is)
	ctemp = bern(1)*z**2/2
	CYPSp1 = z + ctemp 
        i = 0
	temp = 1

	do while ( abs(ctemp) .gt. prec )
	 i = i + 2

c	print*,'From CYPSp1(): i, z, ctemp =', i, z, ctemp

	 if ( i .gt. nn) then 
	  print*,'Error in CYPSp1(): 
     &	  Bernouilli series is too short, return at', nn
	  call YPErrId(1)
	  ctemp = 0
	 
	 else
	  temp = temp*i*(i+1)
	  ctemp = bern(i)*z**(i + 1)/temp	  
	 endif  

	CYPSp1 = CYPSp1 + ctemp
	enddo

	end

************************************************************************
	CVAR function CYPLi2n(x, ix, n)
************************************************************************
* Ref0. Denner and Dittmaier arXiv: 1005.2076
* Computes the function Li2(x1, x2, ..., xn) defined in 
* Ref0. Eqs. (2.13) and (2.14).
* Without ieps input for x1*x2*...*xn
* Use this function for complex x, or in cases where the ieps of 
* ALL the arguments have equal roles.
* See CYPLi2nx() recommended.
************************************************************************
	implicit none
	integer n, i, j
	CVAR x(*), ctemp, sum, ca, cb, cc
	RVAR ix(*), itemp, temp
	CVAR CYPSp, CYPLOG, CYPetax
	external CYPSp, CYPLOG, CYPetax

#include "ff.h"

	CYPLi2n = 0
	ctemp = 1D0
	do i = 1, n
	 ctemp = ctemp*x(i)
	enddo

c	print*,'Lin2n(): ctemp, n=',ctemp, n

	if ( abs(ctemp) .le. 0 ) then
	 CYPLi2n = pi**2/6
	 return
	endif 
	
	itemp = 0
	do i = 1, n
	 temp = 1D0
	 do j = 1, n
	  if ( j .ne. i ) temp = temp*DBLE(x(j))
	 enddo

	 itemp = itemp + temp*ix(i)
	enddo
	if ( itemp .eq. 0 ) itemp = 1D0
	
***[special case
	if ( n .eq. 2 ) then
	 ca = x(1)
	 if ( DIMAG(ca) .eq. 0) ca = DCMPLX(0D0,ix(1))
	 cb = x(2)
	 if ( DIMAG(cb) .eq. 0) cb = DCMPLX(0D0,ix(2))
	 cc = ctemp
	 if ( DIMAG(cc) .eq. 0) cc = DCMPLX(0D0,itemp)
	 
	 CYPLi2n = CYPetax(ca,cb,cc)
	 if ( abs(CYPLi2n) .eq. 0 ) then 
	 CYPLi2n = CYPSp(1 - ctemp, -itemp)	 
	 else

	 if ( abs(ctemp - 1) .le. 0 ) then
	 if ( runmode .ge. 3 ) then
	 print*,'Error in CYPLi2n(): x1*x2*...*xn = 1, return 0'
	 call YPErrId(1)
	 endif
	 return
	 endif 

	 CYPLi2n = CYPSp(1 - ctemp, -itemp) + 
     &	  CYPLOG(1 - ctemp, -itemp)*CYPLi2n
	 endif
	 return
	endif
***]special case

	sum = 0
	do i = 1, n
	 sum = sum + CYPLOG(x(i), ix(i)) 
	enddo

	CYPLi2n = CYPSp(1 - ctemp, -itemp)
	ca = CYPLOG(ctemp, itemp) - sum
	if ( abs(ca) .gt. prec ) then

	if ( abs(ctemp - 1) .le. 0 ) then
	 if ( runmode .ge. 3 ) then
	 print*,'Error in CYPLi2n(): x1*x2*...*xn = 1, return 0'
	 call YPErrId(1)
	 endif
	 return
	endif

	cb = CYPLOG(1 - ctemp, -itemp)
	if ( abs(ctemp) - 1 .gt. 0 ) cb = cb -
     &   CYPLOG(-ctemp, -itemp) + (CYPLOG(ctemp, itemp) + sum)/2
	CYPLi2n = CYPLi2n + ca*cb 
	endif

	end

************************************************************************
	CVAR function CYPLi2nx(x, ix, is, n)
************************************************************************
* Ref0. Denner and Dittmaier arXiv: 1005.2076
* Computes the function Li2(x1, x2, ..., xn) defined in 
* Ref0. Eqs. (2.13) and (2.14).
* Same as CYPLi2n() BUT with ieps input (is) for x1*x2*...*xn
* Use this in any cases if possible.
************************************************************************
	implicit none
	integer n, i, j
	CVAR x(*), ctemp, sum, ca, cb, cc
	RVAR ix(*), is, temp
	CVAR CYPSp, CYPLOG, CYPetax
	external CYPSp, CYPLOG, CYPetax

#include "ff.h"

	CYPLi2nx = 0
	ctemp = 1D0
	do i = 1, n
	 ctemp = ctemp*x(i)
	enddo

c	print*,'Lin2n(): ctemp, n=',ctemp, n

	if ( abs(ctemp) .le. 0 ) then
	 CYPLi2nx = pi**2/6
	 return
	endif 
	
***[special case
	if ( n .eq. 2 ) then
	 ca = x(1)
	 if ( DIMAG(ca) .eq. 0) ca = DCMPLX(0D0,ix(1))
	 cb = x(2)
	 if ( DIMAG(cb) .eq. 0) cb = DCMPLX(0D0,ix(2))
	 cc = ctemp
	 if ( DIMAG(cc) .eq. 0) cc = DCMPLX(0D0,is)
	 
	 CYPLi2nx = CYPetax(ca,cb,cc)
	 if ( abs(CYPLi2nx) .eq. 0 ) then 
	 CYPLi2nx = CYPSp(1 - ctemp, -is)	 
	 else

	 if ( abs(ctemp - 1) .le. 0 ) then
	 if ( runmode .ge. 3 ) then
	 print*,'Error in CYPLi2nx(): x1*x2*...*xn = 1, return 0'
	 call YPErrId(1)
	 endif
	 return
	 endif 

	 CYPLi2nx = CYPSp(1 - ctemp, -is) + 
     &	  CYPLOG(1 - ctemp, -is)*CYPLi2nx
	 endif
	 return
	endif
***]special case

	sum = 0
	do i = 1, n
	 sum = sum + CYPLOG(x(i), ix(i)) 
	enddo

	CYPLi2nx = CYPSp(1 - ctemp, -is)
	ca = CYPLOG(ctemp, is) - sum
	if ( abs(ca) .gt. prec ) then

	 if ( abs(ctemp - 1) .le. 0 ) then
	 if ( runmode .ge. 3 ) then
	 print*,'Error in CYPLi2nx(): x1*x2*...*xn = 1, return 0'
	 call YPErrId(1)
	 endif
	 return
	 endif

	cb = CYPLOG(1 - ctemp, -is)
	if ( abs(ctemp) - 1 .gt. 0 ) cb = cb -
     &   CYPLOG(-ctemp, -is) + (CYPLOG(ctemp, is) + sum)/2
	CYPLi2nx = CYPLi2nx + ca*cb 
	endif

	end

************************************************************************      
      CVAR function CYPLambda(x,y,z)
      implicit none
      CVAR x,y,z
      CYPLambda=x**2+y**2+z**2-2*x*y-2*x*z-2*z*y
      end

************************************************************************
      subroutine CYPKam(p2, rm2, rm3, xs, ixs, ixs2)
************************************************************************
* Ref0. Denner and Dittmaier arXiv: 1005.2076
* Computes xij defined by Ref0. Eq. (2.7).
* Input: p2, rm2, rm3; with p2 = p2 + i*eps, eps > 0. 
* Output: xij = xs; ixs and ixs2 are the ieps part of xs and xs**2. 
* Author: LE Duc Ninh (2011).
************************************************************************
      implicit none
      CVAR p2, rm2, rm3, xs
      RVAR ixs, ixs2

#include "ff.h"

	if ( abs(p2 - (rm2 - rm3)**2) .le. 0 ) then
	 xs = 1
         ixs = 1 
	else
	 xs = sqrt(1 - 4*rm2*rm3/(p2 - (rm2 - rm3)**2))
         ixs = DBLE(rm2*rm3)*DBLE(xs)
         if ( ixs .eq. 0 ) ixs = 1
	 xs = (xs - 1)/(xs + 1) 
	endif 

	ixs2 = ixs*DBLE(xs)
        if ( ixs2 .eq. 0 ) ixs2 = ixs	

      end

************************************************************************
	CVAR function YPLi2(x)
************************************************************************
* Computes the dilogarithm (Li2 or Sp) for real x.
* to be used in softphoton.F in FormCalc. 
* This unecessary routine is not used in this library.
* Author: LE Duc Ninh (2011).
************************************************************************
	implicit none
	CVAR CYPSp
	RVAR x
	external CYPSp

	YPLi2 = x
	YPLi2 = CYPSp(YPLi2, 0D0)
        end
************************************************************************
      subroutine CYPSetM0(p,n)
************************************************************************
* To set small non-vanishing masses to zero when it is collinear safe.
* Limit: up to 6-point
* Author: LE Duc Ninh (2011).
************************************************************************
      implicit none
      integer n
      CVAR p(n)
      RVAR apprx
* input: choose 10 for up to 10-point, 
* increase this for higher points
      RVAR x(10)

#include "ff.h"

***[2-point integrals
	if ( n .eq. 3 ) then
c       apprx = precd
       apprx = prec

      if ( abs(p(2)) .gt. 0 ) then
       x(1) = max(abs(p(1)),abs(p(3)))
       if ( abs(p(2)) .le. apprx*x(1) ) p(2) = 0
      endif

      if ( abs(p(3)) .gt. 0 ) then
       x(1) = max(abs(p(1)),abs(p(2)))
       if ( abs(p(3)) .le. apprx*x(1) ) p(3) = 0
      endif

	return
	endif
***]2-point integrals

***[3-point integrals
	if ( n .eq. 6 ) then
c       apprx = precd
       apprx = stab0

      if ( abs(p(4)) .gt. 0 ) then
       x(1) = max(abs(p(1)),abs(p(5)))
       x(2) = max(abs(p(3)),abs(p(6)))
       x(2) = min(x(1),x(2))

       if ( abs(p(4)) .le. apprx*x(2) ) p(4) = 0
      endif

      if ( abs(p(5)) .gt. 0 ) then
       x(1) = max(abs(p(1)),abs(p(4)))
       x(2) = max(abs(p(2)),abs(p(6)))
       x(2) = min(x(1),x(2))

       if ( abs(p(5)) .le. apprx*x(2) ) p(5) = 0
      endif

      if ( abs(p(6)) .gt. 0 ) then
       x(1) = max(abs(p(2)),abs(p(5)))
       x(2) = max(abs(p(3)),abs(p(4)))
       x(2) = min(x(1),x(2))

       if ( abs(p(6)) .le. apprx*x(2) ) p(6) = 0
      endif

	return
	endif
***]3-point integrals

***[4-point integrals
      if ( n .eq. 10 ) then
       apprx = stab0

      if ( abs(p(7)) .gt. 0 ) then
       x(1) = max(abs(p(1)),abs(p(8)))
       x(2) = max(abs(p(4)),abs(p(10)))
       x(3) = max(abs(p(5)),abs(p(9)))
       x(3) = min(x(1),x(2),x(3))

       if ( abs(p(7)) .le. apprx*x(3) ) p(7) = 0
      endif

      if ( abs(p(8)) .gt. 0 ) then
       x(1) = max(abs(p(1)),abs(p(7)))
       x(2) = max(abs(p(2)),abs(p(9)))
       x(3) = max(abs(p(6)),abs(p(10)))
       x(3) = min(x(1),x(2),x(3))

       if ( abs(p(8)) .le. apprx*x(3) ) p(8) = 0
      endif

      if ( abs(p(9)) .gt. 0 ) then
       x(1) = max(abs(p(2)),abs(p(8)))
       x(2) = max(abs(p(3)),abs(p(10)))
       x(3) = max(abs(p(5)),abs(p(7)))
       x(3) = min(x(1),x(2),x(3))

       if ( abs(p(9)) .le. apprx*x(3) ) p(9) = 0
      endif

      if ( abs(p(10)) .gt. 0 ) then
       x(1) = max(abs(p(3)),abs(p(9)))
       x(2) = max(abs(p(4)),abs(p(7)))
       x(3) = max(abs(p(6)),abs(p(8)))
       x(3) = min(x(1),x(2),x(3))

       if ( abs(p(10)) .le. apprx*x(3) ) p(10) = 0
      endif
	return
      endif
***]4-point integrals

***[5-point integrals
	if ( n .eq. 15 ) then 
       apprx = stab0

      if ( abs(p(11)) .gt. 0 ) then
       x(1) = max(abs(p(1)),abs(p(12)))
       x(2) = max(abs(p(5)),abs(p(15)))
       x(3) = max(abs(p(6)),abs(p(13)))
       x(4) = max(abs(p(9)),abs(p(14)))
       x(4) = min(x(1),x(2),x(3),x(4))

       if ( abs(p(11)) .le. apprx*x(4) ) p(11) = 0
      endif

      if ( abs(p(12)) .gt. 0 ) then
       x(1) = max(abs(p(1)),abs(p(11)))
       x(2) = max(abs(p(2)),abs(p(13)))
       x(3) = max(abs(p(10)),abs(p(15)))
       x(4) = max(abs(p(7)),abs(p(14)))
       x(4) = min(x(1),x(2),x(3),x(4))

       if ( abs(p(12)) .le. apprx*x(4) ) p(12) = 0
      endif

      if ( abs(p(13)) .gt. 0 ) then
       x(1) = max(abs(p(2)),abs(p(12)))
       x(2) = max(abs(p(3)),abs(p(14)))
       x(3) = max(abs(p(6)),abs(p(11)))
       x(4) = max(abs(p(8)),abs(p(15)))
       x(4) = min(x(1),x(2),x(3),x(4))

       if ( abs(p(13)) .le. apprx*x(4) ) p(13) = 0
      endif

      if ( abs(p(14)) .gt. 0 ) then
       x(1) = max(abs(p(3)),abs(p(13)))
       x(2) = max(abs(p(4)),abs(p(15)))
       x(3) = max(abs(p(7)),abs(p(12)))
       x(4) = max(abs(p(9)),abs(p(11)))
       x(4) = min(x(1),x(2),x(3),x(4))

       if ( abs(p(14)) .le. apprx*x(4) ) p(14) = 0
      endif

      if ( abs(p(15)) .gt. 0 ) then
       x(1) = max(abs(p(4)),abs(p(14)))
       x(2) = max(abs(p(5)),abs(p(11)))
       x(3) = max(abs(p(8)),abs(p(13)))
       x(4) = max(abs(p(10)),abs(p(12)))
       x(4) = min(x(1),x(2),x(3),x(4))

       if ( abs(p(15)) .le. apprx*x(4) ) p(15) = 0
      endif

	return
	endif
***]5-point integrals

***[6-point integrals
	if ( n .eq. 21 ) then
       apprx = stab0

      if ( abs(p(16)) .gt. 0 ) then
       x(1) = max(abs(p(1)),abs(p(17)))
       x(2) = max(abs(p(7)),abs(p(18)))
       x(3) = max(abs(p(13)),abs(p(19)))
       x(4) = max(abs(p(11)),abs(p(20)))
       x(5) = max(abs(p(6)),abs(p(21)))
       x(5) = min(x(1),x(2),x(3),x(4),x(5))

       if ( abs(p(16)) .le. apprx*x(5) ) p(16) = 0
      endif

      if ( abs(p(17)) .gt. 0 ) then
       x(1) = max(abs(p(1)),abs(p(16)))
       x(2) = max(abs(p(2)),abs(p(18)))
       x(3) = max(abs(p(8)),abs(p(19)))
       x(4) = max(abs(p(14)),abs(p(20)))
       x(5) = max(abs(p(12)),abs(p(21)))
       x(5) = min(x(1),x(2),x(3),x(4),x(5))

       if ( abs(p(17)) .le. apprx*x(5) ) p(17) = 0
      endif

      if ( abs(p(18)) .gt. 0 ) then
       x(1) = max(abs(p(7)),abs(p(16)))
       x(2) = max(abs(p(2)),abs(p(17)))
       x(3) = max(abs(p(3)),abs(p(19)))
       x(4) = max(abs(p(9)),abs(p(20)))
       x(5) = max(abs(p(15)),abs(p(21)))
       x(5) = min(x(1),x(2),x(3),x(4),x(5))

       if ( abs(p(18)) .le. apprx*x(5) ) p(18) = 0
      endif

      if ( abs(p(19)) .gt. 0 ) then
       x(1) = max(abs(p(13)),abs(p(16)))
       x(2) = max(abs(p(8)),abs(p(17)))
       x(3) = max(abs(p(3)),abs(p(18)))
       x(4) = max(abs(p(4)),abs(p(20)))
       x(5) = max(abs(p(10)),abs(p(21)))
       x(5) = min(x(1),x(2),x(3),x(4),x(5))

       if ( abs(p(19)) .le. apprx*x(5) ) p(19) = 0
      endif

      if ( abs(p(20)) .gt. 0 ) then
       x(1) = max(abs(p(11)),abs(p(16)))
       x(2) = max(abs(p(14)),abs(p(17)))
       x(3) = max(abs(p(9)),abs(p(18)))
       x(4) = max(abs(p(4)),abs(p(19)))
       x(5) = max(abs(p(5)),abs(p(21)))
       x(5) = min(x(1),x(2),x(3),x(4),x(5))

       if ( abs(p(20)) .le. apprx*x(5) ) p(20) = 0
      endif

      if ( abs(p(21)) .gt. 0 ) then
       x(1) = max(abs(p(6)),abs(p(16)))
       x(2) = max(abs(p(12)),abs(p(17)))
       x(3) = max(abs(p(15)),abs(p(18)))
       x(4) = max(abs(p(10)),abs(p(19)))
       x(5) = max(abs(p(5)),abs(p(20)))
       x(5) = min(x(1),x(2),x(3),x(4),x(5))

       if ( abs(p(21)) .le. apprx*x(5) ) p(21) = 0
      endif
	return
	endif
***]6-point integrals

***[7-point integrals
	if ( n .gt. 21 ) then
	 print*,'SetM0() in func_basic.F: 7-point not done'
	 stop
	 return
	endif
***]7-point integrals
      end
***********************************************************************
	integer function YPirk(ir,k)
***********************************************************************
* Index function for tensor reduction
* Ref0: Denner and Dittmaier hep-ph/0509141
* Ref0. Eq. (2.6)
* Author: LE Duc Ninh (2011)
***********************************************************************
	implicit none
	integer ir, k
	if (k .gt. ir ) then
	 YPirk = ir
	else if ( k .lt. ir ) then
	 YPirk = ir - 1
	else
c	 print*,'Error in YPirk(): invalid input ir = k' 
c	 stop
* can be anything here
	 YPirk = 0
	endif
	end
***********************************************************************
	integer function YPDltB(i,j)
***********************************************************************
* Index function for tensor reduction
* Ref0: Denner and Dittmaier hep-ph/0509141
* Ref0. after Eq. (2.9)
* Author: LE Duc Ninh (2011)
***********************************************************************
	implicit none
	integer i, j
	if (i .eq. j ) then
	 YPDltB = 0
	else
	 YPDltB = 1
	endif
	end
************************************************************************
      subroutine YPErrId(id)
************************************************************************
* Error Identification routine.
* Input: id to know the error level (not used by now).
* The error signals: 2 for B0, 3 for C0, 5 for D0.
* Author: LE Duc Ninh (2011)
************************************************************************
      implicit none
      integer id
#include "ff.h"
* WARNING: iers has to be smaller than the integer limit.
      if ( iers .lt. 5**9 ) then
	if ( MOD(iers,2) .eq. 0 ) iers = iers*2
	if ( MOD(iers,3) .eq. 0 ) iers = iers*3
	if ( MOD(iers,5) .eq. 0 ) iers = iers*5
      endif
      end

