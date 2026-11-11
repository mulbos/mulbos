* Ccoeffs.F

#include "defs.h"

*#define SCALAR_COLL_QUAD
*#define DEBUG_Ccoeffs

************************************************************************
	subroutine CYPCcoeff(para, C, ldpara)
************************************************************************
* Method: If abs(detG) is large enough we use double precision,
* otherwise the coefficients C(ccIJK...) are calculated 
* in quadruple precision.   
* Possible approximation here to avoid numerical cancellation for 
* m << p. See CYPSetM0(). This is a very good approximation.
* Author: Le Duc Ninh (2011)
************************************************************************
	implicit none
	integer ldpara

	CVAR para(ldpara,Pcc)
	CVAR paraC(Pcc)
	CVAR C(Ncc)

	CVAR32 paraCQP(Pcc)
	CVAR32 CQP(Ncc)

	CVAR G(2,2), detG
	integer i
	logical dump

#include "lt.h"

***[DEBUG
	dump = ibits(debugkey, DebugC, 1) .ne. 0
	if( dump ) call CYPDumpPara(3, para, ldpara, "CYPCcoeff")
***]DEBUG

	do i = 1, Ncc
	 C(i) = 0
	enddo

* Ninh check
c	return

	
	do i = 1, Pcc
	 paraC(i) = para(1, i)
	enddo	
	
***[set small masses to 0 when it's safe
	call CYPSetM0(paraC, Pcc)
***]set small masses to 0 when it's safe	

	G(1,1) = 2*paraC(1)
	G(2,2) = 2*paraC(3)
	G(1,2) = paraC(1) + paraC(3) - paraC(2)
	G(2,1) = G(1,2)

	detG = G(1,1)*G(2,2)-G(1,2)*G(1,2)
        pmax = max(abs(paraC(1)),abs(paraC(2)),abs(paraC(3)))

        if ( abs(detG) .gt. stab3*4*pmax**2 .or. abs(detG) .eq. 0 ) then
***[using default precision
	call CYPCscalar(paraC, C(cc0))	
	call CYPCtensor(paraC, C)
***]using default precision
	else

***[using quadruple precision
	do i = 1, Pcc
	 paraCQP(i) = para(1, i)
	enddo	

	call CscalarCQP(paraCQP, CQP(cc0))	
	call CtensorCQP(paraCQP, CQP)

	do i = 1, Ncc
	 C(i) = CQP(i)
	enddo	
***]using quadruple precision

	endif

***[DEBUG
	if( dump ) call CYPDumpCoeff(3, C)
***]DEBUG

	end
************************************************************************

	subroutine CYPCscalar(para, fx)
************************************************************************
* For Symmetry see Stuart, Comp. Phys. Comm. 48 (1988) 367
* Ref0. Dao Thi Nhung and Le Duc Ninh arXiv: 0902.0325
* Ref1. Denner and Dittmaier arXiv: 1005.2076
* Author: Le Duc Ninh
************************************************************************
	implicit none
	CVAR para(Pcc), mm, fx, fy
	CVAR p1, p2, p1p2, m1, m2, m3, ctemp
        RVAR raa,rbb,rcc,rgg,rhh,signf,signj
        CVAR cdd,cee,cff,cjj
	integer fail

#ifdef SCALAR_COLL_QUAD
        RVAR16 raa1,rbb1,rcc1,rgg1,rhh1,signf1,signj1
        CVAR32 cdd1,cee1,cff1,cjj1,C0QP
#endif

#include "lt.h"

***[versions
	integer key
	logical check
	key = ibits(versionkey, KeyC0, 2)
C NInh checki
c        key = 1

	check = .FALSE.
***]versions

	fx = 0
	fail = 1
	iers = 3

cccccccccccccccccccccccccccccccccccccccc
#ifdef NINH_CHECK_TEMP
C NInh checki (use OneLoop and LoopTools res here)
c	print*,'C0.F hello, key =', key
        if ( abs(para(4)+para(5)+para(6) - 3*173.10D0**2)/(3*173.10D0**2) .lt. 1D-13 ) then
c	print*,'C0.F hello, key =', abs(para(4)) - 173.10D0**2, 
c     & (DBLE(para(1)+para(2)+para(3)) + 347309.684406793094240128993988037D0), 
c     & (DBLE(para(1)+para(2)+para(3)) + 422.916107644501607865095138549805D0)
c        stop
        IF ( abs(DBLE(para(1)+para(2)+para(3)) + 347309.684406793094240D0)/347309.68D0 .lt. 1D-13 
     & .and. abs(para(1)+para(2)) .lt. 1D-10 ) then
         fx = -9.761565580424066d-006
         return
        ENDIF

        IF ( abs(DBLE(para(1)+para(2)+para(3)) + 422.916107644501607865095D0)/422.9D0 .lt. 1D-13 
     & .and. abs(para(2)+para(3)) .lt. 1D-10 ) then
         fx = -1.666731770297139d-005
         return
        ENDIF
        endif
#endif
cccccccccccccccccccccccccccccccccccccccc        

	p1 = para(1)
	p2 = para(2)
	p1p2 = para(3)
	m1 = para(4)
	m2 = para(5)
	m3 = para(6)

        pmmax = max(abs(p1),abs(p2),abs(p1p2),abs(m1),abs(m2),abs(m3))
        pmax = max(abs(p1),abs(p2),abs(p1p2))
c        pmin = min(abs(p1),abs(p2),abs(p1p2))

***[soft
	if( abs(m1) .le. 0 .and.
     &      (abs(p1 - m2) + abs(p1p2 - m3)) .le. 0 ) then

	 if( abs(m2) .le. abs(m3) ) then
	  call CYPC0ir(fail, p1, p2, p1p2, m1, m2, m3, fx)
         else
	  call CYPC0ir(fail, p1p2, p2, p1, m1, m3, m2, fx)
	 endif 
	 goto 301
	endif 

	if( abs(m2) .le. 0 .and.
     &      (abs(p1 - m1) + abs(p2 - m3)) .le. 0 ) then

	 if( abs(m1) .le. abs(m3) ) then
	  call CYPC0ir(fail, p1, p1p2, p2, m2, m1, m3, fx)
         else
	  call CYPC0ir(fail, p2, p1p2, p1, m2, m3, m1, fx)
	 endif
	 goto 301
	endif 

	if( abs(m3) .le. 0 .and.
     &      (abs(p2 - m2) + abs(p1p2 - m1)) .le. 0 ) then

	 if( abs(m1) .le. abs(m2) ) then
	  call CYPC0ir(fail, p1p2, p1, p2, m3, m1, m2, fx)
	 else
	  call CYPC0ir(fail, p2, p1, p1p2, m3, m2, m1, fx)	 
	 endif
	 goto 301
	endif 
***]soft

***[coll
* m1 = 0
	if( abs(m1) .le. 0 .and. abs(p1 - m2) .le. 0 
     &	.and. abs(m2) .le. stab0*abs(pmmax) ) then
         call CYPC0ir(fail, p1, p2, p1p2, m1, m2, m3, fx)
	 goto 301
	endif 

	if( abs(m1) .le. 0 .and. abs(p1p2 - m3) .le. 0 
     &	.and. abs(m3) .le. stab0*abs(pmmax) ) then
         call CYPC0ir(fail, p1p2, p2, p1, m1, m3, m2, fx)
	 goto 301
	endif 

* m3 = 0
	if( abs(m3) .le. 0 .and. abs(p1p2 - m1) .le. 0 
     &	.and. abs(m1) .le. stab0*abs(pmmax) ) then
         call CYPC0ir(fail, p1p2, p1, p2, m3, m1, m2, fx)
	 goto 301
	endif 

	if( abs(m3) .le. 0 .and. abs(p2 - m2) .le. 0 
     &	.and. abs(m2) .le. stab0*abs(pmmax) ) then
         call CYPC0ir(fail, p2, p1, p1p2, m3, m2, m1, fx)
	 goto 301
	endif 

* m2 = 0
	if( abs(m2) .le. 0 .and. abs(p2 - m3) .le. 0 
     &	.and. abs(m3) .le. stab0*abs(pmmax) ) then
	 call CYPC0ir(fail, p2, p1p2, p1, m2, m3, m1, fx)
	 goto 301
	endif 

	if( abs(m2) .le. 0 .and. abs(p1 - m1) .le. 0 
     &	.and. abs(m1) .le. stab0*abs(pmmax) ) then
	 call CYPC0ir(fail, p1, p1p2, p2, m2, m1, m3, fx)
	 goto 301
	endif 

* p1 = 0
* we order (p1, p2, p1p2) where abs(p2) > abs(p1p2)
	if( abs(p1) .le. 0 .and.
     &   abs(m1 - m2) .le. 0 .and. abs(m2) .le. stab0*abs(pmmax) ) then

	 if ( abs(p2) .ge. abs(p1p2) ) then
	  call CYPC0ir(fail, p1, p2, p1p2, m1, m2, m3, fx)
	 else
	  call CYPC0ir(fail, p1, p1p2, p2, m2, m1, m3, fx)
	 endif

	 goto 301
	endif  

* p1p2 = 0
	if( abs(p1p2) .le. 0 .and.
     &   abs(m3 - m1) .le. 0 .and. abs(m1) .le. stab0*abs(pmmax) ) then

	 if ( abs(p1) .ge. abs(p2) ) then
	  call CYPC0ir(fail, p1p2, p1, p2, m3, m1, m2, fx)
	 else
	  call CYPC0ir(fail, p1p2, p2, p1, m1, m3, m2, fx)
	 endif

	 goto 301
	endif  

* p2 = 0
	if( abs(p2) .le. 0 .and.
     &   abs(m2 - m3) .le. 0 .and. abs(m3) .le. stab0*abs(pmmax) ) then

	 if ( abs(p1p2) .ge. abs(p1) ) then
	  call CYPC0ir(fail, p2, p1p2, p1, m2, m3, m1, fx)
	 else
	  call CYPC0ir(fail, p2, p1, p1p2, m3, m2, m1, fx)
	 endif

	 goto 301
	endif  
***]coll

 301    continue
	if ( fail .eq. 0 ) goto 302

***[finite cases
	if ( lambda .lt. 0 ) goto 302

***[finite with zero masses
	if ( pmmax .gt. 0 ) then
* the infinity mass
	mm = pmmax/prec

	if ( abs(m1) + abs(m2) + abs(m3) .le. 0 ) then
* pinch Mbox1, set Pbox1 = 0
	call CYPD0M0(fail,czero,p2,p1p2,p1,p2,p1,mm,m1,m2,m3,fx)

	elseif ( abs(m2) + abs(m3) .le. 0 ) then
* pinch Mbox4, set Pbox4 = 0
	call CYPD0M0(fail,p1,p2,p1p2,czero,p1p2,p1,m1,m2,m3,mm,fx)

	elseif ( abs(m1) + abs(m3) .le. 0 ) then
* pinch Mbox4, set Pbox4 = 0
	call CYPD0M0(fail,p2,p1p2,czero,p1,p1,p1p2,m2,m3,m1,mm,fx)

	elseif ( abs(m1) + abs(m2) .le. 0 ) then
* pinch Mbox4, set Pbox4 = 0
	call CYPD0M0(fail,p1p2,p1,p1,p1p2,p2,czero,m3,m1,m2,mm,fx)

	elseif ( abs(m3) .le. 0 ) then
* pinch Mbox4, set Pbox4 = 0
	call CYPD0M0(fail,p1,p2,p1p2,czero,p1p2,p1,m1,m2,m3,mm,fx)

	elseif ( abs(m1) .le. 0 ) then
* pinch Mbox4, set Pbox4 = 0
	call CYPD0M0(fail,p2,p1p2,czero,p1,p1,p1p2,m2,m3,m1,mm,fx)

	elseif ( abs(m2) .le. 0 ) then
* pinch Mbox4, set Pbox4 = 0
	call CYPD0M0(fail,p1p2,p1,p1,p1p2,p2,czero,m3,m1,m2,mm,fx)

	endif

	if ( fail .eq. 0 ) then
	fx = -mm*fx
	goto 302
	endif

	endif
***]finite with zero masses

***[finite without zero masses
	if( key .ne. 1 ) then
***[finite without zero masses, method of Ref1.
	if ( pmmax .gt. 0 ) then
* the infinity mass
	mm = pmmax/prec

* pinch Mbox3, set p1p2 = 0
	call CYPD0M0(fail,p1,p1,p1p2,p1p2,czero,p2,m1,m2,mm,m3,fx)

c	print*,'C0 check DD: fail, fx =', fail, fx
cccccccccccccccccccccccccccccccccccccccc
#ifdef NINH_CHECK_TEMP
***[equal masses
c        if ( (abs(m1 - m2) + abs(m1 - m3))/pmmax .lt. prec*1D2 ) then
c        endif
***]equal masses

C NInh checki (use OneLoop and LoopTools res here)
c	print*,'C0.F hello, key =', key
        if ( abs(para(4)+para(5)+para(6) - 3*173.10D0**2)/(3*173.10D0**2) .lt. 1D-13 ) then
c	print*,'C0.F hello, key =', abs(para(4)) - 173.10D0**2, 
c     & (DBLE(para(1)+para(2)+para(3)) + 347309.684406793094240128993988037D0), 
c     & (DBLE(para(1)+para(2)+para(3)) + 422.916107644501607865095138549805D0)
c        stop
        IF ( abs(DBLE(para(1)+para(2)+para(3)) + 347309.684406793094240D0)/347309.68D0 .lt. 1D-13 
     & .and. abs(para(1)+para(2)) .lt. 1D-10 ) then
c         fx = -9.761565580424066d-006
	print*,'C0 check DD: fail, fx =', fail, fx
         stop
         return
        ENDIF

        IF ( abs(DBLE(para(1)+para(2)+para(3)) + 422.916107644501607865095D0)/422.9D0 .lt. 1D-13 
     & .and. abs(para(2)+para(3)) .lt. 1D-10 ) then
c         fx = -1.666731770297139d-005
	print*,'C0 check DD: fail, fx =', fail, fx
         stop
         return
        ENDIF
        endif
#endif
cccccccccccccccccccccccccccccccccccccccc        

	if ( fail .eq. 0 ) then
	 fy = -mm*fx
	 if( key .eq. 0 ) then
	  fx = fy
	  goto 302
	 else
	  check = .TRUE.
	 endif
	endif

	endif
***]finite without zero masses, method of Ref1.
	endif

	if( key .ne. 0 .or. fail .ne. 0 ) then
***[finite without zero masses, method of Ref0.
***[rotate to a good position
	 mmax = max(abs(m1), abs(m2), abs(m3))
	 
c	print*,'C0 fin hello 0'
c	print*, p1, p2, p1p2, m1, m2, m3


	 if ( abs(p1p2)*abs(p2) .gt. 0 .and. abs(p1) .le. 0 ) then 
	  ctemp = p1
	  p1 = p1p2
	  p1p2 = p2
	  p2 = ctemp
	  
	  ctemp = m1
	  m1 = m3
	  m3 = m2
	  m2 = ctemp

	 else
	  if ( abs(m2) .eq. mmax ) then
	   ctemp = p2
	   p2 = p1p2
	   p1p2 = ctemp
	   ctemp = m1
	   m1 = m2
	   m2 = ctemp
	  elseif ( abs(m3) .eq. mmax ) then
	   ctemp = p1
	   p1 = p1p2
	   p1p2 = p2
	   p2 = ctemp
	  
	   ctemp = m1
	   m1 = m3
	   m3 = m2
	   m2 = ctemp
	  endif

	 endif
***]rotate to a good position

c	print*,'C0 fin hello 1'
c	print*, p1, p2, p1p2, m1, m2, m3

* possible option:
***[pseudo-coll singular cases left
#ifdef SCALAR_COLL_QUAD
* use quadruple prec
	if ( fail .eq. 2 ) then
	  raa1 = DBLE(p1p2)
	  rbb1 = DBLE(p2)
	  rcc1 = DBLE(p1 - p2 - p1p2)
	  rgg1 = 0
	  rhh1 = 0
	  cdd1 = m3 - m1 - p1p2
	  cee1 = m2 - m3 + p1p2 - p1
	  cff1 = m1
	  cjj1 = -1
	  signf1 = -1
	  signj1 = signf1

	  call fTCQP(raa1,rbb1,rcc1,rgg1,rhh1,cdd1,cee1,cff1,cjj1,
     >     signf1,signj1,C0QP)
	  fx = C0QP
	  goto 302
	endif
#endif
***]pseudo-coll singular cases left

	  raa = DBLE(p1p2)
	  rbb = DBLE(p2)
	  rcc = DBLE(p1 - p2 - p1p2)
	  rgg = 0
	  rhh = 0
	  cdd = m3 - m1 - p1p2
	  cee = m2 - m3 + p1p2 - p1
	  cff = m1
	  cjj = -1
	  signf = -1
	  signj = signf

	  call CYPfT(raa,rbb,rcc,rgg,rhh,cdd,cee,cff,cjj,
     >     signf,signj,fx)

	 if( key .gt. 1 .and. check .and. 
     &        abs(fx - fy) .gt. maxdev*abs(fx) ) then
cc	  if ( mod(iers,9) .eq. 0 ) then 
	   print *, "Discrepancy in C0:"
	   print*,'C0DD  =',fy
	   print*,'C0abc =',fx
	   print *, "    p1   =", para(1)
	   print *, "    p2   =", para(2)
	   print *, "    p1p2 =", para(3)
	   print *, "    m1 =", para(4)
	   print *, "    m2 =", para(5)
	   print *, "    m3 =", para(6)

c	   fx = fy
cc	  endif
	 endif

***]finite without zero masses, method of Ref0.
	endif
***]finite without zero masses

c	print*,'From CYPCscalar(): C0 =', fx
***]finite cases

 302    continue

***[check for errors
	if ( mod(iers,9) .eq. 0 ) then 
	print *, "CYPCscalar in Ccoeffs.F:"
	print *, "    C0 =", fx
	print *, "    p1   =", para(1)
	print *, "    p2   =", para(2)
	print *, "    p1p2 =", para(3)
	print *, "    m1 =", para(4)
	print *, "    m2 =", para(5)
	print *, "    m3 =", para(6)
	endif
***]check for errors
	iers = 1

	end

************************************************************************
	subroutine CYPC0ir(fail, p1, p2, p1p2, m1, m2, m3, C0)
************************************************************************
* Calculate the infrared (soft and collinear) scalar 3-point integrals.
* Ref0. Dittmaier hep-ph/0308246
* Ref0a. Denner and Dittmaier arXiv: 1005.2076
* Ref1. Beenakker and Denner Nucl. Phys. B338, 349
* Ref2. Ellis and Zanderighi arXiv: 0712.1851v3 
* Output: fail = 0 (OK), fail = 2 (nothing done).
* Important: the correct positions are m1 = 0 and m1 = m2 = very small.
* Convention: 
* same as Ref0a (similar to Ref0) for mass regularization (MR), 
* same as Ref2 for dimensional regularization (DR). 
* Translation rules from Ref0a to Ref2 for DR: 
* Ref2 = z0 + z1/eps + z2/eps^2,
* Ref0a = d0 + d1/eps + d2/eps^2,
* then: z2 = d2, z1 = d1, z0 = d0 + d2*pi^2/6. 
* Author: Le Duc Ninh
************************************************************************
	implicit none
	CVAR p1, p2, p1p2, m1, m2, m3, xs
	CVAR rm1, rm2, rm3
	CVAR C0
	RVAR ixs, ixs2

	integer fail

#include "lt.h"

	C0 = 0
	fail = 0

	if( abs(m1) .gt. 0 .and. abs(m1 - m2) .gt. 0 ) then
	 print*,'Error in CYPC0ir(): invalid input, check rotation'
	 stop
	endif
	
	rm1 = sqrt(m1)
	rm2 = sqrt(m2)
	rm3 = sqrt(m3)

c	print*,'C0 IR hello'
c	print*, p1, p2, p1p2, m1, m2, m3


***[pure soft and soft + coll
* m2 < m3 and cannot be zero here

	if( abs(m1) .le. 0 .and.
     &      (abs(p1 - m2) + abs(p1p2 - m3)) .le. 0 ) then


c	print*,'C0 soft hello'

	if ( abs(m2) .le. 0 .and. abs(m3) .le. 0 ) then
* see Ref0. (B.16)
	if ( lambda .lt. 0 ) then
	 if ( lambda .eq. -2 ) C0 = 1/p2
	 if ( lambda .eq. -1 ) C0 = CYPLOG(-mudim/p2,1D0)/p2 
	 return
	endif	

	C0 = CYPLOG(-mudim/p2,1D0)**2/2
	C0 = C0/p2
	return
	endif

	if ( abs(m2) .le. 0 ) then
* see Ref0. (B.12)
	if ( lambda .lt. 0 ) then
	 if ( lambda .eq. -2 ) C0 = 1/(p2 - m3)/2
	 if ( lambda .eq. -1 ) C0 = (CYPLOG(mudim/(m3-p2),1D0) - 
     &	  log(mudim/m3)/2)/(p2-m3) 
	 return
	endif	

	C0 = (CYPLOG(mudim/(m3-p2),1D0)**2 - 
     &	  log(mudim/m3)**2/2)/2 - CYPSp(p2/(p2-m3),-1D0)
     &   + pi**2/12
	C0 = C0/(p2-m3)
	return
	endif

***[both m2 and m3 are small
* see Ref1. (C.6)
	if( abs(m2) .lt. abs(p2)*stab0 .and.
     &      abs(m3) .lt. abs(p2)*stab0 ) then
	 if( lambda .lt. 0 ) then
	  if( lambda .eq. -1 ) C0 = CYPLOG(-rm2*rm3/p2,1D0)/p2
          return
	 endif 
	
	 C0 = (-CYPLOG(-rm2*rm3/p2,1D0)*log(rm2*rm3) - pi**2/6 +
     &    CYPLOG(-m2/p2,1D0)*CYPLOG(-m3/p2,1D0)/2)/p2
	 if( lambda .eq. 0 ) then
          C0 = C0 + log(mudim)*CYPLOG(-rm2*rm3/p2,1D0)/p2

	 else
          C0 = C0 + log(lambda)*CYPLOG(-rm2*rm3/p2,1D0)/p2
	 endif  
	return
	endif
***]both m2 and m3 are small

***[m2 is small
* see Ref0. (B.8)
* m2 < m3 and cannot be zero here
	if( abs(m2) .lt. abs(p2 - m3)*stab0 ) then
	 if( lambda .lt. 0 ) then
	 if( lambda .eq. -1 ) C0 = -CYPLOG((m3-p2)/rm2/rm3,-1D0)/(p2-m3)
         return
	 endif 
	
	 C0 = (CYPSp(p2/m3,1D0) + CYPLOG((m3-p2)/rm2/rm3,-1D0)*
     &	  CYPLOG(rm2*(m3-p2)/rm3,-1D0))/(p2-m3)

	 if( lambda .eq. 0 ) then
          C0 = C0 - log(mudim)*CYPLOG((m3-p2)/rm2/rm3,-1D0)/(p2-m3)

	 else
          C0 = C0 - log(lambda)*CYPLOG((m3-p2)/rm2/rm3,-1D0)/(p2-m3)
	 endif  
	return
	endif
***]m2 is small

        call CYPKam(p2, rm2, rm3, xs, ixs, ixs2)
	
	if ( abs(xs - 1) .le. 0 ) then 
* Ref1. and Ref2.
*** dimensional regularization 
	 if( lambda .lt. 0 ) then
	  if( lambda .eq. -1 ) C0 = 1/rm2/rm3/2
          return
	 endif 

	 if ( abs(rm3 - rm2) .le. 0 ) then
* take the limit rm2/rm3 -> 1
         C0 = -log(rm2*rm3)
	 else
         C0 = -log(rm2*rm3) - 2 - (rm3 + rm2)/(rm3 - rm2)*log(rm2/rm3)
	 endif
         C0 = C0/2/rm2/rm3
*** dimensional regularization 
	 if( lambda .eq. 0 ) then
          C0 = C0 + log(mudim)/2/rm2/rm3

*** mass regularization
	 else
          C0 = C0 + log(lambda)/2/rm2/rm3
	 endif  

c	print*,'C0 IR special 1, C0 =', C0	
	else
*** dimensional regularization 
	 if( lambda .lt. 0 ) then
	  if( lambda .eq. -1 ) 
     &	  C0 = -CYPLOG(xs,ixs)*xs/rm2/rm3/(1 - xs**2)
          return
	 endif  

*** the finite part
	 C0 = log(rm2*rm3)*CYPLOG(xs,ixs) - 
     &	 CYPLOG(xs,ixs)**2/2 + 2*CYPLOG(xs,ixs)*CYPLOG(1 - xs**2,-ixs2)
     &   + log(rm2/rm3)**2/2 - pi**2/6 + CYPSp(xs**2,ixs2) + 
     &   CYPSp(1 - xs*rm2/rm3,-ixs*DBLE(rm2/rm3)) +
     &   CYPSp(1 - xs*rm3/rm2,-ixs*DBLE(rm3/rm2))
         C0 = C0*xs/rm2/rm3/(1 - xs**2)

*** dimensional regularization 
	 if( lambda .eq. 0 ) then
          C0 = C0 - log(mudim)*CYPLOG(xs,ixs)*xs/rm2/rm3/(1 - xs**2)

*** mass regularization
	 else
          C0 = C0 - log(lambda)*CYPLOG(xs,ixs)*xs/rm2/rm3/(1 - xs**2)
	 endif  

        endif 

	return
	endif 
***]pure soft and soft + coll

c	print*,'C0 coll hello'
c	print*, p1, p2, p1p2, m1, m2, m3

***[pure coll
* Ref. Dittmaier hep-ph/0308246v1
	if( abs(m1) .le. 0 .and. 
     &   abs(p1) .le. 0 .and. abs(m2) .le. 0 ) then
* see (B.4)
	 if ( lambda .lt. 0 ) then
	 if ( lambda .eq. -1 ) 
     &	  C0 = (CYPLOG(m3 - p1p2,-1D0) 
     &    - CYPLOG(m3 - p2,-1D0))/(p2 - p1p2)
	 return
	 endif

	if ( abs(m3) .le. 0 ) then
	 C0 = (CYPLOG(p2,1D0) - CYPLOG(p1p2,1D0))*
     &	 (CYPLOG(-p2/mudim,-1D0) + CYPLOG(-p1p2/mudim,-1D0))/2
	else
	 C0 = log(mudim/m3)*(CYPLOG(m3 - p1p2,-1D0) 
     &    - CYPLOG(m3 - p2,-1D0))
     &   + CYPLOG((m3-p2)/m3,-1D0)**2 - CYPLOG((m3-p1p2)/m3,-1D0)**2
     &   + CYPSp(p2/m3,1D0) - CYPSp(p1p2/m3,1D0)
	endif
	 C0 = C0/(p2 - p1p2)
	 return
	endif

	if( abs(p1) .le. 0 .and. abs(m1 - m2) .le. 0 .and. 
     &   abs(m2 - m3) .le. 0 .and. abs(p1p2) .le. 0 .and. 
     &   abs(m1) .lt. stab0*abs(p2) ) then
* see (B.15)
	C0 = CYPLOG(-m1/p2,1D0)**2/p2/2
	return
	endif

	if( abs(p1) .le. 0 .and. 
     &   abs(m1 - m2) .le. 0 .and. 
     &   abs(m2) .lt. stab0*abs(m3 - p2) .and. 
     &   abs(m3 - p1p2) .le. 0 ) then
* see (B.11)
	 C0 = (CYPLOG((m3 - p2)/(rm3*rm1),-1D0)**2 + pi**2/12 
     &    + CYPSp(p2/m3,1D0))/(p2 - m3)
	 return
	endif

	if( abs(p1) .le. 0 .and. 
     &   abs(m1 - m2) .le. 0 .and. 
     &   abs(m2) .lt. stab0*abs(m3 - p2) .and. 
     &   abs(m2) .lt. stab0*abs(m3 - p1p2) ) then
* see (B.3)
	 if ( abs(m3) .le. 0 ) then
	 C0 = (CYPLOG(p2,1D0) - CYPLOG(p1p2,1D0))*
     &	 (CYPLOG(-p2/m1,-1D0) + CYPLOG(-p1p2/m1,-1D0))/2
	 else
	 C0 = -CYPLOG((m3 - p1p2)/(rm1*rm3),-1D0)**2 
     &   + CYPLOG((m3 - p2)/(rm1*rm3),-1D0)**2
     &   + CYPSp(p2/m3,1D0) - CYPSp(p1p2/m3,1D0)
	 endif
	 C0 = C0/(p2 - p1p2)
	 return
	endif

	if( abs(m1) .le. 0 .and. 
     &   abs(p1 - m2) .le. 0 .and. 
     &   abs(m2) .lt. stab0*abs(m3-p1p2) ) then
	 if ( abs(p2 - m3) .le. 0 ) then
* see (B.10)
	  C0 = (CYPLOG((m3-p1p2)/(rm3*rm2),-1D0)**2 + 5*pi**2/12 
     &     + CYPSp(p1p2/m3,1D0))/(p1p2 - m3)	  
	  return
	 endif
	endif

	if( abs(m1) .le. 0 .and. 
     &   abs(p1 - m2) .le. 0 .and. 
     &   abs(m2) .lt. stab0*abs(m3-p2) .and. 
     &   abs(m2) .lt. stab0*abs(m3-p1p2) ) then
* see (B.2)
	 ixs = DBLE(p2 - p1p2)
	 if ( ixs .eq. 0 ) ixs = 1D0

	if ( abs(m3) .le. 0 ) then
	 C0 = CYPLOG(p2/p1p2,-ixs)*(CYPLOG(-p2/m2,-1D0) + 
     &	 CYPLOG(-p1p2/m2,-1D0))/2 - 2*CYPSp((p2-p1p2)/(-p1p2),ixs)

	else
	 C0 = CYPLOG((m3-p2)/m2,-1D0)*CYPLOG((m3-p2)/m3,-1D0)
     &    - CYPLOG((m3-p1p2)/m2,-1D0)*CYPLOG((m3-p1p2)/m3,-1D0)
     &    + CYPSp(p2/m3,1D0)	 
     &    - CYPSp(p1p2/m3,1D0)
     &    - 2*CYPSp((p2-p1p2)/(m3-p1p2),ixs) 
	endif
	 C0 = C0/(p2 - p1p2)
	 return
	endif

	if( abs(m1) .le. 0 .and. 
     &   abs(p1 - m2) .le. 0 .and. abs(m2) .lt. stab0*abs(p1p2) ) then
	 if ( abs(m3) .le. 0 .and. abs(p2 - m2) .le. 0 ) then
* see (B.13)
	  C0 = (CYPLOG(-m2/p1p2,1D0)**2/2 + 2*pi**2/3)/p1p2
	  return
	 endif 

	 if ( abs(m3 - m2) .le. 0 .and. abs(p2) .le. 0 ) then
* see (B.14)
	  C0 = (CYPLOG(-m2/p1p2,1D0)**2/2 + pi**2/3)/p1p2
	  return
	 endif 

	endif 
***]pure coll

	fail = 2
#ifdef DEBUG_Ccoeffs
	print*,'From CYPC0ir(): pseudo-coll singular cases remain'
	print*,'From CYPC0ir(): C0 =', C0
	print *, "    stab0   =", stab0
	print *, "    p1   =", p1
	print *, "    p2   =", p2
	print *, "    p1p2 =", p1p2
	print *, "    m1 =", m1
	print *, "    m2 =", m2
	print *, "    m3 =", m3
#endif
 
	end

************************************************************************
	subroutine CYPCtensor(para, C)
************************************************************************
* Method: we use Passarino-Veltman reduction method for DetGram .ne. 0.  
* This subroutine is taken and adapted from LoopTools of Thomas Hahn.  
************************************************************************
	implicit none
	integer i
	CVAR para(Pcc)
	CVAR p1, p2, p1p2, m1, m2, m3, f1, f2
	CVAR G(2,2), detG, alpha(0:2)
	CVAR C(Ncc)
	CVAR bsum, b1sum, b00sum, b11sum, in(2)
	integer B12, B23, B13, Di(0:2)

!	CVAR aa

	integer XYPBget
	external XYPBget

#include "lt.h"

	integer perm(2)
#define IN(i) in(perm(i))
#define SOLVE_SETUP CYPDecomp(2, G,2, perm)
#define SOLVE(b) CYPSolve(2, G,2, b)

	do i = 2, Ncc
	 C(i) = 0
	enddo

	p1 = para(1)
	p2 = para(2)
	p1p2 = para(3)
	m1 = para(4)
	m2 = para(5)
	m3 = para(6)

	f1 = m2 - m1 - p1
	f2 = m3 - m1 - p1p2

	G(1,1) = 2*p1
	G(2,2) = 2*p1p2
	G(1,2) = p1 + p1p2 - p2
	G(2,1) = G(1,2)

	detG = G(1,1)*G(2,2)-G(1,2)*G(1,2)
!        pmax = Max(abs(p1),abs(p2),abs(p1p2))

       if ( abs(detG) .eq. 0D0 ) then
        return
       endif

***[using Passarino-Veltman method

	B12 = XYPBget(p1, m1, m2)
	B23 = XYPBget(p2, m2, m3)
	B13 = XYPBget(p1p2, m1, m3)
	Di(0)=B23
	Di(1)=B13
	Di(2)=B12

	bsum = BvalYP(bb0,B23) + BvalYP(bb1,B23)
	b1sum = BvalYP(bb1,B23) + BvalYP(bb11,B23)
	b00sum = BvalYP(bb00,B23) + BvalYP(bb001,B23)
	b11sum = BvalYP(bb11,B23) + BvalYP(bb111,B23)

	call SOLVE_SETUP

	IN(1) = f1*C(cc0) - BvalYP(bb0,B23) + BvalYP(bb0,B13)
	IN(2) = f2*C(cc0) - BvalYP(bb0,B23) + BvalYP(bb0,B12)
	call SOLVE(in)
	C(cc1) = in(1)
	C(cc2) = in(2)

	C(cc00) = .5D0*(m1*C(cc0) -
     &    .5D0*(f1*C(cc1) + f2*C(cc2) - BvalYP(bb0,B23) - 1*xrat))

	IN(1) = f1*C(cc1) + bsum - 2*C(cc00)
	IN(2) = f2*C(cc1) + bsum + BvalYP(bb1,B12)
	call SOLVE(in)
	C(cc11) = in(1)
	C(cc12) = in(2)

	IN(1) = f1*C(cc2) - BvalYP(bb1,B23) + BvalYP(bb1,B13)
	IN(2) = f2*C(cc2) - BvalYP(bb1,B23) - 2*C(cc00)
	call SOLVE(in)
	C(cc12) = .5D0*(C(cc12) + in(1))
	C(cc22) = in(2)

!!![rank3
       if ( rankmax3 .ge. 3 ) then

	C(cc001) = 1/3D0*(m1*C(cc1) -
     &    .5D0*(f1*C(cc11) + f2*C(cc12) + bsum + 1/3D0*xrat))
	C(cc002) = 1/3D0*(m1*C(cc2) -
     &    .5D0*(f1*C(cc12) + f2*C(cc22) - BvalYP(bb1,B23) 
     &    + 1/3D0*xrat))

	bsum = bsum + b1sum
	IN(1) = f1*C(cc11) - bsum - 4*C(cc001)
	IN(2) = f2*C(cc11) - bsum + BvalYP(bb11,B12)
	call SOLVE(in)
	C(cc111) = in(1)
	C(cc112) = in(2)

	IN(1) = f1*C(cc22) - BvalYP(bb11,B23) + BvalYP(bb11,B13)
	IN(2) = f2*C(cc22) - BvalYP(bb11,B23) - 4*C(cc002)
	call SOLVE(in)
	C(cc122) = in(1)
	C(cc222) = in(2)

        endif
!!!]rank3

!!![rank4
       if ( rankmax3 .ge. 4 ) then
	C(cc0000) = 1/4D0*(m1*C(cc00) -
     &    .5D0*(f1*C(cc001) + f2*C(cc002) - BvalYP(bb00,B23) -
     &          (m1 + m2 + m3 - .25D0*(p1 + p2 + p1p2))/6D0*xrat))

***{Not good, contain 1/detG

c	IN(1) = f1*C(cc001) + b00sum - 2*C(cc0000)
c	IN(2) = f2*C(cc001) + b00sum + BvalYP(bb001,B12)
c	call SOLVE(in)
c	C(cc0011) = in(1)
c	C(cc0012) = in(2)

c	IN(1) = f1*C(cc002) - BvalYP(bb001,B23) + BvalYP(bb001,B13)
c	IN(2) = f2*C(cc002) - BvalYP(bb001,B23) - 2*C(cc0000)
c	call SOLVE(in)
c	C(cc0012) = .5D0*(C(cc0012) + in(1))
c	C(cc0022) = in(2)

***}Not good, contain 1/detG      

***{Ninh change, checked
       C(cc0011)=-BvalYP(bb0,Di(0))-2*BvalYP(bb1,Di(0))
     >  -BvalYP(bb11,Di(0))-2*C(cc11)*m1+C(cc111)*f1+C(cc112)*f2
       C(cc0011)=(-C(cc0011)+1/6D0*xrat)/8D0
       
       C(cc0012)=BvalYP(bb1,Di(0))+BvalYP(bb11,Di(0))
     >  -2*C(cc12)*m1+C(cc112)*f1 + C(cc122)*f2       
       C(cc0012)=(-C(cc0012)+1/12D0*xrat)/8D0
       
       C(cc0022)=-BvalYP(bb11,Di(0)) - 2*C(cc22)*m1
     >  + C(cc122)*f1 + C(cc222)*f2
       C(cc0022)=(-C(cc0022)+1/6D0*xrat)/8D0     
***}Ninh change, checked

	bsum = bsum + b1sum + b11sum
	IN(1) = f1*C(cc111) + bsum - 6*C(cc0011)
	IN(2) = f2*C(cc111) + bsum + BvalYP(bb111,B12)
	call SOLVE(in)
	C(cc1111) = in(1)
	C(cc1112) = in(2)

	IN(1) = f1*C(cc222) - BvalYP(bb111,B23) + BvalYP(bb111,B13)
	IN(2) = f2*C(cc222) - BvalYP(bb111,B23) - 6*C(cc0022)
	call SOLVE(in)
	C(cc1222) = in(1)
	C(cc2222) = in(2)

	IN(1) = f1*C(cc122) + b11sum - 2*C(cc0022)
	IN(2) = f2*C(cc122) + b11sum - 4*C(cc0012)
	call SOLVE(in)
	C(cc1122) = in(1)
	C(cc1222) = .5D0*(C(cc1222) + in(2))

       endif
!!!]rank4
	
***[PV rank5
	if ( rankmax3 .ge. 5 ) then
***{ Ninh change
* Ref: Denner and Dittmaier hep-ph/0509141 and hep-ph/0212259
      
      C(cc00001)=BvalYP(bb00,Di(0))+BvalYP(bb001,Di(0))-2*C(cc001)*m1
     >  + C(cc0011)*f1 + C(cc0012)*f2	
      C(cc00001)=-C(cc00001)
     > +(2*p2-5*m1+(p1-5*m2)*2+(p1p2-5*m3))/120D0*xrat
      C(cc00001)=C(cc00001)/10D0
      
      C(cc00002)=-BvalYP(bb001,Di(0))-2*C(cc002)*m1+C(cc0012)*f1
     >  + C(cc0022)*f2
      C(cc00002)=-C(cc00002)
     > +(2*p2-5*m1+(p1-5*m2)+(p1p2-5*m3)*2)/120D0*xrat
      C(cc00002)=C(cc00002)/10D0
      
***{cross-check      
c      IN(1)=-(BvalYP(bb0000,Di(0))-BvalYP(bb0000,Di(1))-C(cc0000)*f1)
c      IN(2)=-(BvalYP(bb0000,Di(0))-BvalYP(bb0000,Di(2))-C(cc0000)*f2)
c	call SOLVE(in)
c	C(cc00001) = in(1)
c	C(cc00002) = in(2)

c      print*,'Check cc00001 =',C(cc00001)
c      print*,'Check cc00002 =',C(cc00002)      	         
***}cross-check      

      C(cc00111)=BvalYP(bb0,Di(0))+3*BvalYP(bb1,Di(0)) 
     > + 3*BvalYP(bb11,Di(0))
     > +BvalYP(bb111,Di(0))-2*C(cc111)*m1+ C(cc1111)*f1+ C(cc1112)*f2
      C(cc00111)=(-C(cc00111)-1/10D0*xrat)/10D0
      
      C(cc00112)=-BvalYP(bb1,Di(0))- 2*BvalYP(bb11,Di(0))
     > - BvalYP(bb111,Di(0))-2*C(cc112)*m1+ C(cc1112)*f1
     > + C(cc1122)*f2
      C(cc00112)=(-C(cc00112)-1/30D0*xrat)/10D0
      
      C(cc00122)=BvalYP(bb11,Di(0))+ BvalYP(bb111,Di(0))
     > -2*C(cc122)*m1+ C(cc1122)*f1+ C(cc1222)*f2
      C(cc00122)=(-C(cc00122)-1/30D0*xrat)/10D0
      
      C(cc00222)=-BvalYP(bb111,Di(0)) - 2*C(cc222)*m1
     > + C(cc1222)*f1 + C(cc2222)*f2      
      C(cc00222)=(-C(cc00222)-1/10D0*xrat)/10D0
      
***{cross-check
c       IN(1)=BvalYP(bb00,Di(0))+ 2*BvalYP(bb001,Di(0))
c     > + BvalYP(bb0011,Di(0))- C(cc0011)*f1 + 4*C(cc00001)
c       IN(1)=-IN(1)
c       IN(2)=BvalYP(bb00,Di(0))+ 2*BvalYP(bb001,Di(0))
c     > + BvalYP(bb0011,Di(0))- BvalYP(bb0011,Di(2))- C(cc0011)*f2
c       IN(2)=-IN(2)  
        
c	call SOLVE(in)
c	C(cc00111) = in(1)
c	C(cc00112) = in(2)
	
c       IN(1)=-(BvalYP(bb0011,Di(0))-BvalYP(bb0011,Di(1))-C(cc0022)*f1)
c       IN(2)=-(BvalYP(bb0011,Di(0))-C(cc0022)*f2+4*C(cc00002))
c       	call SOLVE(in)
c	C(cc00122) = in(1)
c	C(cc00222) = in(2)

c      print*,'Check cc00111 =',C(cc00111)
c      print*,'Check cc00112 =',C(cc00112)      	               
c      print*,'Check cc00122 =',C(cc00122)
c      print*,'Check cc00222 =',C(cc00222)      	                     

***}cross-check            
       IN(1)=BvalYP(bb0,Di(0))+ 4*BvalYP(bb1,Di(0))
     > + 6*BvalYP(bb11,Di(0)) + 4*BvalYP(bb111,Di(0))
     > + BvalYP(bb1111,Di(0)) - C(cc1111)*f1
     > + 8*C(cc00111)
       IN(2)=BvalYP(bb0,Di(0)) + 4*BvalYP(bb1,Di(0))
     > + 6*BvalYP(bb11,Di(0)) + 4*BvalYP(bb111,Di(0))
     > + BvalYP(bb1111,Di(0)) - BvalYP(bb1111,Di(2))
     > - C(cc1111)*f2
       IN(1)=-IN(1)
       IN(2)=-IN(2)
        call SOLVE(in)
	C(cc11111) = in(1)
	C(cc11112) = in(2)
	
       IN(1)= BvalYP(bb11,Di(0)) + 2*BvalYP(bb111,Di(0))
     > + BvalYP(bb1111,Di(0))- C(cc1122)*f1 + 4*C(cc00122)
       IN(2)=BvalYP(bb11,Di(0))+ 2*BvalYP(bb111,Di(0))
     > + BvalYP(bb1111,Di(0))- C(cc1122)*f2+ 4*C(cc00112)
       IN(1)=-IN(1)
       IN(2)=-IN(2)
        call SOLVE(in)
	C(cc11122) = in(1)
	C(cc11222) = in(2)     
     
       IN(1)=BvalYP(bb1111,Di(0))- BvalYP(bb1111,Di(1))- C(cc2222)*f1
       IN(2)=BvalYP(bb1111,Di(0))- C(cc2222)*f2 + 8*C(cc00222)
       IN(1)=-IN(1)
       IN(2)=-IN(2)
        call SOLVE(in)
	C(cc12222) = in(1)
	C(cc22222) = in(2)
***} Ninh change
	endif
***]PV rank5

***]using Passarino-Veltman method
	
	end
	
	
	
