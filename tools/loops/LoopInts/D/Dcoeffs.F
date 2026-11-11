* Dcoeffs.F
#include "defs.h"

* use also Denner-Dittmaier (DD) method (working)
* If this method fails then use D0C
#define D0_DD_MASS

*#define DEBUG_Dcoeffs
**#define CHECK_Root_Kam

c#define DEBUG_CYPD0ir
c#define DEBUG_CYPD0M0


************************************************************************
	subroutine CYPDcoeff(para, D, ldpara)
************************************************************************
* Method: If abs(detG) is large enough we use double precision,
* otherwise the coefficients D(ddIJK...) are calculated 
* in quadruple precision.   
************************************************************************
	implicit none
	integer ldpara

	CVAR para(ldpara,Pdd)
	CVAR paraC(Pdd)
	CVAR G(3,3),detG
	CVAR D(Ndd)
	CVAR32 paraCQP(Pdd)
	CVAR32 DQP(Ndd)
	integer i
	logical dump

#include "lt.h"

***[DEBUG
	dump = ibits(debugkey, DebugD, 1) .ne. 0
	if( dump ) call CYPDumpPara(4, para, ldpara, "CYPDcoeff")
***]DEBUG

	do i = 1, Ndd
	 D(i) = 0
	enddo

* Ninh check
c	return
	
	do i = 1, Pdd
	 paraC(i) = para(1, i)
	enddo
	
***[set small masses to 0 when it's safe
	call CYPSetM0(paraC, Pdd)
***]set small masses to 0 when it's safe	
	
	G(1,1) = 2*paraC(1)
	G(2,2) = 2*paraC(5)
	G(3,3) = 2*paraC(4)
	G(1,2) = paraC(1) + paraC(5) - paraC(2)
	G(2,1) = G(1,2)
	G(1,3) = paraC(1) + paraC(4) - paraC(6)
	G(3,1) = G(1,3)
	G(2,3) = paraC(5) - paraC(3) + paraC(4)
	G(3,2) = G(2,3)	

	detG=G(1,1)*G(2,2)*G(3,3)+2*G(1,2)*G(2,3)*G(1,3)
     >  -G(1,3)*G(2,2)*G(1,3)-G(1,1)*G(2,3)*G(2,3)
     >  -G(3,3)*G(1,2)*G(1,2)
     
        pmax=Max(abs(paraC(1)),abs(paraC(2)),abs(paraC(3)),
     >	abs(paraC(4)),abs(paraC(5)),abs(paraC(6)))
	
        if ( abs(detG) .gt. stab4*8*pmax**3 .or. abs(detG) .eq. 0 ) then
***[using default precision
	call CYPDscalar(paraC, D(dd0))	
	call CYPDtensor(paraC, D)
***]using default precision
	else
***[using quadruple precision
	do i = 1, Pdd
	 paraCQP(i) = paraC(i)
	enddo	

	call DscalarCQP(paraCQP, DQP(dd0))	
	call DtensorCQP(paraCQP, DQP)

	do i = 1, Ndd
	 D(i) = DQP(i)
	enddo	
***]using quadruple precision
	endif

***[DEBUG
	if( dump ) call CYPDumpCoeff(4, D)
***]DEBUG

	end
***********************************************************************
	subroutine CYPDscalar(para, D0)
***********************************************************************
* Calculate the scalar 4-point integrals.
* Ref0. Dao Thi Nhung and Le Duc Ninh arXiv: 0902.0325
* Ref1. Denner and Dittmaier arXiv: 1005.2076
* Decision tree: 
*	if (soft-coll with quasi)
*	if (soft-coll)
*	if (finite with quasi)
*	if (finite with zero masses)
*	if (finite)
* todo: D0M0 (rotate if a = 0), D0ir (more quasi cases)
* Author: Le Duc Ninh (2011)
***********************************************************************
	implicit none
	CVAR para(Pdd)
        CVAR p1, p2, p3, p4, p1p2, p2p3
        CVAR m1, m2, m3, m4
	CVAR D0 
	CVAR ctemp
        CVAR cd,ce,ck,cfx
        RVAR ra,rb,rg,rc,rh,rj,signf
	RVAR pMmmax, pMmmid
        RVAR ap(6), am(4), apMm(12), am2Mm1, am3Mm1, am4Mm1, am3Mm2, 
     &  am4Mm2, am4Mm3
	RVAR ap1Mm1, ap4Mm1, ap1p2Mm1, ap1Mm2, ap2Mm2, ap2p3Mm2, 
     &  ap2Mm3, ap3Mm3, ap1p2Mm3, ap3Mm4, ap4Mm4, ap2p3Mm4
	equivalence (ap1Mm1, apMm(1)), (ap4Mm1, apMm(2)), 
     &  (ap1p2Mm1, apMm(3)), 
     &	(ap1Mm2, apMm(4)), (ap2Mm2, apMm(5)), (ap2p3Mm2, apMm(6)),
     &	(ap2Mm3, apMm(7)), (ap3Mm3, apMm(8)), (ap1p2Mm3, apMm(9)),
     &	(ap3Mm4, apMm(10)), (ap4Mm4, apMm(11)), (ap2p3Mm4, apMm(12))

	integer i, irot, mbase(4), mcase, ii(10)
	integer isoft, icoll, iqua, icolq
	integer fail
	data mbase /8, 4, 2, 1/
#include "lt.h"

***[versions
	integer key
	logical check
	key = ibits(versionkey, KeyD0, 2)
	check = .FALSE.
***]versions
	
	D0 = 0
	fail = 1
	iers = 5

c	print*,'D0.F hello, key =', key

	am(1) = abs(para(7))
	am(2) = abs(para(8))
	am(3) = abs(para(9))
	am(4) = abs(para(10))

***[detect zero masses and sort
* default mass positions:
* m1m2m3m4: 0000 0001 0011 0111 1111
* mcase   : 40   30   20   10   0
* binary trick: 8(1000) + 4(0100) + 2(0010) + 1(0001) = 15(1111)
	mcase = 0
	do i = 1, 4
         if ( am(i) .gt. 0) mcase = mcase + mbase(i)
        enddo
* move zero masses to the front
        call CYPD0per(mcase, ii)
	p1 = para(ii(1))
	p2 = para(ii(2))
	p3 = para(ii(3))
	p4 = para(ii(4))
	p1p2 = para(ii(5))
	p2p3 = para(ii(6))
	m1 = para(ii(7))
	m2 = para(ii(8))
	m3 = para(ii(9))
	m4 = para(ii(10))
***]detect zero masses and sort

	am(1) = abs(m1)
	am(2) = abs(m2)
	am(3) = abs(m3)
	am(4) = abs(m4)

	ap(1) = abs(p1)
	ap(2) = abs(p2)
	ap(3) = abs(p3)
	ap(4) = abs(p4)
	ap(5) = abs(p1p2)
	ap(6) = abs(p2p3)

	am2Mm1 = abs(m2 - m1)
	am3Mm1 = abs(m3 - m1)
	am4Mm1 = abs(m4 - m1)
	am3Mm2 = abs(m3 - m2)
	am4Mm2 = abs(m4 - m2)
	am4Mm3 = abs(m4 - m3)

***[These are all possibilities 
	apMm(1) = abs(p1 - m1)
	apMm(2) = abs(p4 - m1)
	apMm(3) = abs(p1p2 - m1)

	apMm(4) = abs(p1 - m2)
	apMm(5) = abs(p2 - m2)
	apMm(6) = abs(p2p3 - m2)

	apMm(7) = abs(p2 - m3)
	apMm(8) = abs(p3 - m3)
	apMm(9) = abs(p1p2 - m3)

	apMm(10) = abs(p3 - m4)
	apMm(11) = abs(p4 - m4)
	apMm(12) = abs(p2p3 - m4)
***]These are all possibilities
c	pMmmin = min(apMm(1),apMm(2),apMm(3),apMm(4),apMm(5),apMm(6),
c     &  apMm(7),apMm(8),apMm(9),apMm(10),apMm(11),apMm(12))

	pMmmax = max(apMm(1),apMm(2),apMm(3),apMm(4),apMm(5),apMm(6),
     &  apMm(7),apMm(8),apMm(9),apMm(10),apMm(11),apMm(12))
	pMmmid = pMmmax
	do i = 1, 12
	 if ( apMm(i) .lt. pMmmid .and. apMm(i) .gt. stab0*pMmmax ) 
     &   pMmmid = apMm(i)
	enddo

* Ninh checking (these are needed?)
	mmin = min(am(1),am(2),am(3),am(4))
	mmax = max(am(1),am(2),am(3),am(4))
c	mmid = mmax
c	do i = 1, 4
c	 if ( am(i) .lt. mmid .and. am(i) .gt. 0 ) mmid = am(i)
c	enddo

***[count number of soft and collinear singularities
	isoft = 0
	if( am(1) .le. 0 ) then
	  if( ap1Mm2 + ap4Mm4 .le. 0 .or. 
     &   ap1Mm2 + ap1p2Mm3 .le. 0 .or. 
     &   ap4Mm4 + ap1p2Mm3 .le. 0 ) isoft = isoft + 1
	endif

	if( am(2) .eq. 0 ) then
	  if( ap1Mm1 + ap2Mm3 .le. 0 .or. 
     &    ap1Mm1 + ap2p3Mm4 .le. 0 .or. 
     &    ap2Mm3 + ap2p3Mm4 .le. 0 ) isoft = isoft + 1
	endif
	
	if( am(3) .eq. 0 ) then
	  if( ap2Mm2 + ap3Mm4 .le. 0 .or. 
     &   ap2Mm2 + ap1p2Mm1 .le. 0 .or. 
     &   ap3Mm4 + ap1p2Mm1 .le. 0 ) isoft = isoft + 1
	endif
	
	if( am(4) .eq. 0 ) then
	  if( ap4Mm1 + ap3Mm3 .le. 0 .or. 
     &   ap4Mm1 + ap2p3Mm2 .le. 0 .or. 
     &   ap3Mm3 + ap2p3Mm2 .le. 0 ) isoft = isoft + 1
	endif

	icoll = 0
	iqua = 0
      if (am(1) + am(2) + ap(1) .le. 0) then
	icoll = icoll + 1
      else if (am(1) + am(2) + ap(1) .le. 2*stab0*pMmmid) then
	iqua = iqua + 1
      endif

      if (am(1) + am(3) + ap(5) .le. 0) then
	icoll = icoll + 1
      else if (am(1) + am(3) + ap(5) .le. 2*stab0*pMmmid) then
	iqua = iqua + 1
      endif

      if (am(1) + am(4) + ap(4) .le. 0) then
	icoll = icoll + 1
      else if (am(1) + am(4) + ap(4) .le. 2*stab0*pMmmid) then
	iqua = iqua + 1
      endif

      if (am(2) + am(3) + ap(2) .le. 0) then
	icoll = icoll + 1
      else if (am(2) + am(3) + ap(2) .le. 2*stab0*pMmmid) then
	iqua = iqua + 1
      endif

      if (am(2) + am(4) + ap(6) .le. 0) then
	icoll = icoll + 1
      else if (am(2) + am(4) + ap(6) .le. 2*stab0*pMmmid) then
	iqua = iqua + 1
      endif

      if (am(3) + am(4) + ap(3) .le. 0) then
	icoll = icoll + 1
      else if (am(3) + am(4) + ap(3) .le. 2*stab0*pMmmid) then
	iqua = iqua + 1
      endif
      icolq = icoll + iqua 
***]count number of soft and collinear singularities

***[detect quasi singularity: strong necessary conditions
	mcase = 0
	do i = 1, 4
         if ( am(i) .le. 0 ) then
	  mcase = mcase + 10
	 elseif ( am(i) .le. stab0*pMmmid ) then
	  mcase = mcase + 1
	 endif
        enddo
***]detect quasi singularity: strong necessary conditions

***[soft-coll (with quasi)
	if ( isoft + icolq .gt. 0 .and. iqua .gt. 0 ) then

* Ref1. (4.72): 3-soft 4-coll (2-quasi)
	if ( mcase .eq. 31 .and. isoft .eq. 3 .and. icolq .eq. 4 ) then
* default mass position: 0001
* possible rotations: 1, 3, 7, 9=3, 13=7, 15=1
	 if ( ap(1) + ap(2) .le. 0 .and. ap3Mm4 + ap4Mm4 .le. 0 ) then 
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,1)

 	  call CYPD0ir_3m0q_2p0a(fail,p3,p1p2,p2p3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif

	if ( ap(5) + ap(2) .le. 0 .and. ap2p3Mm4 + ap4Mm4 .le. 0 ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,3)

 	  call CYPD0ir_3m0q_2p0a(fail,p3,p1p2,p2p3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif

	if ( ap(1) + ap(5) .le. 0 .and. ap3Mm4 + ap2p3Mm4 .le. 0 ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,7)

 	  call CYPD0ir_3m0q_2p0a(fail,p3,p1p2,p2p3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif
	endif

* Ref1. (4.61): 2-soft 3-coll (1-quasi)
	if ( mcase .eq. 31 .and. isoft .eq. 2 .and. icolq .eq. 3 ) then
* default mass position: 0001
* possible rotations: 1, 3, 7, 9, 13, 15
	 if ( ap(1) + ap(2) .le. 0 .and. ap4Mm4 .le. 0 ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,1)

 	  call CYPD0ir_3m0q_2p0b(fail,p3,p4,p1p2,p2p3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap(5) + ap(2) .le. 0 .and. ap4Mm4 .le. 0 ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,3)

 	  call CYPD0ir_3m0q_2p0b(fail,p3,p4,p1p2,p2p3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap(1) + ap(5) .le. 0 .and. ap2p3Mm4 .le. 0 ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,7)

 	  call CYPD0ir_3m0q_2p0b(fail,p3,p4,p1p2,p2p3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap(2) + ap(5) .le. 0 .and. ap2p3Mm4 .le. 0 ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,9)

 	  call CYPD0ir_3m0q_2p0b(fail,p3,p4,p1p2,p2p3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap(5) + ap(1) .le. 0 .and. ap3Mm4 .le. 0 ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,13)

 	  call CYPD0ir_3m0q_2p0b(fail,p3,p4,p1p2,p2p3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap(1) + ap(2) .le. 0 .and. ap3Mm4 .le. 0 ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,15)

 	  call CYPD0ir_3m0q_2p0b(fail,p3,p4,p1p2,p2p3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif
	endif

* Ref1. (4.47): 2-soft 2-coll (2-quasi)
	if ( mcase .eq. 21 .and. isoft .eq. 2 .and. icolq .eq. 2 .and. 
     &  ap2Mm3 + ap1p2Mm3 + ap4Mm4 + ap2p3Mm4 .le. 0 ) then
* default mass position: 0101
* here p2 = p1p2 = m3, p4 = p2p3 = m4
* possible rotations: 3, 5, (9=3, 11=5)
	if ( ap(5) .le. stab0*pMmmid ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,3)

	 call CYPD0ir_2m0q_0p0a(fail,p1,p3,p1p2,p2p3,D0)
         if ( fail .eq. 0 ) goto 412
	endif

	if ( ap(4) .le. stab0*pMmmid ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,5)

	 call CYPD0ir_2m0q_0p0a(fail,p1,p3,p1p2,p2p3,D0)
         if ( fail .eq. 0 ) goto 412
	endif
	endif

* Ref1. (4.75): 1-soft 4-coll (4-quasi) 
	if ( mcase .eq. 13 .and. isoft .eq. 1 .and. icolq .eq. 4
     &  .and. am3Mm2 + am4Mm3 .le. 0 ) then
* default mass position: 0111
* possible rotations: 1, 2, 3, (4=2, 5=3, 6=1)
	 if ( ap(2) + ap(3) .le. 0 .and. ap1Mm2 + ap4Mm4 .le. 0 ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,1)

 	  call CYPD0ir_1m0q_2p0a(fail,p1,p1p2,p2p3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif

	if ( ap(6) + ap(3) .le. 0 .and. ap1Mm2 + ap1p2Mm3 .le. 0 ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,2)

 	  call CYPD0ir_1m0q_2p0a(fail,p1,p1p2,p2p3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif

	if ( ap(2) + ap(6) .le. 0 .and. ap1p2Mm3 + ap4Mm4 .le. 0 ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,3)

 	  call CYPD0ir_1m0q_2p0a(fail,p1,p1p2,p2p3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif
	endif

* Ref1. (4.39): 1-soft 2-coll (2-quasi)
	if ( mcase .eq. 21 .and. isoft .eq. 1 .and. icolq .eq. 2 ) then
* default mass position: 0101
* possible rotations: 3, 5, 9, 11
	if ( ap1p2Mm3 + ap2Mm3 + ap4Mm4 .le. 0 .and. 
     &  ap(5) .le. stab0*pMmmid ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,3)
 	  call CYPD0ir_2m0q_0p0b(fail,p1,p3,p4,p1p2,p2p3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif

	if ( ap4Mm4 + ap2p3Mm4 + ap1p2Mm3 .le. 0 .and. 
     &  ap(4) .le. stab0*pMmmid ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,5)
 	  call CYPD0ir_2m0q_0p0b(fail,p1,p3,p4,p1p2,p2p3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif

	if ( ap2Mm3 + ap1p2Mm3 + ap2p3Mm4 .le. 0 .and. 
     &  ap(2) .le. stab0*pMmmid ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,9)
 	  call CYPD0ir_2m0q_0p0b(fail,p1,p3,p4,p1p2,p2p3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif

	if ( ap2p3Mm4 + ap4Mm4 + ap2Mm3 .le. 0 .and. 
     &  ap(6) .le. stab0*pMmmid ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,11)
 	  call CYPD0ir_2m0q_0p0b(fail,p1,p3,p4,p1p2,p2p3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif
	endif

* Ref1. (4.31): 1-soft 2-coll (2-quasi) 
* for m3 = 0
	if ( mcase .eq. 22 .and. isoft .eq. 1 .and. icolq .eq. 2 ) then
* default mass position: 0101
* possible rotations: 3, 5=3, 9, 11=9
	 if ( ap1p2Mm3 + ap4Mm4 .le. 0 ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,3)

 	  call CYPD0ir_1m0q_0p0a(fail,p1,p2,p3,p4,p1p2,p2p3,m3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif

	if ( ap2Mm3 + ap2p3Mm4 .le. 0 ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,9)
 	  call CYPD0ir_1m0q_0p0a(fail,p1,p2,p3,p4,p1p2,p2p3,m3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif
	endif

* for m3 > 0
	if ( mcase .eq. 12 .and. isoft .eq. 1 .and. icolq .eq. 2 ) then
* default mass position: 0111
* possible rotations: 1, 2, 3, (4=2, 5=3, 6=1)
	 if ( ap1Mm2 + ap4Mm4 .le. 0 .and. 
     &   am(2) + am(4) .le. 2*stab0*pMmmid ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,1)

 	  call CYPD0ir_1m0q_0p0a(fail,p1,p2,p3,p4,p1p2,p2p3,m3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap1Mm2 + ap1p2Mm3 .le. 0 .and. 
     &   am(2) + am(3) .le. 2*stab0*pMmmid ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,2)

 	  call CYPD0ir_1m0q_0p0a(fail,p1,p2,p3,p4,p1p2,p2p3,m3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap1p2Mm3 + ap4Mm4 .le. 0 .and. 
     &   am(3) + am(4) .le. 2*stab0*pMmmid ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,3)

 	  call CYPD0ir_1m0q_0p0a(fail,p1,p2,p3,p4,p1p2,p2p3,m3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif
	endif

* Ref1. (4.20): 1-soft 1-coll (1-quasi)
* for m3 = 0
	if ( mcase .eq. 21 .and. isoft .eq. 1 .and. icolq .eq. 1 ) then
* default mass position: 0101
* possible rotations: 3, 5, 9, 11
	 if ( ap1p2Mm3 + ap4Mm4 .le. 0 .and. 
     &  ap(5) .le. stab0*pMmmid ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,3)

 	  call CYPD0ir_1m0q_0p0b(fail,p1,p2,p3,p4,p1p2,p2p3,m3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif

	if ( ap4Mm4 + ap1p2Mm3 .le. 0 .and. 
     &  ap(4) .le. stab0*pMmmid ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,5)
 	  call CYPD0ir_1m0q_0p0b(fail,p1,p2,p3,p4,p1p2,p2p3,m3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif

	if ( ap2Mm3 + ap2p3Mm4 .le. 0 .and. 
     &  ap(2) .le. stab0*pMmmid ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,9)
 	  call CYPD0ir_1m0q_0p0b(fail,p1,p2,p3,p4,p1p2,p2p3,m3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif

	if ( ap2p3Mm4 + ap2Mm3 .le. 0 .and. 
     &  ap(6) .le. stab0*pMmmid ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,11)
 	  call CYPD0ir_1m0q_0p0b(fail,p1,p2,p3,p4,p1p2,p2p3,m3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif
	endif 

* for m3 > 0
	if ( mcase .eq. 11 .and. isoft .eq. 1 .and. icolq .eq. 1 ) then
* default mass position: 0111
* possible rotations: 1, 2, 3, 4, 5, 6
	 if ( ap1Mm2 + ap4Mm4 .le. 0 .and. 
     &   am(2) .le. stab0*pMmmid ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,1)

 	  call CYPD0ir_1m0q_0p0b(fail,p1,p2,p3,p4,p1p2,p2p3,m3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap1Mm2 + ap1p2Mm3 .le. 0 .and. 
     &   am(2) .le. stab0*pMmmid ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,2)

 	  call CYPD0ir_1m0q_0p0b(fail,p1,p2,p3,p4,p1p2,p2p3,m3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif

	if ( ap1p2Mm3 + ap4Mm4 .le. 0 .and. 
     &   am(3) .le. stab0*pMmmid ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,3)

 	  call CYPD0ir_1m0q_0p0b(fail,p1,p2,p3,p4,p1p2,p2p3,m3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif

	if ( ap1p2Mm3 + ap1Mm2 .le. 0 .and. 
     &   am(3) .le. stab0*pMmmid ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,4)

 	  call CYPD0ir_1m0q_0p0b(fail,p1,p2,p3,p4,p1p2,p2p3,m3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif

	if ( ap4Mm4 + ap1p2Mm3 .le. 0 .and. 
     &   am(4) .le. stab0*pMmmid ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,5)

 	  call CYPD0ir_1m0q_0p0b(fail,p1,p2,p3,p4,p1p2,p2p3,m3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif

	if ( ap4Mm4 + ap1Mm2 .le. 0 .and. 
     &   am(4) .le. stab0*pMmmid ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,6)

 	  call CYPD0ir_1m0q_0p0b(fail,p1,p2,p3,p4,p1p2,p2p3,m3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif
	endif

	endif
***]soft-coll (with quasi)

***[soft-coll
	if ( isoft + icoll .gt. 0 ) then

c Ninh checking
C CYPD0per(mcase, ii) not working correctly if the masses are already in correct position 0011.
***[already done above, redundant?
#ifdef NINH_REDUNDANT_XXX
	mcase = 0
	do i = 1, 4
         if ( am(i) .gt. 0) mcase = mcase + mbase(i)
        enddo
* move zero masses to the front
        call CYPD0per(mcase, ii)
	p1 = para(ii(1))
	p2 = para(ii(2))
	p3 = para(ii(3))
	p4 = para(ii(4))
	p1p2 = para(ii(5))
	p2p3 = para(ii(6))
	m1 = para(ii(7))
	m2 = para(ii(8))
	m3 = para(ii(9))
	m4 = para(ii(10))

	mcase = 0
	do i = 1, 4
         if ( am(i) .le. 0 ) mcase = mcase + 10
        enddo
***]already done above, redundant?
#endif

* Ref1. (4.71): 4-soft 4-coll (1)
	if ( mcase .eq. 40 .and. isoft .eq. 4 .and. icoll .eq. 4 ) then
* default mass position: 0000
* possible rotations: 1, 2, 3, (4 = 2, 5 = 3, 6 = 1, 7 = 2 ..., 24 = 1)
         if ( ap(1) + ap(2) + ap(3) + ap(4) .le. 0 ) then
* WARNING: pi and mi change after this call.
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,1)

         call CYPD0ir_4m0_4p0a(fail,p1p2,p2p3,D0)
         if ( fail .eq. 0 ) goto 412
         endif

         if ( ap(1) + ap(5) + ap(3) + ap(6) .le. 0 ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,2)

         call CYPD0ir_4m0_4p0a(fail,p1p2,p2p3,D0)
         if ( fail .eq. 0 ) goto 412
         endif

         if ( ap(2) + ap(5) + ap(4) + ap(6) .le. 0 ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,3)

         call CYPD0ir_4m0_4p0a(fail,p1p2,p2p3,D0)
         if ( fail .eq. 0 ) goto 412
         endif
	endif

* Ref1. (4.45): 3-soft 2-coll (2)
	if ( mcase .eq. 30 .and. isoft .eq. 3 .and. icoll .eq. 2 ) then
* default mass position: 0001
* possible rotations: 1, 3, 7, (9=3, 13=7, 15=1)
         if ( ap(1) + ap(2) .le. 0 .and. ap4Mm4 + ap3Mm4 .le. 0 ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,1)

         call CYPD0ir_3m0_2p0a(fail,p3,p1p2,p2p3,D0)
         if ( fail .eq. 0 ) goto 412
         endif

        if ( ap(5) + ap(2) .le. 0 .and. ap4Mm4 + ap2p3Mm4 .le. 0 ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,3)

         call CYPD0ir_3m0_2p0a(fail,p3,p1p2,p2p3,D0)
         if ( fail .eq. 0 ) goto 412
         endif

        if ( ap(5) + ap(1) .le. 0 .and. ap3Mm4 + ap2p3Mm4 .le. 0 ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,7)

         call CYPD0ir_3m0_2p0a(fail,p3,p1p2,p2p3,D0)
         if ( fail .eq. 0 ) goto 412
         endif
	endif

* Ref1. (4.60): 2-soft 3-coll (3)
	if ( mcase .eq. 40 .and. isoft .eq. 2 .and. icoll .eq. 3 ) then
* default mass position: 0000
* possible rotations: 1 ... 24 (24 = 1, 23 = 7, ...)
         if ( ap(1) + ap(2) + ap(3) .le. 0 ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,1)

         call CYPD0ir_4m0_3p0a(fail,p4,p1p2,p2p3,D0)
         if ( fail .eq. 0 ) goto 412
	 endif

         if ( ap(1) + ap(4) + ap(3) .le. 0 ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,8)

         call CYPD0ir_4m0_3p0a(fail,p4,p1p2,p2p3,D0)
         if ( fail .eq. 0 ) goto 412
	 endif

         if ( ap(1) + ap(5) + ap(3) .le. 0 ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,7)

         call CYPD0ir_4m0_3p0a(fail,p4,p1p2,p2p3,D0)
         if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap(1) + ap(6) + ap(3) .le. 0 ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,2)

         call CYPD0ir_4m0_3p0a(fail,p4,p1p2,p2p3,D0)
         if ( fail .eq. 0 ) goto 412
	 endif

         if ( ap(2) + ap(1) + ap(4) .le. 0 ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,15)

         call CYPD0ir_4m0_3p0a(fail,p4,p1p2,p2p3,D0)
         if ( fail .eq. 0 ) goto 412
	 endif

         if ( ap(2) + ap(3) + ap(4) .le. 0 ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,10)

         call CYPD0ir_4m0_3p0a(fail,p4,p1p2,p2p3,D0)
         if ( fail .eq. 0 ) goto 412
	 endif

         if ( ap(2) + ap(5) + ap(4) .le. 0 ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,9)

         call CYPD0ir_4m0_3p0a(fail,p4,p1p2,p2p3,D0)
         if ( fail .eq. 0 ) goto 412
	 endif

         if ( ap(2) + ap(6) + ap(4) .le. 0 ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,16)

         call CYPD0ir_4m0_3p0a(fail,p4,p1p2,p2p3,D0)
         if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap(5) + ap(1) + ap(6) .le. 0 ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,13)

         call CYPD0ir_4m0_3p0a(fail,p4,p1p2,p2p3,D0)
         if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap(5) + ap(2) + ap(6) .le. 0 ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,3)

         call CYPD0ir_4m0_3p0a(fail,p4,p1p2,p2p3,D0)
         if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap(5) + ap(3) + ap(6) .le. 0 ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,12)

         call CYPD0ir_4m0_3p0a(fail,p4,p1p2,p2p3,D0)
         if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap(5) + ap(4) + ap(6) .le. 0 ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,11)

         call CYPD0ir_4m0_3p0a(fail,p4,p1p2,p2p3,D0)
         if ( fail .eq. 0 ) goto 412
	 endif
	endif

* Ref1. (4.36): 2-soft 2-coll (4)
	if ( mcase .eq. 30 .and. isoft .eq. 2 .and. icoll .eq. 2 ) then
* default mass position: 0001
* possible rotations: 1, 3, 7, 9, 13, 15
	 if ( ap(1) + ap(2) .le. 0 .and. ap4Mm4 .le. 0 ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,1)

         call CYPD0ir_3m0_2p0b(fail,p3,p4,p1p2,p2p3,D0)
         if ( fail .eq. 0 ) goto 412
	 endif

         if ( ap(5) + ap(2) .le. 0 .and. ap4Mm4 .le. 0 ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,3)

         call CYPD0ir_3m0_2p0b(fail,p3,p4,p1p2,p2p3,D0)
         if ( fail .eq. 0 ) goto 412
         endif

         if ( ap(5) + ap(1) .le. 0 .and. ap2p3Mm4 .le. 0 ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,7)

         call CYPD0ir_3m0_2p0b(fail,p3,p4,p1p2,p2p3,D0)
         if ( fail .eq. 0 ) goto 412
         endif

	 if ( ap(2) + ap(5) .le. 0 .and. ap2p3Mm4 .le. 0 ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,9)

         call CYPD0ir_3m0_2p0b(fail,p3,p4,p1p2,p2p3,D0)
         if ( fail .eq. 0 ) goto 412
         endif

         if ( ap(5) + ap(1) .le. 0 .and. ap3Mm4 .le. 0 ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,13)

         call CYPD0ir_3m0_2p0b(fail,p3,p4,p1p2,p2p3,D0)
         if ( fail .eq. 0 ) goto 412
         endif

	 if ( ap(1) + ap(2) .le. 0 .and. ap3Mm4 .le. 0 ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,15)

         call CYPD0ir_3m0_2p0b(fail,p3,p4,p1p2,p2p3,D0)
         if ( fail .eq. 0 ) goto 412
	 endif
	endif

* Ref1. (4.25): 2-soft 1-coll (5)
	if ( mcase .eq. 20 .and. isoft .eq. 2 .and. icoll .eq. 1 ) then
* default mass position: 0011
* possible rotations: 1, 2, (7=2, 8=1)
	if ( ap(1) .le. 0 .and. ap2Mm3 + ap4Mm4 .le. 0 ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,1)

	 call CYPD0ir_2m0_1p0a(fail,p2,p3,p4,p1p2,p2p3,D0)
         if ( fail .eq. 0 ) goto 412
	endif

	if ( ap(1) .le. 0 .and. ap2p3Mm4 + ap1p2Mm3 .le. 0 ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,2)

	 call CYPD0ir_2m0_1p0a(fail,p2,p3,p4,p1p2,p2p3,D0)
         if ( fail .eq. 0 ) goto 412
	endif
	endif

* Ref1. (4.9): 2-soft 0-coll (6)
	if ( mcase .eq. 20 .and. isoft .eq. 2 .and. icoll .eq. 0 ) then
* default mass position: 0101
* possible rotations: 3, (5=3, 9=3, 11=3)
	if ( ap2Mm3 + ap1p2Mm3 + ap4Mm4 + ap2p3Mm4 .le. 0 ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,3)

	 call CYPD0ir_2m0_0p0a(fail,p1,p3,p1p2,p2p3,D0)
         if ( fail .eq. 0 ) goto 412
	 endif
	endif

* Ref1. (4.29): 1-soft 2-coll (7)
* for m4 = 0
	if ( mcase .eq. 40 .and. isoft .eq. 1 .and. icoll .eq. 2 ) then
* default mass position: 0000
* possible rotations: 1, ..., 24
	 if ( ap(1) + ap(2) .le. 0 ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,1)

	 call CYPD0ir_3m0_2p0c(fail,p3,p4,p1p2,p2p3,m4,D0)
         if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap(1) + ap(4) .le. 0 ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,8)

	 call CYPD0ir_3m0_2p0c(fail,p3,p4,p1p2,p2p3,m4,D0)
         if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap(1) + ap(5) .le. 0 ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,7)

	 call CYPD0ir_3m0_2p0c(fail,p3,p4,p1p2,p2p3,m4,D0)
         if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap(1) + ap(6) .le. 0 ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,2)

	 call CYPD0ir_3m0_2p0c(fail,p3,p4,p1p2,p2p3,m4,D0)
         if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap(2) + ap(3) .le. 0 ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,10)

	 call CYPD0ir_3m0_2p0c(fail,p3,p4,p1p2,p2p3,m4,D0)
         if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap(2) + ap(5) .le. 0 ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,9)

	 call CYPD0ir_3m0_2p0c(fail,p3,p4,p1p2,p2p3,m4,D0)
         if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap(2) + ap(6) .le. 0 ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,16)

	 call CYPD0ir_3m0_2p0c(fail,p3,p4,p1p2,p2p3,m4,D0)
         if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap(3) + ap(4) .le. 0 ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,6)

	 call CYPD0ir_3m0_2p0c(fail,p3,p4,p1p2,p2p3,m4,D0)
         if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap(3) + ap(5) .le. 0 ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,23)

	 call CYPD0ir_3m0_2p0c(fail,p3,p4,p1p2,p2p3,m4,D0)
         if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap(3) + ap(6) .le. 0 ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,18)

	 call CYPD0ir_3m0_2p0c(fail,p3,p4,p1p2,p2p3,m4,D0)
         if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap(4) + ap(5) .le. 0 ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,20)

	 call CYPD0ir_3m0_2p0c(fail,p3,p4,p1p2,p2p3,m4,D0)
         if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap(4) + ap(6) .le. 0 ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,11)

	 call CYPD0ir_3m0_2p0c(fail,p3,p4,p1p2,p2p3,m4,D0)
         if ( fail .eq. 0 ) goto 412
	 endif
	endif

* for m4 > 0
	if ( mcase .eq. 30 .and. isoft .eq. 1 .and. icoll .eq. 2 ) then
* default mass position: 0001
* possible rotations: 1, 3, 7, 9=3, 13=7, 15=1
	 if ( ap(1) + ap(2) .le. 0 ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,1)

	 call CYPD0ir_3m0_2p0c(fail,p3,p4,p1p2,p2p3,m4,D0)
         if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap(5) + ap(2) .le. 0 ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,3)

	 call CYPD0ir_3m0_2p0c(fail,p3,p4,p1p2,p2p3,m4,D0)
         if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap(1) + ap(5) .le. 0 ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,7)

	 call CYPD0ir_3m0_2p0c(fail,p3,p4,p1p2,p2p3,m4,D0)
         if ( fail .eq. 0 ) goto 412
	 endif
	endif

* Ref1. (4.19): 1-soft 1-coll (8)
* for m3 = 0
	if ( mcase .eq. 30 .and. isoft .eq. 1 .and. icoll .eq. 1 ) then
* default mass position: 0001
* possible rotations: 1, 3, 7, 9, 13, 15
	 if ( ap(1) .le. 0 .and. ap4Mm4 .le. 0 ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,1)

	  call CYPD0ir_2m0_1p0b(fail,p2,p3,p4,p1p2,p2p3,m3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap(5) .le. 0 .and. ap4Mm4 .le. 0 ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,3)

	  call CYPD0ir_2m0_1p0b(fail,p2,p3,p4,p1p2,p2p3,m3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap(1) .le. 0 .and. ap2p3Mm4 .le. 0 ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,7)

	  call CYPD0ir_2m0_1p0b(fail,p2,p3,p4,p1p2,p2p3,m3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap(2) .le. 0 .and. ap2p3Mm4 .le. 0 ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,9)

	  call CYPD0ir_2m0_1p0b(fail,p2,p3,p4,p1p2,p2p3,m3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap(5) .le. 0 .and. ap3Mm4 .le. 0 ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,13)

	  call CYPD0ir_2m0_1p0b(fail,p2,p3,p4,p1p2,p2p3,m3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap(2) .le. 0 .and. ap3Mm4 .le. 0 ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,15)

	  call CYPD0ir_2m0_1p0b(fail,p2,p3,p4,p1p2,p2p3,m3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif
	endif

* for m3 > 0
	if ( mcase .eq. 20 .and. isoft .eq. 1 .and. icoll .eq. 1 ) then
* default mass position: 0011
* possible rotations (m1 = m2 = 0): 1, 2, 7, 8
	 if ( ap(1) .le. 0 .and. ap4Mm4 .le. 0 ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,1)

	  call CYPD0ir_2m0_1p0b(fail,p2,p3,p4,p1p2,p2p3,m3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap(1) .le. 0 .and. ap1p2Mm3 .le. 0 ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,2)

	  call CYPD0ir_2m0_1p0b(fail,p2,p3,p4,p1p2,p2p3,m3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap(1) .le. 0 .and. ap2p3Mm4 .le. 0 ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,7)

	  call CYPD0ir_2m0_1p0b(fail,p2,p3,p4,p1p2,p2p3,m3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap(1) .le. 0 .and. ap2Mm3 .le. 0 ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,8)

	  call CYPD0ir_2m0_1p0b(fail,p2,p3,p4,p1p2,p2p3,m3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif
	endif

* Ref1. (4.4): 1-soft 0-coll (9)
* for m3 = 0
	if ( mcase .eq. 20 .and. isoft .eq. 1 .and. icoll .eq. 0 ) then
* default mass position: 0101
* possible rotations: 3, 5=3, 9, 11=9
	 if ( ap1p2Mm3 + ap4Mm4 .le. 0 ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,3)

 	  call CYPD0ir_1m0_0p0a(fail,p1,p2,p3,p4,p1p2,p2p3,m3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap2Mm3 + ap2p3Mm4 .le. 0 ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,9)

 	  call CYPD0ir_1m0_0p0a(fail,p1,p2,p3,p4,p1p2,p2p3,m3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif
	endif

* for m3 > 0
	if ( mcase .eq. 10 .and. isoft .eq. 1 .and. icoll .eq. 0 ) then
* default mass position: 0111
* possible rotations (m1 = 0): 1, 2, 3, (4=2, 5=3, 6=1)
	 if ( ap1Mm2 + ap4Mm4 .le. 0 ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,1)

 	  call CYPD0ir_1m0_0p0a(fail,p1,p2,p3,p4,p1p2,p2p3,m3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap1Mm2 + ap1p2Mm3 .le. 0 ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,2)

 	  call CYPD0ir_1m0_0p0a(fail,p1,p2,p3,p4,p1p2,p2p3,m3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap1p2Mm3 + ap4Mm4 .le. 0 ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,3)

 	  call CYPD0ir_1m0_0p0a(fail,p1,p2,p3,p4,p1p2,p2p3,m3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif
	endif

* Ref1. (4.52): 0-soft 2-coll (10)
	if ( mcase .eq. 40 .and. isoft .eq. 0 .and. icoll .eq. 2 ) then
* default mass position: 0000
* possible rotations (p1p3, p2p4, p5p6): 1, 3, 5
	 if ( ap(1) + ap(3) .le. 0 ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,1)

 	  call CYPD0ir_4m0_2p0a(fail,p2,p4,p1p2,p2p3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap(2) + ap(4) .le. 0 ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,5)

 	  call CYPD0ir_4m0_2p0a(fail,p2,p4,p1p2,p2p3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap(5) + ap(6) .le. 0 ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,3)

 	  call CYPD0ir_4m0_2p0a(fail,p2,p4,p1p2,p2p3,D0)
          if ( fail .eq. 0 ) goto 412
	 endif
	endif

* Ref1. (4.13): 0-soft 1-coll (11)
* for m3 = m4 = 0
	if ( mcase .eq. 40 .and. isoft .eq. 0 .and. icoll .eq. 1 ) then
* default mass position: 0000
* possible rotations: 1, ..., 24
	 if ( ap(1) .le. 0 ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,1)

 	  call CYPD0ir_2m0_1p0c(fail,p2,p3,p4,p1p2,p2p3,m3,m4,D0)
          if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap(2) .le. 0 ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,9)

 	  call CYPD0ir_2m0_1p0c(fail,p2,p3,p4,p1p2,p2p3,m3,m4,D0)
          if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap(3) .le. 0 ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,17)

 	  call CYPD0ir_2m0_1p0c(fail,p2,p3,p4,p1p2,p2p3,m3,m4,D0)
          if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap(4) .le. 0 ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,19)

 	  call CYPD0ir_2m0_1p0c(fail,p2,p3,p4,p1p2,p2p3,m3,m4,D0)
          if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap(5) .le. 0 ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,3)

 	  call CYPD0ir_2m0_1p0c(fail,p2,p3,p4,p1p2,p2p3,m3,m4,D0)
          if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap(6) .le. 0 ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,11)

 	  call CYPD0ir_2m0_1p0c(fail,p2,p3,p4,p1p2,p2p3,m3,m4,D0)
          if ( fail .eq. 0 ) goto 412
	 endif
	endif

* for m3 = 0, m4 > 0
	if ( mcase .eq. 30 .and. isoft .eq. 0 .and. icoll .eq. 1 ) then
* default mass position: 0001
* possible rotations: 1, 3, 7=1, 9, 13=3, 15=9
	 if ( ap(1) .le. 0 ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,1)

 	  call CYPD0ir_2m0_1p0c(fail,p2,p3,p4,p1p2,p2p3,m3,m4,D0)
          if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap(5) .le. 0 ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,3)

 	  call CYPD0ir_2m0_1p0c(fail,p2,p3,p4,p1p2,p2p3,m3,m4,D0)
          if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap(2) .le. 0 ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,9)

 	  call CYPD0ir_2m0_1p0c(fail,p2,p3,p4,p1p2,p2p3,m3,m4,D0)
          if ( fail .eq. 0 ) goto 412
	 endif
	endif

* for m3 > 0, m4 > 0
	if ( mcase .eq. 20 .and. isoft .eq. 0 .and. icoll .eq. 1 ) then
* default mass position: 0011
* possible rotations (m1 = m2 = 0): 1, (2=1, 7=1, 8=1)
	 if ( ap(1) .le. 0 ) then
          call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,1)

 	  call CYPD0ir_2m0_1p0c(fail,p2,p3,p4,p1p2,p2p3,m3,m4,D0)
          if ( fail .eq. 0 ) goto 412
	 endif
	endif

	endif
***]soft-coll

***[finite cases
	if ( isoft + icoll .gt. 0 ) then
	print*,"Error in Dcoeffs.F: IR case not done"
	print*,"isoft =", isoft
	print*,"icoll =", icoll
        print*,"iqua  =", iqua
        print*,"mcase =", mcase
        print*,"stab0 =", stab0
        print*,"pMmmid =", pMmmid
	print*,"p1    =", para(1)
	print*,"p2    =", para(2)
	print*,"p3    =", para(3)
	print*,"p4    =", para(4)
	print*,"p1p2  =", para(5)
	print*,"p2p3  =", para(6)
	print*,"m1    =", para(7)
	print*,"m2    =", para(8)
	print*,"m3    =", para(9)
	print*,"m4    =", para(10)
	stop

	else

	if ( lambda .lt. 0 ) goto 412

c	print*,'D0 fin hello 0'
c	print*, p1, p2, p3, p4, p1p2, p2p3, m1, m2, m3, m4

***[quasi finite
	if ( iqua .gt. 0 ) then

* Ref1. (4.76): 0-soft 4-coll (4-quasi)
	if ( mcase .eq. 4 .and. iqua .eq. 4 
     &  .and. am2Mm1 + am3Mm1 + am4Mm1 .le. 0 ) then
* default mass position: 1111
* possible rotations: 1, 2, 3, (4 = 2, 5 = 3, 6 = 1, 7 = 2 ..., 24 = 1)

         if ( ap(1) + ap(2) + ap(3) + ap(4) .le. 0 ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,1)

	 call CYPD0ir_0m0q_4p0a(fail,p1p2,p2p3,m1,D0)
         if ( fail .eq. 0 ) goto 412
         endif

         if ( ap(1) + ap(5) + ap(3) + ap(6) .le. 0 ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,2)

         call CYPD0ir_0m0q_4p0a(fail,p1p2,p2p3,m1,D0)
         if ( fail .eq. 0 ) goto 412
         endif

         if ( ap(2) + ap(5) + ap(4) + ap(6) .le. 0 ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,3)

         call CYPD0ir_0m0q_4p0a(fail,p1p2,p2p3,m1,D0)
         if ( fail .eq. 0 ) goto 412
         endif
	endif

* Ref1. (4.14): 0-soft 1-coll (1-quasi)
* for m3 = m4 = 0
	if ( mcase .eq. 31 .and. iqua .eq. 1 ) then
* default mass position: 0100
* possible rotations: 5, 6=5, 11, 12=11, 17, 18=17
	 if ( ap4Mm4 .le. 0 .and. am(4) .le. stab0*pMmmid ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,5)
         call CYPD0ir_1m0q_0p0c(fail,p1,p2,p3,p4,p1p2,p2p3,
     &  m3,m4,D0)
         if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap2p3Mm4 .le. 0 .and. am(4) .le. stab0*pMmmid ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,11)
         call CYPD0ir_1m0q_0p0c(fail,p1,p2,p3,p4,p1p2,p2p3,
     &  m3,m4,D0)
         if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap3Mm4 .le. 0 .and. am(4) .le. stab0*pMmmid ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,17)
         call CYPD0ir_1m0q_0p0c(fail,p1,p2,p3,p4,p1p2,p2p3,
     &  m3,m4,D0)
         if ( fail .eq. 0 ) goto 412
	 endif
	endif

* for m3 = 0
	if ( mcase .eq. 21 .and. iqua .eq. 1 ) then
* default mass position: 0101
* possible rotations: 3, 5, 9, 11
	 if ( ap1p2Mm3 .le. 0 .and. am(3) .le. stab0*pMmmid ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,3)
         call CYPD0ir_1m0q_0p0c(fail,p1,p2,p3,p4,p1p2,p2p3,
     &  m3,m4,D0)
         if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap4Mm4 .le. 0 .and. am(4) .le. stab0*pMmmid ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,5)
         call CYPD0ir_1m0q_0p0c(fail,p1,p2,p3,p4,p1p2,p2p3,
     &  m3,m4,D0)
         if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap2Mm3 .le. 0 .and. am(3) .le. stab0*pMmmid ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,9)
         call CYPD0ir_1m0q_0p0c(fail,p1,p2,p3,p4,p1p2,p2p3,
     &  m3,m4,D0)
         if ( fail .eq. 0 ) goto 412
	 endif

	 if ( ap2p3Mm4 .le. 0 .and. am(4) .le. stab0*pMmmid ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,11)
         call CYPD0ir_1m0q_0p0c(fail,p1,p2,p3,p4,p1p2,p2p3,
     &  m3,m4,D0)
         if ( fail .eq. 0 ) goto 412
	 endif
	endif

* for m3 > 0
	if ( mcase .eq. 11 .and. iqua .eq. 1 ) then
* possible rotations: 1, 2=1, 3, 4=3, 5, 6=5
         if ( ap1Mm2 .le. 0 .and. ap(1) .le. stab0*pMmmid ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,1)

         call CYPD0ir_1m0q_0p0c(fail,p1,p2,p3,p4,p1p2,p2p3,
     &  m3,m4,D0)
         if ( fail .eq. 0 ) goto 412
         endif

	 if ( ap1p2Mm3 .le. 0 .and. ap(5) .le. stab0*pMmmid ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,3)

         call CYPD0ir_1m0q_0p0c(fail,p1,p2,p3,p4,p1p2,p2p3,
     &  m3,m4,D0)
         if ( fail .eq. 0 ) goto 412
         endif

	 if ( ap4Mm4 .le. 0 .and. ap(4) .le. stab0*pMmmid ) then
         call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,5)

         call CYPD0ir_1m0q_0p0c(fail,p1,p2,p3,p4,p1p2,p2p3,
     &  m3,m4,D0)
         if ( fail .eq. 0 ) goto 412
         endif
	endif

	endif
***]quasi finite

***[finite with zero masses
	if ( mmin .le. 0 ) then

* Important:
* choose the first from the top of the symmetry table.
	irot = 1
	if ( abs(m1) + abs(m2) + abs(m3) + abs(m4) .le. 0 ) then
	
	elseif ( abs(m2) + abs(m3) + abs(m4) .le. 0 ) then
	
	elseif ( abs(m1) + abs(m3) + abs(m4) .le. 0 ) then
	 irot = 7 
	elseif ( abs(m1) + abs(m2) + abs(m4) .le. 0 ) then
	 irot = 13
	elseif ( abs(m1) + abs(m2) + abs(m3) .le. 0 ) then
	 irot = 19
	elseif ( abs(m2) + abs(m3) .le. 0 ) then
	
	elseif ( abs(m2) + abs(m4) .le. 0 ) then
	 irot = 2
	elseif ( abs(m1) + abs(m2) .le. 0 ) then
	 irot = 13
	elseif ( abs(m1) + abs(m3) .le. 0 ) then
	 irot = 7
	elseif ( abs(m1) + abs(m4) .le. 0 ) then
	 irot = 8
	elseif ( abs(m3) + abs(m4) .le. 0 ) then
	 irot = 4
	elseif ( abs(m3) .le. 0 ) then
	
	elseif ( abs(m4) .le. 0 ) then
	 irot = 2
	elseif ( abs(m1) .le. 0 ) then
	 irot = 9
	elseif ( abs(m2) .le. 0 ) then
	 irot = 3
	endif
	
* WARNING: pi and mi change after this call.
	 call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,irot)		
	
	call CYPD0M0(fail,p1,p2,p3,p4,p1p2,p2p3,m1,m2,m3,m4,D0)

* to cross-check with D0C
cc	if ( fail .eq. 0 ) then
cc	 ctemp = D0
cc	endif

	if ( fail .eq. 0 ) goto 412
	endif
***]finite with zero masses

***[finite without zero masses
	if( key .ne. 1 ) then
***[finite without zero masses, method of Ref1.
	call CYPRotateP5(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,irot)

	if ( irot .eq. 25 ) then 
	 print*,'Error in CYPDscalar(): 
     >	 all lambda(pi, pj, pk) < 0 not done.'
	 iers = 5**7
	 goto 412

	else
* WARNING: pi and mi change after this call.
	 call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,irot)

***[Ninh checking DD1
* The idea is to use Ref1. Eq. (3.61) as far as possible. 
* The point is that this Eq. is valid for real and positive r20, 
* which is usually(rarely) the case for real(complex) masses. 
* The problem is that the present implementation is only OK 
* for vanishing p1p2 (ie. r20 = 1).
* For zero p1p2, r31 (r20 in Ref1) is real and positive:
c	 call CYPD0M0(fail,p1,p2,p3,p4,p1p2,p2p3,m1,m2,m3,m4,D0)
c	 if ( fail .ne. 0 ) 
c     &   call CYPD0DD(fail,p1,p2,p3,p4,p1p2,p2p3,m1,m2,m3,m4,D0)
c	 if ( fail .eq. 0 ) goto 412
***]Ninh checking DD1

	if ( abs(p1p2) .le. 0 ) then
	 call CYPD0M0(fail,p1,p2,p3,p4,p1p2,p2p3,m1,m2,m3,m4,D0)
* Use Ref1. Eq. (3.65)
	else
         call CYPD0DD(fail,p1,p2,p3,p4,p1p2,p2p3,m1,m2,m3,m4,D0)
	endif

	if ( fail .eq. 0 ) then
	 ctemp = D0
	 if( key .eq. 0 ) then
	  goto 412
	 else
	  check = .TRUE.
	 endif
	endif

	endif
***]finite without zero masses, method of Ref1.
	endif

	if( key .ne. 0 .or. fail .ne. 0 ) then
***[finite without zero masses, method of Ref0.
	call CYPRotate(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,irot)

	if ( irot .eq. 25 ) then 
	 print*,'Error in CYPDscalar(): 
     >	 all lambda(pi, pj, pk) < 0 not done.'
	 iers = 5**7
	 goto 412

	else
* WARNING: pi and mi change after this call.
	 call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,irot)		

         ra = DBLE(p3)
         rb = DBLE(p2)
         rg = DBLE(p1)
         rc = DBLE(p2p3 - p2 - p3)
         rh = DBLE(p2 + p4 - p1p2 - p2p3)
         rj = DBLE(p1p2 - p1 - p2)
         cd = m3 - m4 - p3
         ce = m2 - m3 + p3 - p2p3
         ck = m1 - m2 + p2p3 - p4
         cfx = m4
	 signf = -1	 

         call CYPD0_abc(ra,rb,rg,rc,rh,rj,cd,ce,ck,cfx,signf,D0)

	 if( key .gt. 1 .and. check .and. 
     &        abs(D0 - ctemp) .gt. maxdev*abs(D0) ) then
	   print *, "Discrepancy in D0:"
	   print*,'D0DD  =',ctemp
	   print*,'D0abc =',D0
      	   print *, "    p1   =", para(1)
	   print *, "    p2   =", para(2)
	   print *, "    p3   =", para(3)
	   print *, "    p4   =", para(4)
	   print *, "    p1p2 =", para(5)
	   print *, "    p2p3 =", para(6)
	   print *, "    m1 =", para(7)
	   print *, "    m2 =", para(8)
	   print *, "    m3 =", para(9)
	   print *, "    m4 =", para(10)
c	   D0 = ctemp
	 endif	

	endif
***]finite without zero masses, method of Ref0.
	endif

***]finite without zero masses

c	print*,'From CYPDscalar(): D0 =', D0

	endif
***]finite cases

 412    continue
***[check for errors
	if ( mod(iers,25) .eq. 0 ) then 
	print *, "CYPDscalar in Dcoeffs.F:"
	print *, "    p1   =", para(1)
	print *, "    p2   =", para(2)
	print *, "    p3   =", para(3)
	print *, "    p4   =", para(4)
	print *, "    p1p2 =", para(5)
	print *, "    p2p3 =", para(6)
	print *, "    m1 =", para(7)
	print *, "    m2 =", para(8)
	print *, "    m3 =", para(9)
	print *, "    m4 =", para(10)
	print *, "    D0 =", D0
	endif
***]check for errors
	iers = 1
	end
************************************************************************
	subroutine CYPD0M0(fail,q1,q2,q3,q4,q1q2,q2q3,
     >   n1,n2,n3,n4,D0)
************************************************************************
* Calculate the finite scalar 4-point integrals with zero masses.
* Ref0. Denner and Dittmaier arXiv: 1005.2076
* todo: use rotation to avoid log(0) or 1/0 divergences, if it occurs.
* Author: Le Duc Ninh (2011)
************************************************************************
	implicit none
        CVAR q1, q2, q3, q4, q1q2, q2q3
        CVAR n1, n2, n3, n4
        CVAR p1, p2, p3, p4, p1p2, p2p3
        CVAR m1, m2, m3, m4
        CVAR rm1, rm2, rm3, rm4
	CVAR Y(4,4), a, b, c, d
	CVAR D0, ctemp
	CVAR x(2), cx(2), cr(2), dr(2), r24
	CVAR er(2), fr(2), r31, x31

	CVAR x14, x12, x34, x32, x24

	RVAR ix(2), icx(2), icr(2), idr(2), ir24
	RVAR ier(2), ifr(2), ir31, itemp, ix31, ix31s

	RVAR ix14, ix12, ix34, ix32, ix24
	RVAR ix14s, ix12s, ix34s, ix32s, ix24s

	integer k, l, irot, fail

#include "lt.h"

	D0 = 0
	fail = 0

	irot = 1
 101    continue

	p1 = q1
	p2 = q2
	p3 = q3
	p4 = q4
	p1p2 = q1q2
	p2p3 = q2q3
	m1 = n1
	m2 = n2
	m3 = n3
	m4 = n4

* WARNING: pi and mi change after this call.
	 call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,irot)

	rm1 = sqrt(m1)
	rm2 = sqrt(m2)
	rm3 = sqrt(m3)
	rm4 = sqrt(m4)

	Y(1,1) = 2*m1
	Y(2,2) = 2*m2
	Y(3,3) = 2*m3
	Y(4,4) = 2*m4
	Y(1,2) = m1 + m2 - p1
	Y(1,3) = m1 + m3 - p1p2
	Y(1,4) = m1 + m4 - p4 
	Y(2,3) = m2 + m3 - p2
	Y(2,4) = m2 + m4 - p2p3
	Y(3,4) = m3 + m4 - p3
	do k = 1, 3
	 do l = k + 1, 4
	  Y(l,k) = Y(k,l)
	 enddo
	enddo

***[4 zero masses
	if ( abs(m1) + abs(m2) + abs(m3) + abs(m4) .le. 0 ) then 

#ifdef DEBUG_CYPD0M0
	print*,'CYPD0M0: begin 4 zero masses'
#endif

* Ref0. Eq. (3.85)
	 a = Y(2,4)*Y(3,4)
	 b = Y(1,3)*Y(2,4) + Y(1,2)*Y(3,4) - Y(1,4)*Y(2,3)
	 c = Y(1,2)*Y(1,3) 
	 d = Y(2,3)
	 
* a, c and d cannot be zero here
* See Ref0. Eq. (3.49) 
* ax^2 + bx + c + i*eps*d = 0
	 call CYPRoot2(a, b, c, DBLE(d), x, ix)

* check for Landau singularity
	 if ( abs(x(1) - x(2)) .le. 0 ) then 
	 print*,'Error in CYPD0M0(): 4 zero mass and detY = 0'
	 print*,'todo: try rotation'
	 goto 109
	 endif

* Ref0. Eq. (3.84)
	 do k = 1,2 
	 cx(1) = -x(k)
	 icx(1) = -ix(k)
	 cx(2) = Y(3,4)/Y(1,3)
	 icx(2) = -DBLE(Y(1,3) - Y(3,4))
	 if ( icx(2) .eq. 0 ) icx(2) = 1D0

	 itemp = icx(1)*DBLE(cx(2))
	 if ( itemp .eq. 0 ) itemp = icx(1)

	 D0 = D0 - (-1D0)**k*CYPLi2nx(cx, icx, itemp, 2)
	 
	 cx(1) = -1/x(k)
	 icx(1) = ix(k)
	 cx(2) = Y(1,2)/Y(2,4)
	 icx(2) = -DBLE(Y(2,4) - Y(1,2))
	 if ( icx(2) .eq. 0 ) icx(2) = 1D0

	 itemp = icx(1)*DBLE(cx(2))
	 if ( itemp .eq. 0 ) itemp = icx(1)

	 D0 = D0 + (-1D0)**k*CYPLi2nx(cx, icx, itemp, 2)
	 
	 cx(1) = Y(2,4)/Y(2,3)
	 icx(1) = -DBLE(Y(2,3) - Y(2,4))
	 if ( icx(1) .eq. 0 ) icx(1) = 1D0

	 cx(2) = Y(1,3)/Y(1,4)
	 icx(2) = -DBLE(Y(1,4) - Y(1,3))
	 if ( icx(2) .eq. 0 ) icx(2) = 1D0
	 D0 = D0 + (-1D0)**k*CYPLOG(-x(k),-ix(k))*
     &	 (CYPLOG(cx(1),icx(1)) + CYPLOG(cx(2),icx(2)))
	 enddo

	 D0 = D0/a/(x(1) - x(2))

c	 print*,'finite D0 4M0:',D0
	 return
	endif
***]4 zero masses

***[3 zero masses
	if ( abs(m2) + abs(m3) + abs(m4) .le. 0 ) then 
* possible rotations: 1,2,3,4,5,6.
* check for log(0) singularity

#ifdef DEBUG_CYPD0M0
	print*,'CYPD0M0: begin 3 zero masses'
#endif

	 if ( abs(Y(1,2)) + abs(Y(1,3)) .le. 0 ) then 
	 print*,'Error in CYPD0M0(): 3 zero mass and log(0) 1'
	 goto 109
	 endif
	 
	 if ( abs(Y(1,2))*abs(Y(1,3)) .le. 0 ) then
	 if ( abs(m1 - p4) .le. 0 ) then 
	 print*,'Error in CYPD0M0(): 3 zero mass and log(0) 2'
	 goto 109
	 
	 elseif ( abs(Y(1,2)) .le. 0 ) then
	 irot = 4
	 goto 101
	 else
	 irot = 2
	 goto 101 
	 endif
	 endif 

* Ref0. Eq. (3.83)
	 a = Y(2,4)*Y(3,4)
	 b = Y(1,3)*Y(2,4) + Y(1,2)*Y(3,4) - Y(1,4)*Y(2,3)
	 c = Y(1,2)*Y(1,3) - m1*Y(2,3) 
	 d = Y(2,3)
	 
* Condition: Y12*Y13 .ne. 0 
* a and d cannot be zero here
* See Ref0. Eq. (3.49) 
* ax^2 + bx + c + i*eps*d = 0
	 call CYPRoot2(a, b, c, DBLE(d), x, ix)

* check for Landau singularity
	 if ( abs(x(1) - x(2)) .le. 0 ) then 
	 print*,'Error in CYPD0M0(): 3 zero mass and detY = 0'
	 goto 109
	 endif

* Ref0. Eq. (3.82)
	 do k = 1,2
	 cx(1) = -x(k)
	 icx(1) = -ix(k)
	 cx(2) = Y(3,4)/Y(1,3)
	 icx(2) = DBLE(Y(3,4))
	 if ( icx(2) .eq. 0 ) icx(2) = 1D0
	 itemp = icx(1)*DBLE(cx(2))
	 if ( itemp .eq. 0 ) itemp = icx(1)
	 D0 = D0 - (-1D0)**k*CYPLi2nx(cx, icx, itemp, 2)
	 
	 cx(2) = Y(1,4)/m1
	 icx(2) = -DBLE(m1 - Y(1,4))
	 if ( icx(2) .eq. 0 ) icx(2) = 1D0
	 itemp = icx(1)*DBLE(cx(2))
	 if ( itemp .eq. 0 ) itemp = icx(1)
	 D0 = D0 + (-1D0)**k*CYPLi2nx(cx, icx, itemp, 2)
	 
	 cx(2) = Y(2,4)/Y(1,2)
	 icx(2) = DBLE(Y(2,4))
	 if ( icx(2) .eq. 0 ) icx(2) = 1D0
	 itemp = icx(1)*DBLE(cx(2))
	 if ( itemp .eq. 0 ) itemp = icx(1)
	 D0 = D0 - (-1D0)**k*CYPLi2nx(cx, icx, itemp, 2)

	 
	 cx(1) = Y(1,2)/Y(2,3)
	 icx(1) = -DBLE(Y(2,3))
	 if ( icx(1) .eq. 0 ) icx(1) = 1D0

	 cx(2) = Y(1,3)/m1
	 icx(2) = -DBLE(m1 - Y(1,3))
	 if ( icx(2) .eq. 0 ) icx(2) = 1D0
* if x(k) = 0, ie. c = 0, then 2 cases are possible:
* either this log-term vanishes or diverge (3pt Landau). 
* We set the divergent part to 0 if it occurs.
	 if (abs(x(k)) .gt. 0) then 
	 D0 = D0 + (-1D0)**k*CYPLOG(-x(k),-ix(k))*
     &	 (CYPLOG(cx(1),icx(1)) + CYPLOG(cx(2),icx(2)))
	 endif
	 enddo
	 
	 D0 = D0/a/(x(1) - x(2))

c	 print*,'finite D0 3M0:',D0
	 return
	endif
***]3 zero masses

***[2 zero masses
	if ( abs(m2) + abs(m3) .le. 0 ) then 
* possible rotations: 1,3,22,24.
* check for log(0) singularity

#ifdef DEBUG_CYPD0M0
	print*,'CYPD0M0: begin 2 zero masses'
#endif

	 if ( abs(Y(1,2)) + abs(Y(1,3)) .le. 0 ) then 
	 irot = 22
	 goto 101
	 endif
	 
	 if ( abs(Y(1,2))*abs(Y(1,3)) .le. 0 ) then
	 if ( abs(Y(2,4))*abs(Y(3,4)) .gt. 0 ) then
	 irot = 22
	 goto 101
	 elseif ( abs(Y(1,2)) + abs(Y(3,4)) .le. 0 ) then
	 
	 elseif ( abs(Y(1,3)) + abs(Y(2,4)) .le. 0 ) then
	 irot = 22
	 goto 101
	 else
	 print*,'Error in CYPD0M0(): 2 zero mass, soft cannot be'
	 stop
	 endif
	 endif

* Ref0. Eq. (3.79)
	 a = Y(2,4)*Y(3,4) - m4*Y(2,3)
	 b = Y(1,3)*Y(2,4) + Y(1,2)*Y(3,4) - Y(1,4)*Y(2,3)
	 c = Y(1,2)*Y(1,3) - m1*Y(2,3) 
	 d = Y(2,3)
	 
* Condition: Y12*Y13 .ne. 0 
* d cannot be zero here
* See Ref0. Eq. (3.49) 
* ax^2 + bx + c + i*eps*d = 0
	 call CYPRoot2(a, b, c, DBLE(d), x, ix)

* check for Landau singularity
	 if ( abs(x(1) - x(2)) .le. 0 ) then 
	 print*,'Error in CYPD0M0(): 2 zero mass and detY = 0'
	 goto 109
	 endif

* See Ref0. Eq. (2.5) 
* mj^2 x^2 + Yij x + mi2 - i*eps = 0
	 call CYPRoot2(m4, Y(1,4), m1, -1D0, cr, icr)
	 do k = 1,2
	 cr(k) = -1/cr(k)
	 enddo

* Ref0. Eq. (3.80)
	 if ( abs(Y(1,2)) + abs(Y(3,4)) .le. 0 ) then
	 do k = 1,2
	 cx(1) = -x(k)
	 icx(1) = -ix(k)
	 cx(2) = cr(1)
	 icx(2) = icr(1)
	 itemp = icx(1)*DBLE(cx(2))
	 if ( itemp .eq. 0 ) itemp = icx(1)
	 D0 = D0 + (-1D0)**k*CYPLi2nx(cx, icx, itemp, 2)
	 
	 cx(2) = cr(2)
	 icx(2) = icr(2)
	 itemp = icx(1)*DBLE(cx(2))
	 if ( itemp .eq. 0 ) itemp = icx(1)
	 D0 = D0 + (-1D0)**k*CYPLi2nx(cx, icx, itemp, 2)

	 D0 = D0 + (-1D0)**k*CYPLOG(-x(k),-ix(k))**2/2
	 
	 cx(1) = Y(2,4)/Y(2,3)
	 icx(1) = -DBLE(Y(2,3))
	 if ( icx(1) .eq. 0 ) icx(1) = 1D0

	 cx(2) = Y(1,3)/m1
	 icx(2) = -DBLE(m1 - Y(1,3))
	 if ( icx(2) .eq. 0 ) icx(2) = 1D0
* if x(k) = 0, ie. c = 0, then 2 cases are possible:
* either this log-term vanishes or diverge (3pt Landau). 
* We set the divergent part to 0 if it occurs.
	 if (abs(x(k)) .gt. 0) then 
	 D0 = D0 + (-1D0)**k*CYPLOG(-x(k),-ix(k))*
     &	 (CYPLOG(cx(1),icx(1)) + CYPLOG(cx(2),icx(2)))
	 endif
	 enddo
	 
	 D0 = D0/a/(x(1) - x(2))

c	 print*,'finite D0 2M0 2Y0:',D0
	 return	 
	 endif

* Ref0. Eq. (3.78)
	 do k = 1,2
	 cx(1) = -x(k)
	 icx(1) = -ix(k)
	 cx(2) = Y(3,4)/Y(1,3)
	 icx(2) = DBLE(Y(3,4) - Y(1,3))
	 if ( icx(2) .eq. 0 ) icx(2) = 1D0
	 itemp = icx(1)*DBLE(cx(2))
	 if ( itemp .eq. 0 ) itemp = icx(1)
	 D0 = D0 - (-1D0)**k*CYPLi2nx(cx, icx, itemp, 2)
	 
	 cx(2) = cr(1)
	 icx(2) = icr(1)
	 itemp = icx(1)*DBLE(cx(2))
	 if ( itemp .eq. 0 ) itemp = icx(1)
	 D0 = D0 + (-1D0)**k*CYPLi2nx(cx, icx, itemp, 2)
	 
	 cx(2) = cr(2)
	 icx(2) = icr(2)
	 itemp = icx(1)*DBLE(cx(2))
	 if ( itemp .eq. 0 ) itemp = icx(1)
	 D0 = D0 + (-1D0)**k*CYPLi2nx(cx, icx, itemp, 2)

	 cx(2) = Y(2,4)/Y(1,2)
	 icx(2) = DBLE(Y(2,4) - Y(1,2))
	 if ( icx(2) .eq. 0 ) icx(2) = 1D0
	 itemp = icx(1)*DBLE(cx(2))
	 if ( itemp .eq. 0 ) itemp = icx(1)
	 D0 = D0 - (-1D0)**k*CYPLi2nx(cx, icx, itemp, 2)
	 
	 cx(1) = Y(1,2)/Y(2,3)
	 icx(1) = -DBLE(Y(2,3))
	 if ( icx(1) .eq. 0 ) icx(1) = 1D0

	 cx(2) = Y(1,3)/m1
	 icx(2) = -DBLE(m1 - Y(1,3))
	 if ( icx(2) .eq. 0 ) icx(2) = 1D0
* if x(k) = 0, ie. c = 0, then 2 cases are possible:
* either this log-term vanishes or diverge (3pt Landau). 
* We set the divergent part to 0 if it occurs.
	 if (abs(x(k)) .gt. 0) then 
	 D0 = D0 + (-1D0)**k*CYPLOG(-x(k),-ix(k))*
     &	 (CYPLOG(cx(1),icx(1)) + CYPLOG(cx(2),icx(2)))
	 endif
	 enddo
	 
	 D0 = D0/a/(x(1) - x(2))

c	 print*,'finite D0 2M0:',D0
	 return
	endif
***]2 zero masses

***[1 zero mass
	if ( abs(m3) .le. 0 ) then 
* possible rotations: 1,6,7,12,20,22.
* check for log(0) singularity

#ifdef DEBUG_CYPD0M0
	print*,'CYPD0M0: begin 1 zero mass'
#endif

	 if ( abs(Y(1,3)) .le. 0 ) then 
	 irot = 20
	 goto 101
	 endif

* See Ref0. Eq. (2.5) 
* mj^2 x^2 + Yij x + mi2 - i*eps = 0
	 call CYPRoot2(m4, Y(2,4), m2, -1D0, cr, icr)
	 do k = 1,2
	 cr(k) = -1/cr(k)
	 enddo

* choose r24 which cannot be zero
	 r24 = cr(1)
	 ir24 = icr(1)
	 if ( abs(log(abs(cr(2)))) .lt. abs(log(abs(cr(1)))) ) then
          r24 = cr(2)
     	  ir24 = icr(2)
	 endif

* Ref0. Eq. (3.73)
	 a = m2*r24*Y(3,4) - m4*Y(2,3)
	 b = Y(1,3)*(m2*r24 - m4/r24) + Y(1,2)*Y(3,4) - Y(1,4)*Y(2,3)
	 c = Y(1,3)*(Y(1,2) - Y(1,4)/r24) - m1*(Y(2,3) - Y(3,4)/r24) 
	 d = Y(2,3) - Y(3,4)/r24
	 if ( DBLE(d) .eq. 0 ) d = 1D0
* Condition: Y13 .ne. 0 
* See Ref0. Eq. (3.49) 
* ax^2 + bx + c + i*eps*d = 0

* possible rotations: 1,6,7,12,20,22.
	if ( abs(a) .le. 0 ) then 
	 if ( irot .eq. 1 ) then
          irot = 6
	 elseif ( irot .eq. 6 ) then
	  irot = 7
	 elseif ( irot .eq. 7 ) then
	  irot = 12
	 elseif ( irot .eq. 12 ) then
	  irot = 20
	 elseif ( irot .eq. 20 ) then
	  irot = 22
	 elseif ( irot .eq. 22 ) then
	  goto 109
	 endif
	 goto 101
	endif

	 call CYPRoot2(a, b, c, DBLE(d), x, ix)

* check for Landau singularity
	 if ( abs(x(1) - x(2)) .le. 0 ) then 
	 if ( runmode .ge. 2 ) then 
	 print*,'Error in CYPD0M0(): 1 zero mass and detY = 0'
	 endif
	 goto 109
	 endif

* See Ref0. Eq. (2.5) 
* mj^2 x^2 + Yij x + mi2 - i*eps = 0
	 call CYPRoot2(m4, Y(1,4), m1, -1D0, cr, icr)
	 do k = 1,2
	 cr(k) = -1/cr(k)
	 enddo

	 call CYPRoot2(m2, Y(1,2), m1, -1D0, dr, idr)
	 do k = 1,2
	 dr(k) = -1/dr(k)
	 enddo
	
* Ref0. Eq. (3.76)
	 do k = 1,2
	 cx(1) = -x(k)
	 icx(1) = -ix(k)
	 cx(2) = Y(3,4)/Y(1,3)
	 icx(2) = DBLE(Y(3,4) - Y(1,3))
	 if ( icx(2) .eq. 0 ) icx(2) = 1D0
	 itemp = icx(1)*DBLE(cx(2))
	 if ( itemp .eq. 0 ) itemp = icx(1)
	 D0 = D0 - (-1D0)**k*CYPLi2nx(cx, icx, itemp, 2)
	 
	 cx(2) = cr(1)
	 icx(2) = icr(1)
	 itemp = icx(1)*DBLE(cx(2))
	 if ( itemp .eq. 0 ) itemp = icx(1)
	 D0 = D0 + (-1D0)**k*CYPLi2nx(cx, icx, itemp, 2)
	 
	 cx(2) = cr(2)
	 icx(2) = icr(2)
	 itemp = icx(1)*DBLE(cx(2))
	 if ( itemp .eq. 0 ) itemp = icx(1)
	 D0 = D0 + (-1D0)**k*CYPLi2nx(cx, icx, itemp, 2)

	 cx(1) = -x(k)*r24
	 icx(1) = -ix(k)*DBLE(r24)
	 cx(2) = Y(2,3)/Y(1,3)
	 icx(2) = DBLE(Y(2,3) - Y(1,3))
	 if ( icx(2) .eq. 0 ) icx(2) = 1D0
	 itemp = icx(1)*DBLE(cx(2))
	 if ( itemp .eq. 0 ) itemp = icx(1)
	 D0 = D0 + (-1D0)**k*CYPLi2nx(cx, icx, itemp, 2)
	 
	 cx(2) = dr(1)
	 icx(2) = idr(1)
	 itemp = icx(1)*DBLE(cx(2))
	 if ( itemp .eq. 0 ) itemp = icx(1)
	 D0 = D0 - (-1D0)**k*CYPLi2nx(cx, icx, itemp, 2)
	 
	 cx(2) = dr(2)
	 icx(2) = idr(2)
	 itemp = icx(1)*DBLE(cx(2))
	 if ( itemp .eq. 0 ) itemp = icx(1)
	 D0 = D0 - (-1D0)**k*CYPLi2nx(cx, icx, itemp, 2)

	 cx(1) = -x(k)
	 if ( DIMAG(cx(1)) .eq. 0 ) cx(1) = -x(k) - cI*ix(k)
	 cx(2) = r24
	 if ( DIMAG(cx(2)) .eq. 0 ) cx(2) = r24 + cI*ir24
	 ctemp = -x(k)*r24
	 if ( DIMAG(ctemp) .eq. 0 ) ctemp = ctemp - cI*DBLE(r24)*ix(k)
c	 print*,'Eta ctemp =', cx(1), cx(2), ctemp
 	 ctemp = CYPetax(cx(1), cx(2), ctemp)
	 if ( abs(ctemp) .gt. 0 ) then 
	  
* Ref0. Eq. (3.75)
	 cx(1) = (Y(1,2) - Y(1,4)/r24 + (m2*r24 - m4/r24)*x(k))/
     &	 (Y(2,3) - Y(3,4)/r24)
	 icx(1) = ix(k)*DBLE((m2*r24 - m4/r24)/(Y(2,3) - Y(3,4)/r24))
	 if ( icx(1) .eq. 0 ) icx(1) = 1D0

	 cx(2) = Y(1,3)/m1
	 icx(2) = -DBLE(m1 - Y(1,3))
	 if ( icx(2) .eq. 0 ) icx(2) = 1D0

	 D0 = D0 - (-1D0)**k*ctemp*
     &	 (CYPLOG(cx(1),icx(1)) + CYPLOG(cx(2),icx(2)))
	 endif
	 enddo
	 
	 D0 = D0/a/(x(1) - x(2))

c	 print*,'finite D0, r24 1M0:', D0, r24
c	 return	
	
	 return
	endif
***]1 zero mass

* Ninh checking DD1
#ifdef D0_DD_MASS
***[Ref0. Eq. (3.61)

* See Ref0. Eq. (2.5) 
* mj^2 x^2 + Yij x + mi2 - i*eps = 0
* get r20 (r31 here)
* choose r31 which has to be real and positive here
* the ieps of r31 is needed
	r31 = 0
	call CYPRoot2(m3, Y(1,3), m1, -1D0, cr, icr)
	do k = 1,2
	icr(k) = icr(k)/abs(cr(k))**2
	cr(k) = -1/cr(k)
	enddo

c	 print*,'check r31 root2 = ', cr, icr ! Ninh check11

#ifdef CHECK_Root_Kam
* x31 cannot be zero here
	call CYPKam(p1p2, rm1, rm3, x31, ix31, ix31s)
	cr(1) = rm3/rm1*x31
	icr(1) = DBLE(rm3/rm1)*ix31
	cr(2) = rm3/rm1/x31
	icr(2) = -DBLE(rm3/rm1)*ix31/abs(x31)**2

c	 print*,'check r31 Kam = ', cr, icr
#endif

	if ( DIMAG(cr(1)) .eq. 0 .and. DBLE(cr(1)) .gt. 0 ) then 
	  r31 = cr(1)
	  ir31 = icr(1)
	endif

	if ( DIMAG(cr(2)) .eq. 0 .and. DBLE(cr(2)) .gt. 0 ) then 
	 if ( abs(r31) .le. 0 ) then
          r31 = cr(2)
	  ir31 = icr(2)
	 elseif ( abs(log(cr(2))) .lt. abs(log(cr(1))) ) then
          r31 = cr(2)
	  ir31 = icr(2)
	 endif
	endif

* the above steps are needed to get ir31
	if ( abs(p1p2) .le. 0 ) then
	 r31 = 1D0
	 ir31 = 0
	else
* Ninh checking
ccccccccccccccccc
C Ninh check11
	 if ( irot .eq. 1 ) then
          irot = 2
	 elseif ( irot .eq. 2 ) then
	  irot = 3
	 elseif ( irot .eq. 3 ) then
	  irot = 7
	 elseif ( irot .eq. 7 ) then
	  irot = 8
	 elseif ( irot .eq. 8 ) then
	  irot = 14
	 elseif ( irot .eq. 14 ) then
	  goto 109
	 endif
	 goto 101
ccccccccccccccccc
	 if ( runmode .ge. 1 ) then 
	 print*,'WARNING from CYPD0M0(): this case not OK'
	 print*,'It seems to be the ieps problem'
	 print*,'r31, ir31 =', r31, ir31
	 endif
	 goto 109
	endif

* possible rotations: 1,2,3,7,8,14
	if ( abs(r31) .le. 0 ) then 
	 if ( irot .eq. 1 ) then
          irot = 2
	 elseif ( irot .eq. 2 ) then
	  irot = 3
	 elseif ( irot .eq. 3 ) then
	  irot = 7
	 elseif ( irot .eq. 7 ) then
	  irot = 8
	 elseif ( irot .eq. 8 ) then
	  irot = 14
	 elseif ( irot .eq. 14 ) then
	  goto 109
	 endif
	 goto 101
	endif

* See Ref0. Eq. (2.5) 
* mj^2 x^2 + Yij x + mi2 - i*eps = 0
	 call CYPRoot2(m4, Y(2,4), m2, -1D0, cr, icr)
	 do k = 1,2
	 cr(k) = -1/cr(k)
	 enddo

#ifdef CHECK_Root_Kam
	call CYPKam(p2p3, rm2, rm4, x24, ix24, ix24s)
	cr(1) = rm4/rm2*x24
	icr(1) = DBLE(rm4/rm2)*ix24
	cr(2) = rm4/rm2/x24
	icr(2) = -DBLE(rm4/rm2)*ix24/abs(x24)**2
#endif

* choose r24 which cannot be zero
	 r24 = cr(1)
	 ir24 = icr(1)
	 if ( abs(log(abs(cr(2)))) .lt. abs(log(abs(cr(1)))) ) then
          r24 = cr(2)
     	  ir24 = icr(2)
	 endif

c	 print*,'r31, r24 =',r31, r24

* Ref0. Eq. (3.50)
	 a = m2*r24*(Y(3,4) - Y(1,4)/r31) - m4*(Y(2,3) - Y(1,2)/r31)
	 b = (m3*r31 - m1/r31)*(m2*r24 - m4/r24) + Y(1,2)*Y(3,4) 
     &    - Y(1,4)*Y(2,3)
	 c = m3*r31*(Y(1,2) - Y(1,4)/r24) - m1*(Y(2,3) - Y(3,4)/r24) 
	 d = Y(2,3) - Y(1,2)/r31 - Y(3,4)/r24 + Y(1,4)/(r31*r24)
	 if ( DBLE(d) .eq. 0 ) d = 1D0

* See Ref0. Eq. (3.49) 
* ax^2 + bx + c + i*eps*d = 0
* Ninh checking
	 if ( abs(a) .le. 0 ) then 

cccccccccccc
C Ninh check11
* possible rotations: 1,2,3,7,8,14
c	if ( abs(a) .le. 0 ) then 
	 if ( irot .eq. 1 ) then
          irot = 2
	 elseif ( irot .eq. 2 ) then
	  irot = 3
	 elseif ( irot .eq. 3 ) then
	  irot = 7
	 elseif ( irot .eq. 7 ) then
	  irot = 8
	 elseif ( irot .eq. 8 ) then
	  irot = 14
	 elseif ( irot .eq. 14 ) then
	  goto 109
	 endif
	 goto 101
c	endif
cccccccccccc


	 if ( runmode .ge. 1 ) then 
	 print*,'Error in CYPD0M0(): a = 0 (Ref0. Eq. (3.50))'
	 print*,'todo: try rotation'
	 endif
	 goto 109
	 endif

	 call CYPRoot2(a, b, c, DBLE(d), x, ix)
c 	 print*,'two roots =',x, ix

* check for Landau singularity
	 if ( abs(x(1) - x(2)) .le. 0 ) then 
	 if ( runmode .ge. 2 ) then 
	 print*,'Error in CYPD0M0(): All massive and detY = 0'
	 print*,'todo: try rotation'
	 endif
	 goto 109
	 endif

* See Ref0. Eq. (2.5) 
* solve for rij
* mj^2 x^2 + Yij x + mi2 - i*eps = 0
	 call CYPRoot2(m4, Y(1,4), m1, -1D0, cr, icr)
	 do k = 1,2
	 cr(k) = -1/cr(k)
	 enddo

#ifdef CHECK_Root_Kam
	call CYPKam(p4, rm1, rm4, x14, ix14, ix14s)
	cr(1) = rm4/rm1*x14
	icr(1) = DBLE(rm4/rm1)*ix14
	cr(2) = rm4/rm1/x14
	icr(2) = -DBLE(rm4/rm1)*ix14/abs(x14)**2
#endif

	 call CYPRoot2(m2, Y(1,2), m1, -1D0, dr, idr)
	 do k = 1,2
	 dr(k) = -1/dr(k)
	 enddo

#ifdef CHECK_Root_Kam
	call CYPKam(p1, rm1, rm2, x12, ix12, ix12s)
	dr(1) = rm2/rm1*x12
	idr(1) = DBLE(rm2/rm1)*ix12
	dr(2) = rm2/rm1/x12
	idr(2) = -DBLE(rm2/rm1)*ix12/abs(x12)**2
#endif

	 call CYPRoot2(m2, Y(3,2), m3, -1D0, er, ier)
	 do k = 1,2
	 er(k) = -1/er(k)
	 enddo

#ifdef CHECK_Root_Kam
	call CYPKam(p2, rm3, rm2, x32, ix32, ix32s)
	er(1) = rm2/rm3*x32
	ier(1) = DBLE(rm2/rm3)*ix32
	er(2) = rm2/rm3/x32
	ier(2) = -DBLE(rm2/rm3)*ix32/abs(x32)**2
#endif

	 call CYPRoot2(m4, Y(3,4), m3, -1D0, fr, ifr)
	 do k = 1,2
	 fr(k) = -1/fr(k)
	 enddo

#ifdef CHECK_Root_Kam
	call CYPKam(p3, rm3, rm4, x34, ix34, ix34s)
	fr(1) = rm4/rm3*x34
	ifr(1) = DBLE(rm4/rm3)*ix34
	fr(2) = rm4/rm3/x34
	ifr(2) = -DBLE(rm4/rm3)*ix34/abs(x34)**2
#endif

* Ref0. Eq. (3.61)

	 do k = 1,2

	 do l = 1,2
	 cx(1) = -x(k)
	 icx(1) = -ix(k)
	 cx(2) = fr(l)/r31
	 icx(2) = ifr(l)*DBLE(r31) - ir31*DBLE(fr(l))
	 if ( icx(2) .eq. 0 ) icx(2) = ifr(l)

	 itemp = icx(1)*DBLE(cx(2))
	 if ( itemp .eq. 0 ) itemp = icx(1)

c	 print*,'check cx 1 = ', cx, icx, itemp

	 D0 = D0 - (-1D0)**k*CYPLi2nx(cx, icx, itemp, 2)
	 
	 cx(2) = cr(l)
	 icx(2) = icr(l)

	 itemp = icx(1)*DBLE(cx(2))
	 if ( itemp .eq. 0 ) itemp = icx(1)

c	 print*,'check cx 2 = ', cx, icx, itemp

	 D0 = D0 + (-1D0)**k*CYPLi2nx(cx, icx, itemp, 2)

	 cx(1) = -x(k)*r24
	 icx(1) = -ix(k)*DBLE(r24)
	 cx(2) = er(l)/r31
	 icx(2) = ier(l)*DBLE(r31) - ir31*DBLE(er(l))
	 if ( icx(2) .eq. 0 ) icx(2) = ier(l)

	 itemp = icx(1)*DBLE(cx(2))
	 if ( itemp .eq. 0 ) itemp = icx(1)

c	 print*,'check cx 3 = ', cx, icx, itemp

	 D0 = D0 + (-1D0)**k*CYPLi2nx(cx, icx, itemp, 2)
	 
	 cx(2) = dr(l)
	 icx(2) = idr(l)

	 itemp = icx(1)*DBLE(cx(2))
	 if ( itemp .eq. 0 ) itemp = icx(1)

c	 print*,'check cx 4 = ', cx, icx, itemp

	 D0 = D0 - (-1D0)**k*CYPLi2nx(cx, icx, itemp, 2)
	 enddo

	 cx(1) = -x(k)
	 if ( DIMAG(cx(1)) .eq. 0 ) cx(1) = -x(k) - cI*ix(k)
	 cx(2) = r24
	 if ( DIMAG(cx(2)) .eq. 0 ) cx(2) = r24 + cI*ir24
	 ctemp = -x(k)*r24
	 if ( DIMAG(ctemp) .eq. 0 ) ctemp = ctemp - cI*DBLE(r24)*ix(k)

c	 print*,'Eta ctemp =', cx(1), cx(2), ctemp

 	 ctemp = CYPetax(cx(1), cx(2), ctemp)

	 if ( abs(ctemp) .gt. 0 ) then 
	  
* Ref0. Eq. (3.60)
	 cx(1) = (r31*(Y(2,3) - Y(3,4)/r24) + (m2*r24 - m4/r24)*x(k))/
     &	 (Y(1,2) - Y(1,4)/r24 + (m2*r24 - m4/r24)*x(k))/r31**2
	 icx(1) = -ix(k)*DBLE(m2*r24 - m4/r24)*DBLE(
     &    r31*(Y(2,3) - Y(3,4)/r24) - (Y(1,2) - Y(1,4)/r24))
	 if ( icx(1) .eq. 0 ) icx(1) = 1D0

	 cx(2) = m3/m1
	 icx(2) = -DBLE(m1 - m3)
	 if ( icx(2) .eq. 0 ) icx(2) = 1D0

	 D0 = D0 + (-1D0)**k*ctemp*
     &	 (CYPLOG(cx(1),icx(1)) - CYPLOG(cx(2),icx(2)))
	 endif
	 enddo
	 
	 D0 = D0/a/(x(1) - x(2))

c	 print*,'finite 4M D0:', D0 ! NInh check11
	 return
***]Ref0. Eq. (3.61)
#endif

 109    continue
	fail = 2
	return
	end

#ifdef D0_DD_MASS
************************************************************************
	subroutine CYPD0DD(fail,q1,q2,q3,q4,q1q2,q2q3,
     >   n1,n2,n3,n4,D0)
************************************************************************
* Calculate the finite scalar 4-point integrals for non-zero masses and 
* non-zero momenta.
* Ref0. Denner and Dittmaier arXiv: 1005.2076
* Use Ref0. Eq. (3.65) for gamma NOT very close to 0 or 1.
* Author: Le Duc Ninh (2011)
************************************************************************
	implicit none
        CVAR q1, q2, q3, q4, q1q2, q2q3
        CVAR n1, n2, n3, n4
        CVAR p1, p2, p3, p4, p1p2, p2p3
        CVAR m1, m2, m3, m4
        CVAR alpha, alpha1, alpha2
	CVAR rb, rg, rj
	CVAR D0, ctemp, mm, pp
	integer fail, irot

#include "lt.h"

        pmin = min(abs(q1),abs(q2),abs(q3),abs(q4),abs(q1q2),abs(q2q3))
	if ( pmin .le. 0 ) then
         print*,'Error in CYPD0DD(): zero momenta, invalid input'
	 print*,'check rotation'
	 stop
	endif

	D0 = 0
	fail = 0

	irot = 1
 201    continue

	p1 = q1
	p2 = q2
	p3 = q3
	p4 = q4
	p1p2 = q1q2
	p2p3 = q2q3
	m1 = n1
	m2 = n2
	m3 = n3
	m4 = n4

* WARNING: pi and mi change after this call.
	 call CYPRotate24(p1,p2,p3,p4,p1p2,p2p3,
     >   m1,m2,m3,m4,irot)

* Ref0. Eq. (3.64)
***{ calculate alpha
* alpha is one root of: rg*alpha^2+rj*alpha+rb=0
* rg cannot be zero here.
* we do not need the ieps for alpha.

      rb = p4
      rg = p2p3
      rj = p1 - p2p3 - p4

      alpha2 = rj**2-4*rg*rb
      if ( DBLE(alpha2) .lt. 0 ) then 
        print*,'Error in CYPD0DD(): complex roots, invalid input'
	print*,'check rotation'
	stop
      endif

      alpha1 = (-rj-sqrt(alpha2))/2/rg
      alpha2 = (-rj+sqrt(alpha2))/2/rg
      
      if (abs(alpha1).gt.abs(alpha2)) then 
        alpha2 = alpha1
	alpha1 = rb/rg/alpha1
      else
        alpha1 = rb/rg/alpha2
      endif
     		
*** alpha(i) cannot be 0 or 1 here
* Ref0. Eq. (3.67)
      if ( DBLE(alpha1) .gt. 0 .and. DBLE(alpha1) .lt. 1 ) then
       alpha = alpha1
      elseif ( DBLE(alpha2) .gt. 0 .and. DBLE(alpha2) .lt. 1 ) then
       alpha = alpha2
      else
       alpha = 0
      endif 	            	  
***} calculate alpha

c      print*,'irot =',irot
c      print*,'alpha i=',alpha1, alpha2, alpha

* possible rotations: 1, 7, 20
        if ( abs(alpha) .le. 0 ) then
	 if ( irot .eq. 1 ) then
	  irot = 7
	 elseif ( irot .eq. 7 ) then
	  irot = 20
         elseif ( irot .eq. 20 ) then 
          print*,'Error in CYPD0DD(): invalid real roots'
	  print*,'check rotation'
	  stop
         endif
	 goto 201
	endif

* WARNING: if alpha is very close to 0 or 1 then 
* numerical cancellation occurs.
* Here 1D-8 is just a safe cut-off between Ref0. Eq. (3.65) 
* and D0C mthod.
        if ( abs(alpha) .lt. 1D-8 .or. 
     &  abs(1-alpha) .lt. 1D-8 ) then
	 fail = 2
	 goto 209
	endif

c       print*,'Info CYPD0DD() alpha =', alpha

* Ref0. Eq. (3.63)
	mm = alpha*(m2 - p1) + (1 - alpha)*(m4 - p4)
* Ref0. Eq. (3.66)
        pp = alpha*(p2 - p1) + (1 - alpha)*(p3 - p4)

* Ref0. Eq. (3.65)
	call CYPD0M0(fail,p4,alpha**2*p2p3,pp,p1p2,czero,p3,
     &  m1,m4,mm,m3,D0)
	ctemp = D0*alpha

	call CYPD0M0(fail,p1,(1-alpha)**2*p2p3,pp,p1p2,czero,p2,
     &  m1,m2,mm,m3,D0)
	D0 = ctemp + (1-alpha)*D0

c        print*,'End CYPD0DD() =', D0
 209    continue
	return
	end
#endif
***********************************************************************
      subroutine CYPRotate(rp1,rp2,rp3,rp4,rp1p2,rp2p3,
     > cm1,cm2,cm3,cm4,irot)
     
      IMPLICIT NONE
      
      CVAR rp1,rp2,rp3,rp4,rp1p2,rp2p3
      CVAR cm1,cm2,cm3,cm4
      CVAR CYPLambda, temp
      integer irot
#include "lt.h"      

***{rotation

      irot = 0 

***1{detect special cases of 2 lightlike momenta
***{D0C for boxes with at least 2 zero external momenta
       if (abs(rp1).le.0 .and. abs(rp3).le.0) then
	irot = 1
	return

       else if (abs(rp2).le.0 .and. abs(rp4).le.0) then
        irot = 6
        return
       
       
       else if (abs(rp1p2).le.0 .and. abs(rp2p3).le.0) then
	irot = 3
        return

***1}detect special cases of 2 lightlike momenta

***2{detect special cases of 1 lightlike momenta
      elseif (abs(rp1).le.0) then
       irot = 1
       return
      
      else if (abs(rp2).le.0) then
	irot = 9
	return
      
      elseif (abs(rp3).le.0) then
       irot = 17
       return

      else if (abs(rp4).le.0) then
       irot = 5
       return 
       
      else if (abs(rp1p2).le.0) then
       irot = 3
       return 
       
      else if (abs(rp2p3).le.0) then
       irot = 11
       return 

***2}detect special cases of 1 lightlike momenta      
     
***3{detect special cases of lambda(pi,pj,pk)>=0

*** the last case ***
      elseif (abs(CYPLambda(rp1p2,rp1,rp2)).ge.0D0) then 
       irot = 1
       return 

      elseif (abs(CYPLambda(rp2p3,rp2,rp3)).ge.0D0) then
	irot = 10
	return
      
      elseif (abs(CYPLambda(rp1p2,rp3,rp4)).ge.0D0) then             
	irot = 17
	return
	
      elseif (abs(CYPLambda(rp2p3,rp4,rp1)).ge.0D0) then       	
	irot = 19
	return
      
***3}detect special cases of lambda(pi,pj,pk)>=0

***{the case of all lambda(pi,pj,pk)<0
      else
        irot = 25 
        return	
***}the case of all lambda(pi,pj,pk)<0

      endif
***}rotation

      end
************************************************************************
#ifdef D0_DD_MASS
      subroutine CYPRotateP5(rp1,rp2,rp3,rp4,rp1p2,rp2p3,
     > cm1,cm2,cm3,cm4,irot)
************************************************************************
* Rotate to the position P5 = p1p2 = 0, 
* or to the position where lambda(p2p3,p1,p4) > 0. 
* Author: Le Duc Ninh (2011)
************************************************************************
      IMPLICIT NONE
      
      CVAR rp1,rp2,rp3,rp4,rp1p2,rp2p3
      CVAR cm1,cm2,cm3,cm4
      CVAR CYPLambda, temp
      integer irot
#include "lt.h"      

***{rotation

      irot = 0 

***{detect special cases of 1 lightlike momenta
      if (abs(rp1).le.0) then
       irot = 3
       return
      
      else if (abs(rp2).le.0) then
	irot = 7
	return
      
      elseif (abs(rp3).le.0) then
       irot = 14
       return

      else if (abs(rp4).le.0) then
       irot = 2
       return 
       
      else if (abs(rp1p2).le.0) then
       irot = 1
       return 
       
      else if (abs(rp2p3).le.0) then
       irot = 8
       return 

***}detect special cases of 1 lightlike momenta      
     
***{detect special cases of lambda(pi,pj,pk)>=0

*** the last case ***
      elseif (abs(CYPLambda(rp1p2,rp1,rp2)).ge.0D0) then 
       irot = 8
       return 

      elseif (abs(CYPLambda(rp2p3,rp2,rp3)).ge.0D0) then
	irot = 15
	return
      
      elseif (abs(CYPLambda(rp1p2,rp3,rp4)).ge.0D0) then             
	irot = 19
	return
	
      elseif (abs(CYPLambda(rp2p3,rp4,rp1)).ge.0D0) then       	
	irot = 1
	return
      
***}detect special cases of lambda(pi,pj,pk)>=0

***{the case of all lambda(pi,pj,pk)<0
      else
        irot = 25 
        return	
***}the case of all lambda(pi,pj,pk)<0

      endif
***}rotation
      end
#endif
************************************************************************      
      subroutine CYPD0_abc(ra,rb,rg,rc,rh,rj,cd,ce,ck,cfx,signf,csum)
************************************************************************     
* Def. D0 = \int_0^1dx\int_0^xdy\int_0^ydz
* *\fr{1}{(ra*x^2+rb*y^2+rg*z^2+rc*xy+rh*xz+rj*yz+cd*x+ce*y+ck*z+cfx+i*signf)} 
* with signf=-eps, {ra,rb,rg,rc,rh,rj} are real, {cd,ce,ck,cfx} are complex.
* output: d0 = csum
* Method: barycentric technique
* important: variables "signX" is the sign of img(X) in case X becomes real.
* Author: Le Duc Ninh (2010)
************************************************************************
      IMPLICIT NONE

      CVAR csum,ctemp
      CVAR cd,ce,ck,cfx
      CVAR cdd,cee,cff,cjj,cii
      RVAR ra,rb,rg,rc,rh,rj,signf
      RVAR raa,rbb,rcc,rgg,rhh,signj
      RVAR alpha1,alpha2,alpha
      
#include "lt.h"

c      print*,'d0c_abc.F: hello'
c      print*,'ra,rb,rg,rc,rh,rj=',ra,rb,rg,rc,rh,rj

***{ the case of 2 opposite lightlike external momenta
      if (abs(ra).le.0 .and. abs(rg).le.0) then 
      raa=rb+rj
      rcc=rc+rh
      rgg=rb
      rhh=rc
      
      cdd=cd
      cee=ce+ck
      cff=cfx
      cii=ce
      
      call CYPfT13(raa,rcc,rgg,rhh,cdd,cee,cff,cii,
     > signf,csum)
     
      return 
      endif      
***} the case of 2 opposite lightlike external momenta      
     
***{ the case of 1 lightlike external momentum
      if (abs(rg).le.0) then
      
* the first term I1
      raa=ra
      rbb=rb
      rcc=rc
      rgg=rh
      rhh=rj
      
      cdd=cd
      cee=ce
      cff=cfx
      cjj=ck
      
      signj=signf  
      
      call CYPfT(raa,rbb,rcc,rgg,rhh,cdd,cee,cff,cjj,
     > signf,signj,csum)
     
* the second term I2      
      raa=ra
      rbb=rj+rb
      rcc=rc+rh
      rgg=rh
      rhh=rj
      
      cdd=cd
      cee=ce+ck
      cff=cfx
      cjj=ck
      
      signj=signf  
      
      call CYPfT(raa,rbb,rcc,rgg,rhh,cdd,cee,cff,cjj,
     > signf,signj,ctemp)
     
* final result: (I1-I2)     

      csum=csum-ctemp
      return     
      endif 
***} the case of 1 lightlike external momentum

***{ calculate alpha
* alpha is one root of: rg*alpha^2+rj*alpha+rb=0
* rg cannot be zero here.
* we do not need the ieps for alpha.

      if ( abs(rb) .le. 0 ) then 
       alpha = 0
       goto 401
      endif     	

      alpha2 = rj**2-4*rg*rb
      alpha1 = (-rj-sqrt(alpha2))/2/rg
      alpha2 = (-rj+sqrt(alpha2))/2/rg
      
      if (abs(alpha1).gt.abs(alpha2)) then 
        alpha2 = alpha1
	alpha1 = rb/rg/alpha1
      else
        alpha1 = rb/rg/alpha2
      endif
      
cc      print*,'alpha i=',alpha1, alpha2
      		
*** alpha(i) cannot be 0 or 1 here ***
    	
*** which one for alpha ? ***
      if (abs(LOG(abs(alpha1))).lt.abs(LOG(abs(alpha2)))) then
        alpha = alpha1
      else
        alpha = alpha2
      endif 	            	  
***} calculate alpha

 401  continue  

c      print*,'alpha =',alpha  
c      print*,'alpha i',alpha1,alpha2

* the first term: j=1
      raa=ra+rb+rc
      rbb=rg
      rcc=rj+rh
      rgg=-(rc+2D0*rb+(rh+rj)*alpha)
      rhh=-(rj+2D0*alpha*rg)
      
      cdd=cd+ce
      cee=ck
      cff=cfx
      cjj=-(ce+ck*alpha)
      
      signj=-signf
      
      call CYPfT(raa,rbb,rcc,rgg,rhh,cdd,cee,cff,cjj,
     > signf,signj,csum)
     
c      print*,'d0c_abc.F csum1 =',csum 

* the second term: j=2
      if (alpha.eq.1D0) then 
      ctemp = 0
      
      else 
      raa=ra
      rbb=rg+rj+rb
      rcc=rc+rh
      rgg=(rc+rh*alpha)/(1D0-alpha)
      rhh=(rj+2D0*alpha*rg)
      
      cdd=cd
      cee=ce+ck
      cff=cfx
      cjj=(ce+ck*alpha)/(1D0-alpha)
      
      signj=signf/(1D0-alpha)
      
      call CYPfT(raa,rbb,rcc,rgg,rhh,cdd,cee,cff,cjj,
     > signf,signj,ctemp)
     
c      print*,'d0c_abc.F csum2 =',ctemp      
     
      endif
     
      csum=csum+ctemp
      
* the third term: j=3
      if (alpha.eq.0) then
      ctemp=0
      
      else 
      raa=ra
      rbb=rb
      rcc=rc
      rgg=rc/alpha+rh
      rhh=-(rj+2D0*alpha*rg)
      
      cdd=cd
      cee=ce
      cff=cfx
      cjj=ce/alpha+ck
      
      signj=signf/alpha
      
      call CYPfT(raa,rbb,rcc,rgg,rhh,cdd,cee,cff,cjj,
     > signf,signj,ctemp)
      endif
      
c      print*,'d0c_abc.F csum3 =',ctemp            
     
      csum=csum+ctemp

      end

************************************************************************
	subroutine CYPDtensor(para, D)
************************************************************************
* Method: we use Passarino-Veltman reduction method for DetGram .ne. 0.  
* This subroutine is taken and adapted from LoopTools of Thomas Hahn.  
************************************************************************
	implicit none
	integer i
	CVAR para(Pdd)
	CVAR p1, p2, p3, p4, p1p2, p2p3, m1, m2, m3, m4
	CVAR f1, f2, f3
	CVAR G(3,3),detG,detG2,alpha(0:3)
	CVAR D(Ndd)

	CVAR c0sum, c1sum, c2sum, csum
	CVAR c00sum, c11sum, c12sum, c22sum
	CVAR in(3)
	integer C234, C134, C124, C123, Di(0:3)

!	CVAR aa

	integer XYPCget
	external XYPCget

#include "lt.h"

	integer perm(3)
#define IN(i) in(perm(i))
#define SOLVE_SETUP CYPDecomp(3, G,3, perm)
#define SOLVE(b) CYPSolve(3, G,3, b)

	do i = 2, Ndd
	 D(i) = 0
	enddo

c	print*,'CYPDtensor hello'

	p1 = para(1)
	p2 = para(2)
	p3 = para(3)
	p4 = para(4)
	p1p2 = para(5)
	p2p3 = para(6)
	m1 = para(7)
	m2 = para(8)
	m3 = para(9)
	m4 = para(10)

	G(1,1) = 2*p1
	G(2,2) = 2*p1p2
	G(3,3) = 2*p4
	G(1,2) = p1 + p1p2 - p2
	G(2,1) = G(1,2)
	G(1,3) = p1 + p4 - p2p3
	G(3,1) = G(1,3)
	G(2,3) = p1p2 - p3 + p4
	G(3,2) = G(2,3)

	detG=G(1,1)*G(2,2)*G(3,3)+2*G(1,2)*G(2,3)*G(1,3)
     >  -G(1,3)*G(2,2)*G(1,3)-G(1,1)*G(2,3)*G(2,3)
     >  -G(3,3)*G(1,2)*G(1,2)

        if ( abs(detG) .eq. 0d0 ) then
         return
        endif

***[using Passarino-Veltman method
	C234 = XYPCget(p2, p3, p2p3, m2, m3, m4)
	C134 = XYPCget(p1p2, p3, p4, m1, m3, m4)
	C124 = XYPCget(p1, p2p3, p4, m1, m2, m4)
	C123 = XYPCget(p1, p2, p1p2, m1, m2, m3)

	
	c0sum = CvalYP(cc0,C234) + CvalYP(cc1,C234) + CvalYP(cc2,C234)
	c1sum = CvalYP(cc1,C234) + CvalYP(cc11,C234) + CvalYP(cc12,C234)
	c2sum = CvalYP(cc2,C234) + CvalYP(cc12,C234) + CvalYP(cc22,C234)
	csum = c0sum + c1sum + c2sum

	c00sum = CvalYP(cc00,C234) +
     &    CvalYP(cc001,C234) + CvalYP(cc002,C234)
	c11sum = CvalYP(cc11,C234) +
     &    CvalYP(cc111,C234) + CvalYP(cc112,C234)
	c12sum = CvalYP(cc12,C234) +
     &    CvalYP(cc112,C234) + CvalYP(cc122,C234)
	c22sum = CvalYP(cc22,C234) +
     &    CvalYP(cc122,C234) + CvalYP(cc222,C234)

	f1 = m2 - m1 - p1
	f2 = m3 - m1 - p1p2
	f3 = m4 - m1 - p4

	call SOLVE_SETUP

	IN(1) = f1*D(dd0) - CvalYP(cc0,C234) + CvalYP(cc0,C134)
	IN(2) = f2*D(dd0) - CvalYP(cc0,C234) + CvalYP(cc0,C124)
	IN(3) = f3*D(dd0) - CvalYP(cc0,C234) + CvalYP(cc0,C123)

	call SOLVE(in)
	D(dd1) = in(1)
	D(dd2) = in(2)
	D(dd3) = in(3)

	D(dd00) = m1*D(dd0) - .5D0*
     &    (D(dd1)*f1 + D(dd2)*f2 + D(dd3)*f3 - CvalYP(cc0,C234))

	IN(1) = f1*D(dd1) + c0sum - 2*D(dd00)
	IN(2) = f2*D(dd1) + c0sum + CvalYP(cc1,C124)
	IN(3) = f3*D(dd1) + c0sum + CvalYP(cc1,C123)
	call SOLVE(in)
	D(dd11) = in(1)
	D(dd12) = in(2)
	D(dd13) = in(3)

	IN(1) = f1*D(dd2) - CvalYP(cc1,C234) + CvalYP(cc1,C134)
	IN(2) = f2*D(dd2) - CvalYP(cc1,C234) - 2*D(dd00)
	IN(3) = f3*D(dd2) - CvalYP(cc1,C234) + CvalYP(cc2,C123)
	call SOLVE(in)
	D(dd12) = .5D0*(D(dd12) + in(1))
	D(dd22) = in(2)
	D(dd23) = in(3)

	IN(1) = f1*D(dd3) - CvalYP(cc2,C234) + CvalYP(cc2,C134)
	IN(2) = f2*D(dd3) - CvalYP(cc2,C234) + CvalYP(cc2,C124)
	IN(3) = f3*D(dd3) - CvalYP(cc2,C234) - 2*D(dd00)
	call SOLVE(in)
	D(dd13) = .5D0*(D(dd13) + in(1))
	D(dd23) = .5D0*(D(dd23) + in(2))
	D(dd33) = in(3)

c	print*,'Dcoeffs.F 2, D33=',D(dd33)

***[PV rank3
	if ( rankmax4 .ge. 3 ) then

c	IN(1) = f1*D(dd00) - CvalYP(cc00,C234) + CvalYP(cc00,C134)
c	IN(2) = f2*D(dd00) - CvalYP(cc00,C234) + CvalYP(cc00,C124)
c	IN(3) = f3*D(dd00) - CvalYP(cc00,C234) + CvalYP(cc00,C123)
c	call SOLVE(in)
c	D(dd001) = in(1)
c	D(dd002) = in(2)
c	D(dd003) = in(3)
***{Ninh change
        D(dd001)=-CvalYP(cc0,C234) - CvalYP(cc1,C234) -CvalYP(cc2,C234)
        D(dd001)=D(dd001) + 2*m1*D(dd1) - f1*D(dd11) 
     >  - f2*D(dd12) - f3*D(dd13)
        D(dd001)=D(dd001)/4D0  
        
        D(dd002)=CvalYP(cc1,C234)
        D(dd002)=D(dd002) + 2*m1*D(dd2) - f1*D(dd12) 
     >  - f2*D(dd22) - f3*D(dd23)
        D(dd002)=D(dd002)/4D0  
        
        D(dd003)=CvalYP(cc2,C234)
        D(dd003)=D(dd003) + 2*m1*D(dd3) - f1*D(dd13) 
     >  - f2*D(dd23) - f3*D(dd33)
        D(dd003)=D(dd003)/4D0          

c        print*,'dd00i check',D(dd003),in(3)
***}Ninh change	
	

	IN(1) = f1*D(dd11) - csum - 4*D(dd001)
	IN(2) = f2*D(dd11) - csum + CvalYP(cc11,C124)
	IN(3) = f3*D(dd11) - csum + CvalYP(cc11,C123)
	call SOLVE(in)
	D(dd111) = in(1)
	D(dd112) = in(2)
	D(dd113) = in(3)

	IN(1) = f1*D(dd22) - CvalYP(cc11,C234) + CvalYP(cc11,C134)
	IN(2) = f2*D(dd22) - CvalYP(cc11,C234) - 4*D(dd002)
	IN(3) = f3*D(dd22) - CvalYP(cc11,C234) + CvalYP(cc22,C123)
	call SOLVE(in)
	D(dd122) = in(1)
	D(dd222) = in(2)
	D(dd223) = in(3)

	IN(1) = f1*D(dd33) - CvalYP(cc22,C234) + CvalYP(cc22,C134)
	IN(2) = f2*D(dd33) - CvalYP(cc22,C234) + CvalYP(cc22,C124)
	IN(3) = f3*D(dd33) - CvalYP(cc22,C234) - 4*D(dd003)
	call SOLVE(in)
	D(dd133) = in(1)
	D(dd233) = in(2)
	D(dd333) = in(3)

	IN(1) = f1*D(dd13) + c2sum - 2*D(dd003)
	IN(2) = f2*D(dd13) + c2sum + CvalYP(cc12,C124)
	IN(3) = f3*D(dd13) + c2sum - 2*D(dd001)
	call SOLVE(in)
	D(dd113) = .5D0*(D(dd113) + in(1))
	D(dd123) = in(2)
	D(dd133) = .5D0*(D(dd133) + in(3))

        endif
***]PV rank3

***[PV rank4
	if ( rankmax4 .ge. 4 ) then

	D(dd0000) = 1/3D0*(m1*D(dd00) -
     &    .5D0*(f1*D(dd001) + f2*D(dd002) + f3*D(dd003) -
     &      CvalYP(cc00,C234) - 1/6D0*xrat))
	D(dd0011) = 1/3D0*(m1*D(dd11) -
     &    .5D0*(f1*D(dd111) + f2*D(dd112) + f3*D(dd113) - csum))
	D(dd0012) = 1/3D0*(m1*D(dd12) -
     &    .5D0*(f1*D(dd112) + f2*D(dd122) + f3*D(dd123) + c1sum))
	D(dd0013) = 1/3D0*(m1*D(dd13) -
     &    .5D0*(f1*D(dd113) + f2*D(dd123) + f3*D(dd133) + c2sum))
	D(dd0022) = 1/3D0*(m1*D(dd22) -
     &    .5D0*(f1*D(dd122) + f2*D(dd222) + f3*D(dd223) -
     &      CvalYP(cc11,C234)))
	D(dd0023) = 1/3D0*(m1*D(dd23) -
     &    .5D0*(f1*D(dd123) + f2*D(dd223) + f3*D(dd233) -
     &      CvalYP(cc12,C234)))
	D(dd0033) = 1/3D0*(m1*D(dd33) -
     &    .5D0*(f1*D(dd133) + f2*D(dd233) + f3*D(dd333) -
     &      CvalYP(cc22,C234)))

	c1sum = c1sum + c11sum + c12sum
	c2sum = c2sum + c12sum + c22sum
	csum = csum + c1sum + c2sum

	IN(1) = f1*D(dd111) + csum - 6*D(dd0011)
	IN(2) = f2*D(dd111) + csum + CvalYP(cc111,C124)
	IN(3) = f3*D(dd111) + csum + CvalYP(cc111,C123)
	call SOLVE(in)
	D(dd1111) = in(1)
	D(dd1112) = in(2)
	D(dd1113) = in(3)

	IN(1) = f1*D(dd113) - c2sum - 4*D(dd0013)
	IN(2) = f2*D(dd113) - c2sum + CvalYP(cc112,C124)
	IN(3) = f3*D(dd113) - c2sum - 2*D(dd0011)
	call SOLVE(in)
	D(dd1113) = .5D0*(D(dd1113) + in(1))
	D(dd1123) = in(2)
	D(dd1133) = in(3)

	IN(1) = f1*D(dd122) + c11sum - 2*D(dd0022)
	IN(2) = f2*D(dd122) + c11sum - 4*D(dd0012)
	IN(3) = f3*D(dd122) + c11sum + CvalYP(cc122,C123)
	call SOLVE(in)
	D(dd1122) = in(1)
	D(dd1222) = in(2)
	D(dd1223) = in(3)

	IN(1) = f1*D(dd222) - CvalYP(cc111,C234) + CvalYP(cc111,C134)
	IN(2) = f2*D(dd222) - CvalYP(cc111,C234) - 6*D(dd0022)
	IN(3) = f3*D(dd222) - CvalYP(cc111,C234) + CvalYP(cc222,C123)
	call SOLVE(in)
	D(dd1222) = .5D0*(D(dd1222) + in(1))
	D(dd2222) = in(2)
	D(dd2223) = in(3)

	IN(1) = f1*D(dd233) - CvalYP(cc122,C234) + CvalYP(cc122,C134)
	IN(2) = f2*D(dd233) - CvalYP(cc122,C234) - 2*D(dd0033)
	IN(3) = f3*D(dd233) - CvalYP(cc122,C234) - 4*D(dd0023)
	call SOLVE(in)
	D(dd1233) = in(1)
	D(dd2233) = in(2)
	D(dd2333) = in(3)

	IN(1) = f1*D(dd333) - CvalYP(cc222,C234) + CvalYP(cc222,C134)
	IN(2) = f2*D(dd333) - CvalYP(cc222,C234) + CvalYP(cc222,C124)
	IN(3) = f3*D(dd333) - CvalYP(cc222,C234) - 6*D(dd0033)
	call SOLVE(in)
	D(dd1333) = in(1)
	D(dd2333) = .5D0*(D(dd2333) + in(2))
	D(dd3333) = in(3)
	
c	print*,'Dcoeffs.F 3, D3333=',D(dd3333)
	endif
***]PV rank4

***[PV rank5
	if ( rankmax4 .ge. 5 ) then
	
	c00sum = c00sum +
     &  CvalYP(cc001,C234) + CvalYP(cc0011,C234) + CvalYP(cc0012,C234) +
     &  CvalYP(cc002,C234) + CvalYP(cc0012,C234) + CvalYP(cc0022,C234)
	c11sum = c11sum +
     &  CvalYP(cc111,C234) + CvalYP(cc1111,C234) + CvalYP(cc1112,C234) +
     &  CvalYP(cc112,C234) + CvalYP(cc1112,C234) + CvalYP(cc1122,C234)
	c12sum = c12sum +
     &  CvalYP(cc112,C234) + CvalYP(cc1112,C234) + CvalYP(cc1122,C234) +
     &  CvalYP(cc122,C234) + CvalYP(cc1122,C234) + CvalYP(cc1222,C234)
	c22sum = c22sum +
     &  CvalYP(cc122,C234) + CvalYP(cc1122,C234) + CvalYP(cc1222,C234) +
     &  CvalYP(cc222,C234) + CvalYP(cc1222,C234) + CvalYP(cc2222,C234)
	c1sum = c1sum + c11sum + c12sum
	c2sum = c2sum + c12sum + c22sum
	csum = csum + c1sum + c2sum

c	IN(1) = f1*D(dd0000) - CvalYP(cc0000,C234) + CvalYP(cc0000,C134)
c	IN(2) = f2*D(dd0000) - CvalYP(cc0000,C234) + CvalYP(cc0000,C124)
c	IN(3) = f3*D(dd0000) - CvalYP(cc0000,C234) + CvalYP(cc0000,C123)
c	call SOLVE(in)
c	D(dd00001) = in(1)
c	D(dd00002) = in(2)
c	D(dd00003) = in(3)

***{Ninh change
        D(dd00002)=CvalYP(cc001,C234)
        D(dd00002)=D(dd00002) + 2*m1*D(dd002) - f1*D(dd0012) 
     >  - f2*D(dd0022) - f3*D(dd0023) -1/24D0*xrat
        D(dd00002)=D(dd00002)/8D0  
        
        D(dd00003)=CvalYP(cc002,C234)
        D(dd00003)=D(dd00003) + 2*m1*D(dd003) - f1*D(dd0013) 
     >  - f2*D(dd0023) - f3*D(dd0033) -1/24D0*xrat
        D(dd00003)=D(dd00003)/8D0          
        
        D(dd00001)=-CvalYP(cc00,C234)-CvalYP(cc001,C234)
     >	-CvalYP(cc002,C234)
        D(dd00001)=D(dd00001) + 2*m1*D(dd001) - f1*D(dd0011) 
     >  - f2*D(dd0012) - f3*D(dd0013) -1/24D0*xrat
        D(dd00001)=D(dd00001)/8D0          

c        print*,'dd00i check',D(dd00001),in(1)
***}Ninh change	

c	IN(1) = f1*D(dd0011) - c00sum - 4*D(dd00001)
c	IN(2) = f2*D(dd0011) - c00sum + CvalYP(cc0011,C124)
c	IN(3) = f3*D(dd0011) - c00sum + CvalYP(cc0011,C123)
c	call SOLVE(in)
c	D(dd00111) = in(1)
c	D(dd00112) = in(2)
c	D(dd00113) = in(3)

c	IN(1) = f1*D(dd0022) - CvalYP(cc0011,C234) + CvalYP(cc0011,C134)
c	IN(2) = f2*D(dd0022) - CvalYP(cc0011,C234) - 4*D(dd00002)
c	IN(3) = f3*D(dd0022) - CvalYP(cc0011,C234) + CvalYP(cc0022,C123)
c	call SOLVE(in)
c	D(dd00122) = in(1)
c	D(dd00222) = in(2)
c	D(dd00223) = in(3)

c	IN(1) = f1*D(dd0033) - CvalYP(cc0022,C234) + CvalYP(cc0022,C134)
c	IN(2) = f2*D(dd0033) - CvalYP(cc0022,C234) + CvalYP(cc0022,C124)
c	IN(3) = f3*D(dd0033) - CvalYP(cc0022,C234) - 4*D(dd00003)
c	call SOLVE(in)
c	D(dd00133) = in(1)
c	D(dd00233) = in(2)
c	D(dd00333) = in(3)
	
c	IN(1) = f1*D(dd0023) - CvalYP(cc0012,C234) + CvalYP(cc0012,C134)
c	IN(2) = f2*D(dd0023) - CvalYP(cc0012,C234) - 2*D(dd00003)
c	IN(3) = f3*D(dd0023) - CvalYP(cc0012,C234) - 2*D(dd00002)
c	call SOLVE(in)
c	D(dd00123) = in(1)
c	D(dd00223) = .5D0*(D(dd00223) + in(2))
c	D(dd00233) = .5D0*(D(dd00233) + in(3))

***{Ninh change
        D(dd00111)=-CvalYP(cc0,C234)-3*(CvalYP(cc1,C234)
     >	+CvalYP(cc2,C234))
     >  -3*(CvalYP(cc11,C234)+2*CvalYP(cc12,C234)+CvalYP(cc22,C234))
     >  -CvalYP(cc111,C234)-3*CvalYP(cc112,C234)-3*CvalYP(cc122,C234)
     >  -CvalYP(cc222,C234) 
        D(dd00111)=D(dd00111) + 2*m1*D(dd111) - f1*D(dd1111) 
     >  - f2*D(dd1112) - f3*D(dd1113)
        D(dd00111)=D(dd00111)/8D0          
        
        D(dd00112)=CvalYP(cc1,C234)+2*(CvalYP(cc11,C234)
     >	+CvalYP(cc12,C234))
     >  + CvalYP(cc111,C234) + 2*CvalYP(cc112,C234) + CvalYP(cc122,C234)
        D(dd00112)=D(dd00112) + 2*m1*D(dd112) - f1*D(dd1112) 
     >  - f2*D(dd1122) - f3*D(dd1123)
        D(dd00112)=D(dd00112)/8D0          
        
        D(dd00113)=CvalYP(cc2,C234)+2*(CvalYP(cc12,C234)
     >	+CvalYP(cc22,C234))
     >  + CvalYP(cc112,C234) + 2*CvalYP(cc122,C234) + CvalYP(cc222,C234)
        D(dd00113)=D(dd00113) + 2*m1*D(dd113) - f1*D(dd1113) 
     >  - f2*D(dd1123) - f3*D(dd1133)
        D(dd00113)=D(dd00113)/8D0          
        
        D(dd00122)=-CvalYP(cc11,C234)-CvalYP(cc111,C234)
     >	-CvalYP(cc112,C234)
        D(dd00122)=D(dd00122) + 2*m1*D(dd122) - f1*D(dd1122) 
     >  - f2*D(dd1222) - f3*D(dd1223)
        D(dd00122)=D(dd00122)/8D0  
        
        D(dd00123)=-CvalYP(cc12,C234)-CvalYP(cc112,C234)
     >	-CvalYP(cc122,C234)
        D(dd00123)=D(dd00123) + 2*m1*D(dd123) - f1*D(dd1123) 
     >  - f2*D(dd1223) - f3*D(dd1233)
        D(dd00123)=D(dd00123)/8D0          
        
        D(dd00133)=-CvalYP(cc22,C234)-CvalYP(cc122,C234)
     >	-CvalYP(cc222,C234)
        D(dd00133)=D(dd00133) + 2*m1*D(dd133) - f1*D(dd1133) 
     >  - f2*D(dd1233) - f3*D(dd1333)
        D(dd00133)=D(dd00133)/8D0          

        D(dd00222)=CvalYP(cc111,C234)
        D(dd00222)=D(dd00222) + 2*m1*D(dd222) - f1*D(dd1222) 
     >  - f2*D(dd2222) - f3*D(dd2223)
        D(dd00222)=D(dd00222)/8D0  
        
        D(dd00223)=CvalYP(cc112,C234)
        D(dd00223)=D(dd00223) + 2*m1*D(dd223) - f1*D(dd1223) 
     >  - f2*D(dd2223) - f3*D(dd2233)
        D(dd00223)=D(dd00223)/8D0  
        
        D(dd00233)=CvalYP(cc122,C234)
        D(dd00233)=D(dd00233) + 2*m1*D(dd233) - f1*D(dd1233) 
     >  - f2*D(dd2233) - f3*D(dd2333)
        D(dd00233)=D(dd00233)/8D0                  
        
        D(dd00333)=CvalYP(cc222,C234)
        D(dd00333)=D(dd00333) + 2*m1*D(dd333) - f1*D(dd1333) 
     >  - f2*D(dd2333) - f3*D(dd3333)
        D(dd00333)=D(dd00333)/8D0          

c        print*,'dd00i check',D(dd00123),in(1)        
***}Ninh change	
	
	IN(1) = f1*D(dd1111) - csum - 8*D(dd00111)
	IN(2) = f2*D(dd1111) - csum + CvalYP(cc1111,C124)
	IN(3) = f3*D(dd1111) - csum + CvalYP(cc1111,C123)
	call SOLVE(in)
	D(dd11111) = in(1)
	D(dd11112) = in(2)
	D(dd11113) = in(3)

	IN(1) = f1*D(dd2222) - CvalYP(cc1111,C234) + CvalYP(cc1111,C134)
	IN(2) = f2*D(dd2222) - CvalYP(cc1111,C234) - 8*D(dd00222)
	IN(3) = f3*D(dd2222) - CvalYP(cc1111,C234) + CvalYP(cc2222,C123)
	call SOLVE(in)
	D(dd12222) = in(1)
	D(dd22222) = in(2)
	D(dd22223) = in(3)

	IN(1) = f1*D(dd3333) - CvalYP(cc2222,C234) + CvalYP(cc2222,C134)
	IN(2) = f2*D(dd3333) - CvalYP(cc2222,C234) + CvalYP(cc2222,C124)
	IN(3) = f3*D(dd3333) - CvalYP(cc2222,C234) - 8*D(dd00333)
	call SOLVE(in)
	D(dd13333) = in(1)
	D(dd23333) = in(2)
	D(dd33333) = in(3)

	IN(1) = f1*D(dd1122) - c11sum - 4*D(dd00122)
	IN(2) = f2*D(dd1122) - c11sum - 4*D(dd00112)
	IN(3) = f3*D(dd1122) - c11sum + CvalYP(cc1122,C123)
	call SOLVE(in)
	D(dd11122) = in(1)
	D(dd11222) = in(2)
	D(dd11223) = in(3)

	IN(1) = f1*D(dd1133) - c22sum - 4*D(dd00133)
	IN(2) = f2*D(dd1133) - c22sum + CvalYP(cc1122,C124)
	IN(3) = f3*D(dd1133) - c22sum - 4*D(dd00113)
	call SOLVE(in)
	D(dd11133) = in(1)
	D(dd11233) = in(2)
	D(dd11333) = in(3)

	IN(1) = f1*D(dd2233) - CvalYP(cc1122,C234) + CvalYP(cc1122,C134)
	IN(2) = f2*D(dd2233) - CvalYP(cc1122,C234) - 4*D(dd00233)
	IN(3) = f3*D(dd2233) - CvalYP(cc1122,C234) - 4*D(dd00223)
	call SOLVE(in)
	D(dd12233) = in(1)
	D(dd22233) = in(2)
	D(dd22333) = in(3)

	IN(1) = f1*D(dd1123) - c12sum - 4*D(dd00123)
	IN(2) = f2*D(dd1123) - c12sum - 2*D(dd00113)
	IN(3) = f3*D(dd1123) - c12sum - 2*D(dd00112)
	call SOLVE(in)
	D(dd11123) = in(1)
	D(dd11223) = .5D0*(D(dd11223) + in(2))
	D(dd11233) = .5D0*(D(dd11233) + in(3))

	IN(1) = f1*D(dd2223) - CvalYP(cc1112,C234) + CvalYP(cc1112,C134)
	IN(2) = f2*D(dd2223) - CvalYP(cc1112,C234) - 6*D(dd00223)
	IN(3) = f3*D(dd2223) - CvalYP(cc1112,C234) - 2*D(dd00222)
	call SOLVE(in)
	D(dd12223) = in(1)
	D(dd22223) = .5D0*(D(dd22223) + in(2))
	D(dd22233) = .5D0*(D(dd22233) + in(3))

	IN(1) = f1*D(dd2333) - CvalYP(cc1222,C234) + CvalYP(cc1222,C134)
	IN(2) = f2*D(dd2333) - CvalYP(cc1222,C234) - 2*D(dd00333)
	IN(3) = f3*D(dd2333) - CvalYP(cc1222,C234) - 6*D(dd00233)
	call SOLVE(in)
	D(dd12333) = in(1)
	D(dd22333) = .5D0*(D(dd22333) + in(2))
	D(dd23333) = .5D0*(D(dd23333) + in(3))

	endif
***]PV rank5	


***]using Passarino-Veltman method

	end
************************************************************************
      subroutine CYPRotate24(rp1,rp2,rp3,rp4,rp1p2,rp2p3,
     > cm1,cm2,cm3,cm4,i)
************************************************************************
* Ref0. Dao Thi Nhung and Le Duc Ninh arXiv: 0902.0325v2
* Author: Le Duc Ninh (2010)
************************************************************************
      IMPLICIT NONE
      
      CVAR rp1,rp2,rp3,rp4,rp1p2,rp2p3
      CVAR cm1,cm2,cm3,cm4
      CVAR temp
      integer i
#include "lt.h"

* 1
	if ( i .le. 1 ) return

* 2
	if ( i .eq. 2 ) then
	 temp=rp2
	 rp2=rp2p3
	 rp2p3=temp
	 temp=rp4
	 rp4=rp1p2
	 rp1p2=temp
	 
	 temp=cm3
	 cm3=cm4
	 cm4=temp
         return
	endif	 
* 3
	if ( i .eq. 3 ) then
       temp=rp1
       rp1=rp1p2
       rp1p2=temp
       temp=rp3
       rp3=rp2p3
       rp2p3=temp
       
       temp=cm2
       cm2=cm3
       cm3=temp
	return
	endif

* 4
	if ( i .eq. 4 ) then       
       temp=rp2
       rp2=rp3
       rp3=rp2p3
       rp2p3=temp
       temp=rp4 
       rp4=rp1
       rp1=rp1p2
       rp1p2=temp
       
       temp=cm2
       cm2=cm3
       cm3=cm4
       cm4=temp
	return
	endif

** 5
	if ( i .eq. 5 ) then
	 temp = rp1
	 rp1 = rp4
	 rp4 = rp1p2
	 rp1p2 = temp
	 temp = rp2
	 rp2 = rp2p3
	 rp2p3 = rp3
	 rp3 = temp
	 
	 temp = cm2
	 cm2 = cm4
	 cm4 = cm3
	 cm3 = temp 
	return
	endif

* 6
	if ( i .eq. 6 ) then
       temp=rp1
       rp1=rp4
       rp4=temp
       temp=rp2
       rp2=rp3
       rp3=temp
       
       temp=cm2
       cm2=cm4
       cm4=temp
	return
	endif

** 7
	if ( i .eq. 7 ) then
	 temp = rp2
	 rp2 = rp1p2
	 rp1p2 = temp
	 temp = rp4
	 rp4 = rp2p3
	 rp2p3 = temp
	 
	 temp = cm1
	 cm1 = cm2
	 cm2 = temp
	return
	endif
** 8
	if ( i .eq. 8 ) then
	 temp = rp2
	 rp2 = rp4
	 rp4 = temp
	 temp = rp1p2
	 rp1p2 = rp2p3
	 rp2p3 = temp
	 
	 temp = cm1
	 cm1 = cm2
	 cm2 = temp
	 temp = cm3
	 cm3 = cm4
	 cm4 = temp
	return
	endif
* 9
	if ( i .eq. 9 ) then
       temp=rp1
       rp1=rp2
       rp2=rp1p2
       rp1p2=temp
       temp=rp3
       rp3=rp4
       rp4=rp2p3
       rp2p3=temp
       
       temp=cm1
       cm1=cm2
       cm2=cm3
       cm3=temp
	return
	endif
* 10
	if ( i .eq. 10 ) then
       temp=rp1
       rp1=rp2
       rp2=rp3
       rp3=rp4
       rp4=temp
       temp=rp1p2
       rp1p2=rp2p3
       rp2p3=temp
       
       temp=cm1
       cm1=cm2
       cm2=cm3
       cm3=cm4
       cm4=temp
	return
	endif
* 11
	if ( i .eq. 11 ) then
       temp=rp1
       rp1=rp2p3
       rp2p3=rp3
       rp3=rp1p2
       rp1p2=temp
       temp=rp2
       rp2=rp4
       rp4=temp
       
       temp=cm1
       cm1=cm2
       cm2=cm4
       cm4=cm3
       cm3=temp
       return
       endif

* 12
	if ( i .eq. 12 ) then
       temp=rp1
       rp1=rp2p3
       rp2p3=rp4
       rp4=temp
       temp=rp2
       rp2=rp3
       rp3=rp1p2
       rp1p2=temp
       
       temp=cm1
       cm1=cm2
       cm2=cm4
       cm4=temp
	return
	endif

** 13
	if ( i .eq. 13 ) then
	 temp = rp1
	 rp1 = rp1p2
	 rp1p2 = rp2
	 rp2 = temp
	 temp = rp3
	 rp3 = rp2p3
	 rp2p3 = rp4
	 rp4 = temp
	 
	 temp = cm1
	 cm1 = cm3
	 cm3 = cm2
	 cm2 = temp
	 return
	 endif
* 14
	if ( i .eq. 14 ) then
       temp=rp1
       rp1 = rp1p2
       rp1p2 = rp3
       rp3 = rp2p3
       rp2p3 = temp
       temp = rp2
       rp2 = rp4
       rp4 = temp

       temp = cm1
       cm1 = cm3
       cm3 = cm4
       cm4 = cm2
       cm2 = temp
	return
	endif

** 15
	if ( i .eq. 15 ) then
	 temp = rp1
	 rp1 = rp2
	 rp2 = temp
	 temp = rp3
	 rp3 = rp4
	 rp4 = temp
	 
	 temp = cm1
	 cm1 = cm3
	 cm3 = temp
	 return
	 endif
** 16
	if ( i .eq. 16 ) then
	 temp = rp1
	 rp1 = rp2
	 rp2 = rp2p3
	 rp2p3 = temp
	 temp = rp3
	 rp3 = rp4
	 rp4 = rp1p2
	 rp1p2 = temp
	 
	 temp = cm1
	 cm1 = cm3
	 cm3 = cm4
	 cm4 = temp
	return
	endif
* 17
	if ( i .eq. 17 ) then
       temp=rp1
       rp1=rp3
       rp3=temp
       temp=rp2
       rp2=rp4
       rp4=temp

       temp=cm1
       cm1=cm3
       cm3=temp
       temp=cm2
       cm2=cm4
       cm4=temp 
	return
	endif

** 18
	if ( i .eq. 18 ) then
	 temp = rp1
	 rp1 = rp3
	 rp3 = temp
	 temp = rp2
	 rp2 = rp2p3
	 rp2p3 = rp4
	 rp4 = rp1p2
	 rp1p2 = temp
	 
	 temp = cm1
	 cm1 = cm3
	 cm3 = cm2
	 cm2 = cm4
	 cm4 = temp
	 return
	 endif
* 19 
	if ( i .eq. 19 ) then
       temp=rp1
       rp1=rp4
       rp4=rp3
       rp3=rp2
       rp2=temp
       temp=rp1p2
       rp1p2=rp2p3
       rp2p3=temp

       temp=cm1
       cm1=cm4
       cm4=cm3
       cm3=cm2
       cm2=temp     
	return
	endif

** 20
	if ( i .eq. 20 ) then
	 temp = rp1
	 rp1 = rp4
	 rp4 = rp2p3
	 rp2p3 = temp
	 temp = rp2
	 rp2 = rp1p2
	 rp1p2 = rp3
	 rp3 = temp
	 
	 temp = cm1
	 cm1 = cm4
	 cm4 = cm2
	 cm2 = temp
	return
	endif
** 21
	if ( i .eq. 21 ) then
	 temp = rp1
	 rp1 = rp2p3
	 rp2p3 = rp2
	 rp2 = temp
	 temp = rp3
	 rp3 = rp1p2
	 rp1p2 = rp4
	 rp4 = temp
	 
	 temp = cm1
	 cm1 = cm4
	 cm4 = cm3
	 cm3 = temp
	 return
	 endif
** 22
	if ( i .eq. 22 ) then
	 temp = rp1
	 rp1 = rp2p3
	 rp2p3 = temp
	 temp = rp3
	 rp3 = rp1p2
	 rp1p2 = temp
	 
	 temp = cm1
	 cm1 = cm4
	 cm4 = temp	 	 	 	 	 
	 return
	 endif	 	  

* 23 
	if ( i .eq. 23 ) then
	 temp=rp1
	 rp1=rp3
	 rp3=temp
	 temp=rp2
	 rp2=rp1p2
	 rp1p2=rp4
	 rp4=rp2p3
	 rp2p3=temp
	 
	 temp=cm1
	 cm1=cm4
	 cm4=cm2
	 cm2=cm3
	 cm3=temp
	 return
	 endif
* 24
	if ( i .eq. 24 ) then
	 temp=rp1
	 rp1=rp3
	 rp3=temp
	 temp=rp1p2
	 rp1p2=rp2p3
	 rp2p3=temp
	 
	 temp=cm1
	 cm1=cm4
	 cm4=temp
	 temp=cm2
	 cm2=cm3
	 cm3=temp
	 return
	endif
	end 
