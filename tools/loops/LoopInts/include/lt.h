* lt.h
* Internal common blocks for the cache system
* Declaration of popular functions and subroutines

#include "ff.h"

* the cache-pointer structure is (see cache.c):
* 1. int valid
* 2. Node *last
* 3. Node *first
* 4. (not used)

* number of cache layers for 2, 3, 4, 5, 6 point integrals; 
* each integral needs two layers (see RC = 1,2 in defs.h) 
* for real and complex masses.
	integer ncaches
	parameter (ncaches = 10)

#ifndef NO_CACHE
* KIND = 1 for double prec., 2 for quadruple
	integer*8 cacheptr(KIND,4,ncaches)
	integer*8 savedptr(2,ncaches)

	common /ltvars/
     &    cacheptr, savedptr

* 3-KIND here because cache() must have the same 
* storage space as cacheptr()
	CVAR cache(3-KIND,ncaches)
	equivalence (cacheptr, cache)

***[for quadruple precision only
	integer*8 cacheptrQP(4,ncaches)
	integer*8 savedptrQP(2,ncaches)

	common /ltvarsQP/
     &    cacheptrQP, savedptrQP

	CVAR cacheQP(1,ncaches)
	equivalence (cacheptrQP, cacheQP)
***]for quadruple precision only
#endif

#ifdef NO_CACHE
* no cache here
	integer iBget, iCget, iDget, iEget, iFget
	integer iBgetsv, iCgetsv, iDgetsv, iEgetsv, iFgetsv	

	common /YPltiget/ iBget, iCget, iDget, iEget, iFget,
     &  iBgetsv, iCgetsv, iDgetsv, iEgetsv, iFgetsv

* Max number of loop integrals for each layer 
* (ncaches = 10 layers, read above). 
* This must be the same in looptools.h.
	integer nints
	parameter (nints = 10000)	

	CVAR cache(nints,ncaches)
	CVAR cacheQP(nints,ncaches)
	common /ltvars/ cache
	common /ltvarsQP/ cacheQP
#endif

***[popular variables, functions and subroutines
	CVAR CYPLOG, CYPeta, CYPetax, CYPSp, CYPLi2n, CYPLi2nx
	external CYPLOG, CYPeta, CYPetax, CYPSp, CYPLi2n, CYPLi2nx

	external CYPRoot2, CYPKam
	integer YPirk, YPDltB
	external YPirk, YPDltB

	RVAR pmax, mmax, pmmax
	RVAR pmin, mmin, pmmin
***]popular variables, functions and subroutines
