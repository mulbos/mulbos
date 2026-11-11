* ini.F

#include "defs.h"

	subroutine YPclearcache
	implicit none

#include "lt.h"

#ifndef NO_CACHE
	integer i

	do i = 1, ncaches
#ifndef QUADPREC
	  cacheptr(1,1,i) = 0
	  cacheptr(1,2,i) = 0
#else
	  cacheptrQP(1,i) = 0
	  cacheptrQP(2,i) = 0
#endif
	enddo

#else
        iBget = 0
        iCget = 0
        iDget = 0
        iEget = 0
        iFget = 0 
#endif
	end

************************************************************************

	subroutine YPmarkcache
	implicit none

#include "lt.h"

#ifndef NO_CACHE
	integer i

	do i = 1, ncaches
#ifndef QUADPREC
	  savedptr(1,i) = cacheptr(1,1,i)
	  savedptr(2,i) = cacheptr(1,2,i)
#else
	  savedptrQP(1,i) = cacheptrQP(1,i)
	  savedptrQP(2,i) = cacheptrQP(2,i)
#endif
	enddo

#else
        iBgetsv = iBget
        iCgetsv = iCget
        iDgetsv = iDget
        iEgetsv = iEget
        iFgetsv = iFget 	
#endif
	end

************************************************************************

	subroutine YPrestorecache
	implicit none

#include "lt.h"

#ifndef NO_CACHE
	integer i

	do i = 1, ncaches
#ifndef QUADPREC
	  cacheptr(1,1,i) = savedptr(1,i)
	  cacheptr(1,2,i) = savedptr(2,i)
#else
	  cacheptrQP(1,i) = savedptrQP(1,i)
	  cacheptrQP(2,i) = savedptrQP(2,i)
#endif
	enddo

#else
        iBget = iBgetsv
        iCget = iCgetsv
        iDget = iDgetsv
        iEget = iEgetsv
        iFget = iFgetsv
#endif
	end

************************************************************************

	subroutine YPsetmudim(mudim_)
	implicit none
	RVAR mudim_

#include "lt.h"

	mudim = mudim_
	call YPclearcache
	end

************************************************************************

	RVAR function YPgetmudim()
	implicit none

#include "lt.h"

	YPgetmudim = mudim
	end

************************************************************************

	subroutine YPsetdelta(delta_)
	implicit none
	RVAR delta_

#include "lt.h"

	delta = delta_
	call YPclearcache
	end

************************************************************************

	RVAR function YPgetdelta()
	implicit none

#include "lt.h"

	YPgetdelta = delta
	end

************************************************************************

	subroutine YPsetlambda(lambda_)
	implicit none
	RVAR lambda_

#include "lt.h"

	lambda = lambda_
	if ( lambda .lt. 0) xrat = 0
	call YPclearcache
	end

************************************************************************

	RVAR function YPgetlambda()
	implicit none

#include "lt.h"

	YPgetlambda = lambda
	end

************************************************************************
