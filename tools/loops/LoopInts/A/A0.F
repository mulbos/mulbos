* A0.F
* the scalar one-point function

#include "defs.h"

	CVAR function XYPA0(m)
	implicit none
	XVAR m
	CVAR mC

#include "lt.h"

	CVAR CYPA0s
	external CYPA0s

	CVAR res(0:1)
	integer key, ier

	mC = m
	res(0) = CYPA0s(mC)
	
	XYPA0 = res(0)
	end

****************************************************************
	CVAR function XYPA00(m)
	implicit none
	XVAR m

#include "lt.h"

	CVAR CYPA00s, mC
	external CYPA00s

	mC = m
	XYPA00 = CYPA00s(mC)
	end

************************************************************************
* adapter code for C++
