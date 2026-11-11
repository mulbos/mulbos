* Bget.F
* retrieve the two-point tensor coefficients

#include "defs.h"

	integer function XYPBget(p, m1, m2)
	implicit none
	XVAR p, m1, m2

#include "lt.h"

	integer cachelookup, cachelookupqp
	external cachelookup, cachelookupqp, CYPBcoeff

	CVAR para(Pbb)
#ifdef NO_CACHE
	CVAR B(Nbb)
	integer i
#endif

	para(1) = p
	para(2) = m1
	para(3) = m2

#ifndef NO_CACHE

#ifndef QUADPREC
	XYPBget=cachelookup(para, BvalYP(1,0), CYPBcoeff, RC*Pbb, Nbb)
#else
	XYPBget=cachelookupqp(para, BvalYP(1,0), CYPBcoeff, RC*Pbb, Nbb)
#endif

#else
* no cache here
        iBget = iBget + 1 + Nbb
        XYPBget = iBget
	call CYPBcoeff(para, B, 1)
	do i = 1, Nbb
	BvalYP(i, XYPBget) = B(i)
	enddo
#endif

#ifdef DEBUG1
	print *, "XYPBget in Bget.F:"
	print *, "    p =", para(1)
	print *, "    m1 =", para(2)
	print *, "    m2 =", para(3)
#endif
	
	end
************************************************************************

	CVAR function XYPB0i(i, p, m1, m2)
	implicit none
	integer i
	XVAR p, m1, m2

#include "lt.h"

	integer XYPBget
	external XYPBget

	integer b

	b = XYPBget(p, m1, m2)
	XYPB0i = BvalYP(i,b)
	end
