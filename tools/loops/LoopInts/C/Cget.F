* Cget.F
* the three-point integrals

#include "defs.h"

	integer function XYPCget(p1, p2, p1p2, m1, m2, m3)
	implicit none
	XVAR p1, p2, p1p2, m1, m2, m3

#include "lt.h"

	integer cachelookup, cachelookupqp
	external cachelookup, cachelookupqp
	external CYPCcoeff
	CVAR para(Pcc)

#ifdef NO_CACHE
	CVAR C(Ncc)
	integer i
#endif

	para(1) = p1
	para(2) = p2
	para(3) = p1p2
	para(4) = m1
	para(5) = m2
	para(6) = m3

#ifndef NO_CACHE

#ifndef QUADPREC
	XYPCget=cachelookup(para, CvalYP(1,0), CYPCcoeff, RC*Pcc, Ncc)
#else
	XYPCget=cachelookupqp(para, CvalYP(1,0), CYPCcoeff, RC*Pcc, Ncc)
#endif

#else
* no cache here
        iCget = iCget + 1 + Ncc
        XYPCget = iCget
	call CYPCcoeff(para, C, 1)
	do i = 1, Ncc
	CvalYP(i, XYPCget) = C(i)
	enddo		
#endif

#ifdef DEBUG1
	print *, "XYPCget in Cget.F:"
	print *, "    p1   =", para(1)
	print *, "    p2   =", para(2)
	print *, "    p1p2 =", para(3)
	print *, "    m1 =", para(4)
	print *, "    m2 =", para(5)
	print *, "    m3 =", para(6)
#endif

	end


