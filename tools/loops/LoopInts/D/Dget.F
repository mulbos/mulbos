* Dget.F
* retrieve the four-point integrals

#include "defs.h"

	integer function XYPDget(p1, p2, p3, p4, p1p2, p2p3,
     &    m1, m2, m3, m4)
	implicit none
	XVAR p1, p2, p3, p4, p1p2, p2p3
	XVAR m1, m2, m3, m4

#include "lt.h"

	integer cachelookup, cachelookupqp
	external cachelookup, cachelookupqp
	external CYPDcoeff

	CVAR para(Pdd)
#ifdef NO_CACHE
	CVAR D(Ndd)
	integer i
#endif

	para(1) = p1
	para(2) = p2
	para(3) = p3
	para(4) = p4
	para(5) = p1p2
	para(6) = p2p3
	para(7) = m1
	para(8) = m2
	para(9) = m3
	para(10) = m4

#ifndef NO_CACHE

#ifndef QUADPREC
	XYPDget=cachelookup(para, DvalYP(1,0), CYPDcoeff, RC*Pdd, Ndd)
#else
	XYPDget=cachelookupqp(para, DvalYP(1,0), CYPDcoeff, RC*Pdd, Ndd)
#endif

#else
* no cache here
        iDget = iDget + 1 + Ndd
        XYPDget = iDget
	call CYPDcoeff(para, D, 1)
	do i = 1, Ndd
	DvalYP(i, XYPDget) = D(i)
	enddo		
#endif

#ifdef DEBUG1
	print *, "XYPDget in Dget.F:"
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
#endif
	end


