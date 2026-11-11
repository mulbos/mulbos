* Dump.F
* dumps the parameters and coefficients on screen
* this file is part of LoopTools
* last modified 8 Dec 05 th
* adapted by LE Duc Ninh (2011)

#include "defs.h"


	subroutine CYPDumpPara(npoint, para, ldpara, origin)
	implicit none
	integer npoint, ldpara
	XVAR para(ldpara,*)
	character*(*) origin

c#include "lt.h"

	integer i

	character*6 paraname(Pff,2:6)
	common /YPparanames/ paraname

	integer npara(2:6)
	data npara /Pbb, Pcc, Pdd, Pee, Pff/

	print *, origin

	do i = 1, npara(npoint)
	  print *, "  ", paraname(i,npoint), "=", para(1,i)
	enddo

	call flush(6)
	end

************************************************************************

	subroutine CYPDumpCoeff(npoint, coeff)
	implicit none
	integer npoint
	CVAR coeff(*)

c#include "lt.h"

	integer i

	character*8 coeffname(Nff,2:6)
	common /YPcoeffnames/ coeffname

	integer ncoeff(2:6)
	data ncoeff /Nbb, Ncc, Ndd, Nee, Nff/

	do i = 1, ncoeff(npoint)
	  print *, coeffname(i,npoint), "=", coeff(i)
	enddo

	print *, "===================================================="

	call flush(6)
	end
