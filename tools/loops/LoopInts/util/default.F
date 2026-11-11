#include "defs.h"

	subroutine YPDefault

#include "ff.h"

***[DEBUG
        debugkey = 0
***]DEBUG

***[run mode
* possible values: 
* 0 (best run) 
* 1 (print small problems)
* 2 (print small + big problems) 
* ...
* 9 (print all problems)
	runmode = 0
***]run mode

***[version comparison
	versionkey = 0
	maxdev = 1D-10
***[version comparison

***[error tracking
	iers = 1
***]error tracking

***[Max rank for tensor coefficients
* default
	rankmax3 = 5
	rankmax4 = 5
	rankmax5 = 5
	rankmax6 = 5 ! rank 6 not available
***]Max rank for tensor coefficients

***[regulation parameters
	mudim = 1
	delta = 0
	lambda = 1

* switch (1 or 0) for rational terms in tensor coefficients
	xrat = 1
***]regulation parameters

***[common approximation parameters for DP and QP
* can be used in CYPSetM0()
* legacy parameter now
	precd = 1D-15
***]common approximation parameters for DP and QP

***[instability parameters
* approximation parameter for collinear scalar integrals, 
* it can depend on prec:
cc	stab0 = 1D-12
	stab0 = prec
* cut-off (no approximation) parameters:
* stab2: fixed to 1D-6, problem independence
	stab2 = 1D-6
* These cut-offs are problem independent
	stab3 = 1D-3
	stab4 = 1D-3
	stab5 = 1D-10
	stab6 = 1D-12
***]instability parameters

	end

***********************************************************************
      subroutine YPSetPrec
      implicit none
      RVAR xx, yy, zz
      integer i, Nmax

#include "ff.h"

      xx = 1
      zz = 0

      Nmax = 1000

 1    continue

      do i = 1, Nmax
       xx = xx/2

cc       yy = log(1 + xx)
	yy = 1 + xx

       if ( yy .eq. zz) goto 2
       zz = yy 
      enddo

      if ( i .ge. Nmax) goto 1 

 2    continue

* take three bits for safety
      prec = xx*8

      print*,'MESSAGE from YPSetPrec: precision set to', prec

      end

************************************************************************
