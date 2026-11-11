* TInvariant.F
* calculates an t-type invariant, i.e. t_{ij} = (k_i - k_j)^2
* this file is part of FormCalc
* modified by Le Duc Ninh 28 Oct 2025

	double precision function TInvariant(i, j)
	implicit none
	integer i, j
#include "vecset.h"

	double precision T11, T22
	double complex T12

	T11 = DBLE(vec_K(1,1,i)) - DBLE(vec_K(1,1,j))
	T22 = DBLE(vec_K(2,2,i)) - DBLE(vec_K(2,2,j))
	T12 = vec_K(1,2,i) - vec_K(1,2,j)

        TInvariant = T11*T22 - DBLE(T12)**2 - DIMAG(T12)**2
	end

