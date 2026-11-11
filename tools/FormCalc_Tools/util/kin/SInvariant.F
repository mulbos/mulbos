* SInvariant.F
* calculates an s-type invariant, i.e. s_{ij} = (k_i + k_j)^2
* this file is part of FormCalc
* modified by Le Duc Ninh 28 Oct 2025

	double precision function SInvariant(i, j)
	implicit none
	integer i, j
#include "vecset.h"

	double precision S11, S22
	double complex S12

	S11 = DBLE(vec_K(1,1,i)) + DBLE(vec_K(1,1,j))
	S22 = DBLE(vec_K(2,2,i)) + DBLE(vec_K(2,2,j))
	S12 = vec_K(1,2,i) + vec_K(1,2,j)

	SInvariant = S11*S22 - DBLE(S12)**2 - DIMAG(S12)**2
	end

