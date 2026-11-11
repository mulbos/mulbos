* SeS.F
* computes Spinor.eps.Spinor
* this file is part of FormCalc
* modified by Le Duc Ninh 28 Oct 2025

	double complex function SeS(A, B)
	implicit none
	double complex A(2), B(2)

	SeS = A(1)*B(2) - A(2)*B(1)
	end

