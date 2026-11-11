* VeS.F
* computes Vector.eps.Spinor
* this file is part of FormCalc
* modified by Le Duc Ninh 28 Oct 2025

	function VeS(A, B)
	implicit none

	double complex A(2,2), B(2), VeS(2)

	VeS(1) = A(1,1)*B(2) - A(1,2)*B(1)
	VeS(2) = A(2,1)*B(2) - A(2,2)*B(1)

	end

