* SxS.F
* computes Spinor.Spinor
* this file is part of FormCalc
* modified by Le Duc Ninh 28 Oct 2025

	double complex function SxS(A, B)
	implicit none
	double complex A(2), B(2)

	SxS = A(1)*B(1) + A(2)*B(2)
	end

