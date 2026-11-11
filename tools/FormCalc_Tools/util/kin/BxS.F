* BxS.F
* computes Bar[Vector].Spinor
* this file is part of FormCalc
* modified by Le Duc Ninh 28 Oct 2025

	function BxS(A, B)
	implicit none

	double complex A(2,2), B(2)
	double complex BxS(2)

	BxS(1) = A(2,2)*B(1) - A(1,2)*B(2)
	BxS(2) = A(1,1)*B(2) - A(2,1)*B(1)
	end

