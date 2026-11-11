* VxS.F
* computes Vector.Spinor
* this file is part of FormCalc
* modified by Le Duc Ninh 28 Oct 2025

	function VxS(A, B)
	implicit none

	double complex A(2,2), B(2)
	double complex VxS(2)

	VxS(1) = A(1,1)*B(1) + A(1,2)*B(2)
	VxS(2) = A(2,1)*B(1) + A(2,2)*B(2)
	end

