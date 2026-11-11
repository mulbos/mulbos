* Pair.F
* calculates the scalar product of two four-vectors
* this file is part of FormCalc
* modified by Le Duc Ninh 28 Oct 2025

	double complex function Pair(A, B)
	implicit none

	double complex A(2,2), B(2,2)

	Pair = .5D0*(
     &    DBLE(A(1,1))*DBLE(B(2,2)) +
     &    DBLE(A(2,2))*DBLE(B(1,1)) -
     &    A(1,2)*B(2,1) - A(2,1)*B(1,2) )
	end

