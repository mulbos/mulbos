* Eps.F
* calculates -I times the Levi-Civita tensor contracted with
* four four-vectors
* this file is part of FormCalc
* modified by Le Duc Ninh 28 Oct 2025

	double complex function Eps(A, B, C, D)
	implicit none

	double complex A(2,2), B(2,2), C(2,2), D(2,2)

#define diag(i, j) (DBLE(i(1,1))*DBLE(j(2,2)) - DBLE(i(2,2))*DBLE(j(1,1)))
#define offdiag(i, j) (i(2,1)*j(1,2) - i(1,2)*j(2,1))

	Eps = .25D0*(
     &    diag(A, B)*
     &      offdiag(C, D) +
     &    diag(C, D)*
     &      offdiag(A, B) -
     &    diag(A, C)*
     &      offdiag(B, D) -
     &    diag(B, D)*
     &      offdiag(A, C) +
     &    diag(A, D)*
     &      offdiag(B, C) +
     &    diag(B, C)*
     &      offdiag(A, D) )
	end

