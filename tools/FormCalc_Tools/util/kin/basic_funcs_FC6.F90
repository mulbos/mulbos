! adapted from FormCalc-6 by Le Duc Ninh
! simple interface to FormCalc-6 amplitudes, 
! removed the dangerous equivalence (vec, ...)
! last modified by Le Duc Ninh 28 Oct 2025

module basic_funcs_FC6
implicit none
contains
	function VxS(A, B) result(V)
	implicit none

	double complex, intent(in) :: A(2,2), B(2)
	double complex :: V(2)

	V(1) = A(1,1)*B(1) + A(1,2)*B(2)
	V(2) = A(2,1)*B(1) + A(2,2)*B(2)
	end function VxS

	function VeS(A, B) result(V)
	implicit none

	double complex, intent(in) :: A(2,2), B(2)
	double complex :: V(2)

	V(1) = A(1,1)*B(2) - A(1,2)*B(1)
	V(2) = A(2,1)*B(2) - A(2,2)*B(1)

	end function VeS

	function BxS(A, B) result(V)
	implicit none

	double complex, intent(in) :: A(2,2), B(2)
	double complex :: V(2)

	V(1) = A(2,2)*B(1) - A(1,2)*B(2)
	V(2) = A(1,1)*B(2) - A(2,1)*B(1)
	end function BxS

	function BeS(A, B) result(V)
	implicit none

	double complex, intent(in) :: A(2,2), B(2)
	double complex :: V(2)

	V(1) = A(1,2)*B(1) + A(2,2)*B(2)
	V(2) = -A(1,1)*B(1) - A(2,1)*B(2)
	end function BeS

	function SxS(A, B)
	implicit none
	double complex, intent(in) :: A(2), B(2)
	double complex :: SxS

	SxS = A(1)*B(1) + A(2)*B(2)
	end function SxS

	function SeS(A, B)
	implicit none
	double complex, intent(in) :: A(2), B(2)
	double complex :: SeS

	SeS = A(1)*B(2) - A(2)*B(1)
	end function SeS

	function Pair(A, B)
	implicit none

	double complex, intent(in) :: A(2,2), B(2,2)
	double complex :: Pair

	Pair = .5D0*(&
     &    DBLE(A(1,1))*DBLE(B(2,2)) +&
     &    DBLE(A(2,2))*DBLE(B(1,1)) -&
     &    A(1,2)*B(2,1) - A(2,1)*B(1,2) )
	end function Pair

	function Eps(A, B, C, D)
	implicit none

	double complex, intent(in) :: A(2,2), B(2,2), C(2,2), D(2,2)
	double complex :: Eps

#define diag(i, j) (DBLE(i(1,1))*DBLE(j(2,2)) - DBLE(i(2,2))*DBLE(j(1,1)))
#define offdiag(i, j) (i(2,1)*j(1,2) - i(1,2)*j(2,1))

	Eps = .25D0*(&
     &    diag(A, B)*&
     &      offdiag(C, D) +&
     &    diag(C, D)*&
     &      offdiag(A, B) -&
     &    diag(A, C)*&
     &      offdiag(B, D) -&
     &    diag(B, D)*&
     &      offdiag(A, C) +&
     &    diag(A, D)*&
     &      offdiag(B, C) +&
     &    diag(B, C)*&
     &      offdiag(A, D) )
	end function Eps


end module basic_funcs_FC6
