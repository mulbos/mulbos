	subroutine GetPolarVec(i1,epsp,epspc)
!	use global_vars
!	use global_mom_xpdf
! Input: particle intex i1
! Output: polarization vector and its complex conjugate
	implicit none
	integer i1, j1
      	double complex epsp(-1:1,4), epspc(-1:1,4)
      	double complex cI
      	parameter (cI = (0D0, 1D0))

#include "../FormCalc_Tools/vecset.h"

        epsp = (0D0, 0D0)
        epspc = (0D0, 0D0)
        
***[ get the polarization vectors
	epsp(1,1) = (vec_E(1,1,1,i1) + vec_E(2,2,1,i1))/2D0
	epsp(1,4) = (vec_E(1,1,1,i1) - vec_E(2,2,1,i1))/2D0
	epsp(1,2) = (vec_E(2,1,1,i1) + vec_E(1,2,1,i1))/2D0
* convention for eps K(Adot, B)
	epsp(1,3) = cI*(vec_E(2,1,1,i1) - vec_E(1,2,1,i1))/2D0

	epsp(-1,1) = (vec_E(1,1,-1,i1) + vec_E(2,2,-1,i1))/2D0
	epsp(-1,4) = (vec_E(1,1,-1,i1) - vec_E(2,2,-1,i1))/2D0
	epsp(-1,2) = (vec_E(2,1,-1,i1) + vec_E(1,2,-1,i1))/2D0
	epsp(-1,3) = cI*(vec_E(2,1,-1,i1) - vec_E(1,2,-1,i1))/2D0
* the complex conjugate
	do j1 = 1, 4
	 epspc(1,j1) = DCONJG(epsp(1,j1))
	 epspc(-1,j1) = DCONJG(epsp(-1,j1))
	enddo
***] get the polarization vectors
        return
	end subroutine
