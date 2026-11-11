* born.F
* form factors for {Men(offshell W)} -> {e+, nu_e, A} = i1 -> i2 + i3 + i4

	subroutine decayW2lA(amp, i1, i2, i3, i4)
	use basic_funcs_FC6
	implicit none
	double complex amp
	double complex Fx1, Fx3, Fx4, Fx2, Paix6, Paix5, Paix3, Paix2
	double complex Paix1, Paix4, Abx1, Abx2, Abx3, Abx4, AbxSum1
	double complex AbxSum2
        double precision Tx, Sx23

	integer i1, i2, i3, i4

#include "vars.h"
C checking
c#define CHECK_INPUT
#ifdef CHECK_INPUT
        print*,"check 1 =",Pair(k(i1),k(i1))/MW2
        print*,"check 0 =",Pair(k(i2),k(i2))/MW2
        print*,"check 0 =",Pair(k(i3),k(i3))/MW2
        print*,"check 0 =",Pair(k(i4),k(i4))/MW2
        print*,"check 1 =",2*(Pair(k(i2),k(i3))
     & +Pair(k(i2),k(i4))
     & +Pair(k(i3),k(i4)))/MW2
        print*,"check e11 =",Pair(e(i1),k(i1))/MW2
        print*,"check e14 =",Pair(e(i1),k(i4))/MW2
        print*,"check e41 =",Pair(e(i4),k(i1))/MW2
        print*,"check e44 =",Pair(e(i4),k(i4))/MW2
#endif

	Tx = TInvariant(i1,i3)

	Sx23 = SInvariant(i2,i3)

	Fx1 = Hel(i2)*SxS(DottedSpinor(i3,1,7),VxS(e(i1),Spinor(i2,-1,7)))

	Fx3 = Hel(i2)*SxS(DottedSpinor(i3,1,7),VxS(ec(i4),Spinor(i2,-1,7)))

	Fx4 = Hel(i2)*SxS(DottedSpinor(i3,1,7),VxS(k(i1),Spinor(i2,-1,7)))

	Paix4 = Pair(e(i1),ec(i4))

	Paix3 = Pair(e(i1),k(i2))

	Paix2 = Pair(e(i1),k(i3))

	Paix1 = Pair(ec(i4),k(i1))

        Fx2 = Hel(i2)*SxS(DottedSpinor(i3,1,7),
     &    VxS(e(i1),BxS(ec(i4),VxS(k(i1),Spinor(i2,-1,7)))))

	Abx1 = Fx1*Paix1

	Abx2 = Fx3*Paix2

	Abx3 = Fx3*Paix3

	Abx4 = Fx4*Paix4

	AbxSum1 = Abx1 - Abx2 - 1/2.D0*Fx2

	AbxSum2 = Abx1 - Abx2 - Abx3 - Abx4

        amp = 
     &   4*Alfa*Pi*((AbxSum2/(-MW2 + Sx23) + AbxSum1/(-ME2 + Tx))*
     &       Sqrt(2.D0))/SW


#ifdef DEBUG
	DEB 'amp =', amp
#endif
	end

