* born.F
* form factors for {Men(offshell Z)} -> {e+, e-, gamma} = i1 -> i2 + i3 + i4

	subroutine decayZ2lA(amp, i1, i2, i3, i4)
	use basic_funcs_FC6
	implicit none
	double complex amp
	integer i1, i2, i3, i4

	double complex Fx7, Fx2, Fx3, Fx8, Fx5, Fx6, Fx1, Fx4, Paix3, Paix2
	double complex Paix6, Paix5, Paix1, Paix4, Abx1, Abx3, Abx2
	double complex Abx4, Abx5, Abx7, Abx6, Abx8, AbxSum3, AbxSum1
	double complex AbxSum4, AbxSum2

#include "vars.h"
C checking
#ifdef CHECK_INPUT
        print*,"check 1 =",Pair(k(i1),k(i1))/MZ2
        print*,"check 0 =",Pair(k(i2),k(i2))/MZ2
        print*,"check 0 =",Pair(k(i3),k(i3))/MZ2
        print*,"check 0 =",Pair(k(i4),k(i4))/MZ2
        print*,"check 1 =",2*(Pair(k(i2),k(i3))
     & +Pair(k(i2),k(i4))
     & +Pair(k(i3),k(i4)))/MZ2
#endif

C angle
	Paix1 = Pair(k(i1),k(i2))

	Paix4 = Pair(k(i1),k(i3))

C abbr0_hel
	Fx7 = Hel(i2)*SxS(DottedSpinor(i3,1,7),VxS(e(i1),Spinor(i2,-1,7)))

	Fx2 = Hel(i2)*SxS(DottedSpinor(i3,1,7),VxS(ec(i4),Spinor(i2,-1,7)))

	Fx3 = Hel(i2)*SxS(DottedSpinor(i3,1,7),VxS(k(i1),Spinor(i2,-1,7)))

	Fx8 = -(Hel(i2)*SxS(DottedSpinor(i3,1,6),BxS(e(i1),Spinor(i2,-1,6))))

        Fx5 = -(Hel(i2)*SxS(DottedSpinor(i3,1,6),
     &      BxS(ec(i4),Spinor(i2,-1,6))))

	Fx6 = -(Hel(i2)*SxS(DottedSpinor(i3,1,6),BxS(k(i1),Spinor(i2,-1,6))))

	Paix3 = Pair(e(i1),ec(i4))

	Paix2 = Pair(e(i1),k(i2))

	Paix6 = Pair(e(i1),k(i3))

	Paix5 = Pair(ec(i4),k(i1))

        Fx1 = Hel(i2)*SxS(DottedSpinor(i3,1,7),
     &    VxS(e(i1),BxS(ec(i4),VxS(k(i1),Spinor(i2,-1,7)))))

        Fx4 = -(Hel(i2)*SxS(DottedSpinor(i3,1,6),
     &      BxS(e(i1),VxS(ec(i4),BxS(k(i1),Spinor(i2,-1,6))))))

	Abx1 = Fx2*Paix2

	Abx3 = Fx5*Paix2

	Abx2 = Fx3*Paix3

	Abx4 = Fx6*Paix3

	Abx5 = Fx7*Paix5

	Abx7 = Fx8*Paix5

	Abx6 = Fx2*Paix6

	Abx8 = Fx5*Paix6

	AbxSum3 = Abx5 - Abx6 - 1/2.D0*Fx1

	AbxSum1 = -2*Abx1 - 2*Abx2 + Fx1

	AbxSum4 = Abx7 - Abx8 - 1/2.D0*Fx4

	AbxSum2 = -2*Abx3 - 2*Abx4 + Fx4

C the amplitude
        amp =  
     &   (2*Alfa*Pi)/CW*(((1 - 2*SW2)*AbxSum1 - 2*SW2*AbxSum2)/
     &        (MZ2 - 2*Paix1) + 
     &       (2*((-1 + 2*SW2)*AbxSum3 + 2*SW2*AbxSum4))/
     &        (MZ2 - 2*Paix4))/SW

#ifdef DEBUG
	DEB 'amp =', amp
#endif
	end

