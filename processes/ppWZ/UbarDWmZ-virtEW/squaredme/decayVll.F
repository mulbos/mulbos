* form factors for {Mmm(offshell V)} -> {l+, l-} = i1 -> i2 + i3

	subroutine decayVll(ampZ, ampA, i1, i2, i3)
	use basic_funcs_FC6
	implicit none
	double complex ampZ, ampA, Fx1, Fx2
	integer i1, i2, i3

#include "vars.h"

	Fx1 = Hel(i2)*SxS(DottedSpinor(i3,1,7),VxS(e(i1),Spinor(i2,-1,7)))

	Fx2 = -(Hel(i2)*SxS(DottedSpinor(i3,1,6),BxS(e(i1),Spinor(i2,-1,6))))

        ampZ = 1/(2.D0*CW)*(EL*((-1 + 2*SW2)*Fx1 + 2*SW2*Fx2))/SW

	ampA = EL*(Fx1 + Fx2)

#ifdef DEBUG
	DEB 'decayWll: vec =', vec(1,1,1,i2), Hel(i2)
	DEB 'decayWll: amp =', amp
c	stop
#endif

	end
