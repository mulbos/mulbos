* form factors for {Mmm(offshell V)} -> {l+, l-} = i1 -> i2 + i3

	subroutine decayZll_NLO(ampZ, i1, i2, i3)
	use global_virt
	use basic_funcs_FC6
	implicit none
	double complex ampZ(2), Fx1, Fx2
	integer i1, i2, i3

#include "vars.h"

	Fx1 = Hel(i2)*SxS(DottedSpinor(i3,1,7),VxS(e(i1),Spinor(i2,-1,7)))

	Fx2 = -(Hel(i2)*SxS(DottedSpinor(i3,1,6),BxS(e(i1),Spinor(i2,-1,6))))

        ampZ(1) = 1/(2.D0*CW)*(EL*((-1 + 2*SW2)*Fx1 + 2*SW2*Fx2))/SW
        ampZ(2) = 0
C counter term
        ampZ(2) = ampZ(2) + 
     &   1/(4.D0*CW)*(EL*(2*Fx2*SW*DCONJG(dZfR1(2,2,2)) + 
     &        (-1 + 2*SW2)*(Fx1*DCONJG(dZfL1(2,2,2)))/SW + 
     &        2/CW2*(Fx2*(2*dSW1 + 
     &             CW2*(2*dZe1*SW + dZZZ1*SW + SW*dZfR1(2,2,2) + 
     &                CW*dZAZ1))) + 
     &        1/(CW2*SW2)*(Fx1*
     &            (2*SW2*(dSW1*SW) + 
     &              CW2*(2*dSW1*SW + 
     &                 SW2*
     &                  (-dZfL1(2,2,2) + (-2 + 4*SW2)*dZe1 + 
     &                    (-1 + 2*SW2)*dZZZ1 + 2*CW*(dZAZ1*SW) + 
     &                    2*SW2*dZfL1(2,2,2)))))/SW))

C loop
        ampZ(2) = ampZ(2) + 
     &   Alfa/(32.D0*CW*CW2*Pi*SW2)*
     &    (EL*(Fx1*xrat - 2*Fx1*Bval(bb0,iinz3) + 
     &         4*Fx1*Cval(cc00,iinz7) - 2*CW2*(Fx1*xrat) + 
     &         4*CW2**2*(Fx1*xrat) - 6*SW2*(Fx1*xrat) + 
     &         4*CW2*SW2*(Fx1*xrat) + 12*SW2**2*(Fx1*xrat) - 
     &         8*CW2*SW2**2*(Fx1*xrat) - 8*SW2**3*(Fx1*xrat) - 
     &         8*CW2*SW2**2*(Fx2*xrat) - 8*SW2**3*(Fx2*xrat) + 
     &         4*(1 - 2*CW2)*CW2*(Fx1*Bval(bb0,iinz2)) + 
     &         12*SW2*(Fx1*Bval(bb0,iinz3)) - 
     &         24*SW2**2*(Fx1*Bval(bb0,iinz3)) + 
     &         16*SW2**3*(Fx1*Bval(bb0,iinz3)) + 
     &         16*SW2**3*(Fx2*Bval(bb0,iinz3)) - 
     &         8*CW2**2*MW2*(Fx1*Cval(cc0,iinz4)) + 
     &         8*CW2*MZ2*SW2*(Fx1*Cval(cc0,iinz5)) - 
     &         16*CW2*MZ2*SW2**2*(Fx1*Cval(cc0,iinz5)) - 
     &         16*CW2*MZ2*SW2**2*(Fx2*Cval(cc0,iinz5)) - 
     &         16*CW2**2*(Fx1*Cval(cc00,iinz4)) + 
     &         16*CW2*SW2*(Fx1*Cval(cc00,iinz5)) - 
     &         32*CW2*SW2**2*(Fx1*Cval(cc00,iinz5)) - 
     &         32*CW2*SW2**2*(Fx2*Cval(cc00,iinz5)) - 
     &         8*CW2*(Fx1*Cval(cc00,iinz6)) - 
     &         24*SW2*(Fx1*Cval(cc00,iinz7)) + 
     &         48*SW2**2*(Fx1*Cval(cc00,iinz7)) - 
     &         32*SW2**3*(Fx1*Cval(cc00,iinz7)) - 
     &         32*SW2**3*(Fx2*Cval(cc00,iinz7)) + 
     &         8*CW2**2*MZ2*(Fx1*Cval(cc1,iinz4)) + 
     &         8*CW2*MZ2*SW2*(Fx1*Cval(cc1,iinz5)) - 
     &         16*CW2*MZ2*SW2**2*(Fx1*Cval(cc1,iinz5)) - 
     &         16*CW2*MZ2*SW2**2*(Fx2*Cval(cc1,iinz5)) + 
     &         4*CW2*MZ2*(Fx1*Cval(cc1,iinz6)) - 
     &         2*MZ2*(Fx1*Cval(cc1,iinz7)) + 
     &         12*MZ2*SW2*(Fx1*Cval(cc1,iinz7)) - 
     &         24*MZ2*SW2**2*(Fx1*Cval(cc1,iinz7)) + 
     &         16*MZ2*SW2**3*(Fx1*Cval(cc1,iinz7)) + 
     &         16*MZ2*SW2**3*(Fx2*Cval(cc1,iinz7)) + 
     &         4*CW2*MZ2*(Fx1*Cval(cc2,iinz6)) - 
     &         2*MZ2*(Fx1*Cval(cc2,iinz7)) + 
     &         12*MZ2*SW2*(Fx1*Cval(cc2,iinz7)) - 
     &         24*MZ2*SW2**2*(Fx1*Cval(cc2,iinz7)) + 
     &         16*MZ2*SW2**3*(Fx1*Cval(cc2,iinz7)) + 
     &         16*MZ2*SW2**3*(Fx2*Cval(cc2,iinz7)) + 
     &         8*CW2*SW2*(Bval(bb0,iinz1)*
     &            ((-1 + 2*SW2)*Fx1 + 2*SW2*Fx2))))/SW

#ifdef DEBUG
	DEB 'decayWll: vec =', vec(1,1,1,i2), Hel(i2)
	DEB 'decayZll LO: amp(1) =', amp(1)
	DEB 'decayZll NLO: amp(2) =', amp(2)
c	stop
#endif

	end
