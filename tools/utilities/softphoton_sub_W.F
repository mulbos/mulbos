* endpoint.F
* routines for the treatment of IR-sensitive enpoint contribution
* in the framework of dipole subtraction scheme
* Ref0. Dittmaier et al 1507.04676

! Ninh implementation for the Endpoint

	subroutine SoftFactor_decayW(Asq,SoftFac)
	use global_vars
	use global_virt
	implicit none

#include "model.h"

	double precision Wsf, Asq, SoftFac
	external Wsf

	SoftFac = -Wsf()

      if (lNLO_EW) then	
	SoftFac = -Alfa/(2*pi)*SoftFac
      else
	SoftFac = -2*Alfas/(3*pi)*SoftFac
      endif

	SoftFac = SoftFac*Asq
	end
************************************************************************
	double precision function Wsf()
	use global_vars
	use global_virt
	use global_mom_xpdf
	implicit none
	integer i, a

#include "model.h"

	double precision mi, ma
	double precision mi2, ma2, half
        parameter (half = 1D0/2)

	double precision log_ia
        double precision log_ga, temp

	if ( lambda2 .le. 0 ) then
***[for DimReg
         print*,"Error in softphoton_sub_W.F: final-init"
         print*,"DimReg for massive initial-state particle not done"
         stop
***]for DimReg
	endif

* convention: i = final, a = initial ***

*** common kinematics ***
        mi = ME
        ma = MW

	mi2 = mi*mi
	ma2 = ma*ma

c        print*,"mgam, mi, ma =",lambda, mi, ma

c we have P_ia^2 = m_nu^2 = 0, 
C Qia^2 = -ma2 - mi2, gamma = 1

        log_ia = log(mi2/ma2)
        log_ga = log(lambda2/ma2)

c Ref0. eq. (3.23)
        Wsf = log_ga*log_ia + 2*log_ga - half*log_ia**2 
     & + half*log_ia - pi**2/3 + 1.5D0

c        print*,"Wsf sum",Wsf

	end
