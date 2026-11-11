* endpoint.F
* routines for the treatment of IR-sensitive enpoint contribution
* in the framework of dipole subtraction scheme
* Ref0. Dittmaier hep-ph/9904440v2
* Ref1. Catani and Seymour hep-ph/9605323
* Ref2. Beenakker, Dittmaier et al hep-ph/0211352v1
* Ref3. Catani, Dittmaier et al hep-ph/0201036v1

* Option NINH_DITTMAIER is recommended, to match 
* the (N+1) contribution, same for EW and QCD. 
#define NINH_DITTMAIER
* This option NINH_CATANI is used if the (N+1) contribution 
* use the Catani mapping as done in Autodipole
cc#define NINH_CATANI

! Ninh implementation for the Endpoint

	subroutine SoftFactor_decayZ(Asq,SoftFac)
	use global_vars
	use global_virt
	implicit none

#include "model.h"

	double precision Zsf, Asq, SoftFac
	external Zsf

	SoftFac = -Zsf(7,8)-Zsf(8,7)

      if (lNLO_EW) then	
	SoftFac = -Alfa/(2*pi)*SoftFac
      else ! NLOQCD
      	SoftFac = -2*Alfas/(3*pi)*SoftFac
      endif

	SoftFac = SoftFac*Asq
	end
************************************************************************
	double precision function Zsf(i,j)
	use global_vars
	use global_virt
	use global_mom_xpdf
	implicit none
	integer i, j

#include "model.h"
!!![matching loopints.h
	double complex Li2
	external Li2
!!!]matching loopints.h

	double precision mi, mj, Pij2, a1, a2, a3
	double precision mi2, mj2, lamij, rlamij

	double precision LambdaX, vij, rho, Qij2, Qij, mui2, muj2

	double precision rhoi, rhoj, vv(-2:0), gg(-2:0)
	double precision log_fer, log_pho 
	integer i_off, j_off

* convention: i = emitter, j = spectator ***
* Endpoint terms for: massless incoming particles and 
* massive outgoing particles.
* For massless outgoing particles: formulae are simpler
* For massive incoming particles: can be tricky  

*** common kinematics ***
        mi = MM ! muon mass
        mj = mi

	mi2 = mi*mi
	mj2 = mj*mj
! on-shell mom (oext) : q(1) q(2) > W(3)Z(4) > e(5) nu_e(6) mu+(7) mu-(8)
! off-shell mom (pext): q(1) q(2) > e(3) nu_e(4) mu+(5) mu-(6)
        i_off = i - 2
        j_off = j - 2
	
***[final - final

! Other choice: using Pij2 = ma2_off (ma2=MZ2)
! may be OK, but we have to use PiK_off = yij*ma2_off/2 in 
! the corresponding subtraction term in the (n+1) part.
! Off-shell MZ2 choice (see FuncDPsing_decayZ.F):
!       Pij2 = pext(i_off,1)*pext(j_off,1)-pext(i_off,2)*pext(j_off,2) 
!     & - pext(i_off,3)*pext(j_off,3) - pext(i_off,4)*pext(j_off,4)
!       Pij2 = Pij2*2 ! (pi + pj)^2

! On-shell Pij2 = MZ2 choice (default):
!	Pij2 = 2*DBLE(Pair(k(i), k(j)))	! on-shell momenta (Ninh paper WZ 2018)
      Pij2 = oext(i,1)*oext(j,1)-oext(i,2)*oext(j,2)-oext(i,3)*oext(j,3)-oext(i,4)*oext(j,4)
      Pij2 = 2*Pij2

!        print*,"Pij2 =",Pij2
!        print*,"MZ2  =",MZ2
!        stop

!!![dimensional regularization
	if ( lambda2 .le. 0 ) then
#ifdef NINH_CATANI
***[for DimReg Catani
	Qij2 = Pij2 + mi2 + mj2
	Qij = sqrt(Qij2)
	mui2 = mi2/Qij2
	muj2 = mj2/Qij2

* Ref3. Eqs. (6.20), (5.30)
	vij = sqrt(LambdaX(Qij2, mi2, mj2))/Pij2
	rho = (1 - vij)/(1 + vij)
	rho = sqrt(rho)
	rhoi = (1 - vij + 2*mui2/(1 - mui2 - muj2))/
     &  (1 + vij + 2*mui2/(1 - mui2 - muj2))
	rhoi = sqrt(rhoi)
	rhoj = (1 - vij + 2*muj2/(1 - mui2 - muj2))/
     &  (1 + vij + 2*muj2/(1 - mui2 - muj2))
	rhoj = sqrt(rhoj)

	if ( lambda2 .eq. -2 ) then
	 if ( mi2 .gt. 0 .and. mj2 .gt. 0 ) then
         Zsf = 0
	 else
	 print*,'Error in softphoton_sub_Z.F: special case FFM0_D not done'
	 stop
	 endif
	 return
	elseif ( lambda2 .eq. -1 ) then
	 if ( mi2 .gt. 0 .and. mj2 .gt. 0 ) then
* Ref3. Eq. (6.20)
         vv(-1) = log(rho)/vij
	 gg(-1) = 1D0
	 Zsf = vv(-1) + gg(-1)

c	 print*,'FF =', i, j, Zsf, rho, vij

	 else
	 print*,'Error in softphoton_sub_Z.F: special case FFM0_D not done'
	 stop
	 endif
	 return
	else
	 if ( mi2 .gt. 0 .and. mj2 .gt. 0 ) then
* Ref3. Eq. (6.20)
	 vv(0) = (-log(rhoi**2)**2/4 - log(rhoj**2)**2/4 - pi**2/6)/vij 
     &  + log(rho)*log(Qij2/Pij2)/vij
* Ref3. Eq. (6.21), the non-singular part
	 vv(0) = vv(0) + 3D0/2*log(Pij2/Qij2) + 
     &   (log(rho**2)*log(1 + rho**2) + 2*Li2(rho**2) - 
     &    Li2(1 - rhoi**2) - Li2(1 - rhoj**2) - pi**2/6)/vij + 
     &   log((Qij - mj)/Qij) - 2*log(((Qij-mj)**2-mi2)/Qij2) - 
     &   2*mi2/Pij2*log(mi/(Qij-mj)) - mj/(Qij - mj) + 
     &   2*mj*(2*mj - Qij)/Pij2 + pi**2/2

	 Zsf = log(mu_ir2/Pij2)*log(rho)/vij + vv(0) - pi**2/3
	 gg(0) = log(mi2/mu_ir2)/2 - 2
	 Zsf = Zsf + gg(0) + log(mu_ir2/Pij2)*3/2 + 3D0/2 
     &  + 7D0/2 - pi**2/6
	 else
	 print*,'Error in softphoton_sub_Z.F: special case FFM0_D not done'
	 stop
	 endif
	 return
	endif
***]for DimReg Catani
#endif

#ifdef NINH_DITTMAIER
***[for DimReg Dittmaier
	lamij = Pij2**2 - 4*mi2*mj2
        rlamij = sqrt(lamij)

	
	a1 = (Pij2 + 2*mi2 - rlamij)/(Pij2 + 2*mi2 + rlamij)
	a2 = (Pij2 - rlamij)/(Pij2 + rlamij)	
	a3 = mi/(sqrt(Pij2+mi2+mj2) - mj)

	if ( lambda2 .eq. -2 ) then
	 if ( mi2 .gt. 0 .and. mj2 .gt. 0 ) then
         Zsf = 0
	 else
	 print*,'Error in softphoton_sub_Z.F: special case FFM0_D2 not done'
	 stop
	 endif
	 return
	elseif ( lambda2 .eq. -1 ) then
	 if ( mi2 .gt. 0 .and. mj2 .gt. 0 ) then
* Ref3. Eq. (4.10)
         vv(-1) = Pij2/rlamij*log(a1)
	 gg(-1) = 1D0
	 Zsf = vv(-1) + gg(-1)

c	 print*,'FF =', i, j, Zsf, rho, vij

	 else
	 print*,'Error in softphoton_sub_Z.F: special case FFM0_D1 not done'
	 stop
	 endif
	 return
	else
	 if ( mi2 .gt. 0 .and. mj2 .gt. 0 ) then
	 Zsf = 1.5D0+log(mu_ir2*a3**3/mi2)-2*log(1D0-a3**2)+a3**2/2D0
     &   +(Pij2/rlamij)*(log(a1)*log(mu_ir2*mj2/(lamij*a2))
     &   + 2*Li2(a1)+4*Li2(-sqrt(a2/a1))-4*Li2(-sqrt(a1*a2))
     &   + 0.5D0*log(a1)**2-pi**2/3D0)
	 else
	 print*,'Error in softphoton_sub_Z.F: special case FFM0_D0 not done'
	 print*,'lambda2, mi2, mj2 =',lambda2,mi2,mj2
	 stop
	 endif
	return
	endif
***]for DimReg Dittmaier
#endif

	endif
!!!]dimensional regularization

	if ( lambda2 .gt. 0 ) then
***[for MasReg
C Ref0 (Dittmaier 1999) eq. 3.7 
       log_fer = log(mi2/Pij2)
       log_pho = log(lambda2/Pij2)
       Zsf = log_fer*log_pho+log_pho-log_fer**2/2D0+log_fer/2D0
       Zsf = Zsf - pi**2/3D0 + 3D0/2
c          print*,"Mass0 Zsf =",Zsf,log_fer,mi
***]for MasReg
	endif
***]final - final

        return
	end
