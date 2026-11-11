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
	subroutine SoftFactor(Asq,SoftFac)
	use global_vars
	use global_virt
	implicit none

#include "model.h"

	double precision Gij, Asq, SoftFac
	external Gij

	integer i, j, nlegs
        parameter (nlegs = 4)

	double precision charge(nlegs), r
	
      if (lNLO_EW) then	
       charge(1) = charge_qed(1)
       charge(2) = charge_qed(2)
       charge(3) = charge_qed(3) + charge_qed(4)
       charge(4) = charge_qed(5) + charge_qed(6)
      else
       charge(1) = charge_qcd(1)
       charge(2) = charge_qcd(2)
       charge(3) = 0D0
       charge(4) = 0D0
      endif
      
c      print*,"charges =",charge
	
	SoftFac = 0

	do i = 1, nlegs
	  do j = 1, nlegs
	    r = charge(i)*charge(j)
	    if( r .ne. 0 ) SoftFac = SoftFac +
     &        sign(1, 2 - i)*sign(1, 2 - j)*
     &        r*Gij(i, j)

c          print*,'i,j,Gij=',i,j,Gij(i,j),r

	  enddo
	enddo

      if (lNLO_EW) then
	SoftFac = -Alfa/(2*pi)*SoftFac
      else
	SoftFac = -2*Alfas/(3*pi)*SoftFac
      endif

	SoftFac = SoftFac*Asq
	
c	print*,"endpoint =",SoftFac,Asq

	end


************************************************************************
	double precision function Gij(i, j)
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
	integer i_off, j_off

! on-shell mom (oext) : q(1) q(2) > W(3)Z(4) > e(5) nu_e(6) mu+(7) mu-(8)
! off-shell mom (pext): q(1) q(2) > e(3) nu_e(4) mu+(5) mu-(6)

	if( i .eq. j ) then
          Gij = 0 
	  return
	endif

* convention: i = emitter, j = spectator ***
* Endpoint terms for: massless incoming particles and 
* massive outgoing particles.
* For massless outgoing particles: formulae are simpler
* For massive incoming particles: can be tricky  

!*** common kinematics ***
        if ( i .gt. 2 ) then
	mi = mass_VV(i-2)
        else
        mi = MU
        endif
        if ( j .gt. 2 ) then
	mj = mass_VV(j-2)
        else
        mj = MU
        endif

	mi2 = mi*mi
	mj2 = mj*mj

***[initial - initial
! Note: Pij2_off = Pij2_OnShell here

	if (i .le. 2 .and. j .le. 2) then 
        i_off = i
        j_off = j

	if ( lambda2 .le. 0 ) then
***[for DimReg
! off-shell momenta = on-shell, 
! as explained in Denner ZZ 2021 paper
       Pij2 = pext(i_off,1)*pext(j_off,1)-pext(i_off,2)*pext(j_off,2) 
     & - pext(i_off,3)*pext(j_off,3) - pext(i_off,4)*pext(j_off,4)
       Pij2 = Pij2*2 ! (pi + pj)^2

	if ( lambda2 .eq. -2 ) then
	 Gij = 1D0
	elseif ( lambda2 .eq. -1 ) then
	 Gij = 3/2D0 + log(mu_ir2/Pij2)
	else
	 Gij = 3/2D0*log(mu_ir2/Pij2) + 1/2D0*log(mu_ir2/Pij2)**2 
	 Gij = Gij + 5 - pi**2/2
	endif
***]for DimReg
	endif

	if ( lambda2 .gt. 0 ) then
***[for MasReg	
! off-shell momenta, 
! as explained in Denner ZZ 2021 paper
       Pij2 = pext(i_off,1)*pext(j_off,1)-pext(i_off,2)*pext(j_off,2) 
     & - pext(i_off,3)*pext(j_off,3) - pext(i_off,4)*pext(j_off,4)
       Pij2 = Pij2*2 + mi2 + mj2 ! (pi + pj)^2

	Gij=-pi**2/3D0+2D0
	Gij = (log(mi2/Pij2)+1D0)*log(lambda2/Pij2) 
     &  -0.5D0*(log(mi2/Pij2)**2-log(mi2/Pij2))+Gij	
***]for MasReg	
        endif

	return
	endif 
***]initial - initial

***[final (i) - initial (j)
! Note: Pij2_off .ne. Pij2_OnShell here, 
! must use the on-shell one to match the Dipole subtraction 
! in the (n+1) contribution

	if (i .ge. 3 .and. j .le. 2) then 

	if ( lambda2 .le. 0 ) then
***[for DimReg
!	Pij2 = 2*DBLE(Pair(k(i), k(j)))
      Pij2 = oext(i,1)*oext(j,1)-oext(i,2)*oext(j,2)-oext(i,3)*oext(j,3)-oext(i,4)*oext(j,4)
      Pij2 = 2*Pij2
	Qij2 = Pij2 + mi2 + mj2

* Ref3. Eq. (6.52)
	if ( lambda2 .eq. -2 ) then
	 if ( mi2 .gt. 0 ) then
         Gij = 1D0/2
	 else
	 print*,'Error in softphoton_sub.F: special case FIM0_D not done'
	 stop
	 endif
	 return
	elseif ( lambda2 .eq. -1 ) then
	 if ( mi2 .gt. 0 ) then
* Ref3. Eq. (6.20)
         vv(-1) = log(mi2/Pij2)/2
	 gg(-1) = 1D0
	 Gij = vv(-1) + log(mu_ir2/Pij2)/2 + gg(-1)

c	 print*,'FI =', i, j, Gij

	 else
	 print*,'Error in softphoton_sub.F: special case FIM0_D not done'
	 stop
	 endif
	 return
	else
	 if ( mi2 .gt. 0 ) then
* Ref3. Eq. (6.20)
	 vv(0) = -log(mi2/Pij2)**2/4 - pi**2/12 
     &   - log(mi2/Pij2)*log(Pij2/Qij2)/2 - 
     &   log(mi2/Qij2)*log(Pij2/Qij2)/2

* Ref3. Eq. (6.22), the non-singular part
	 vv(0) = vv(0) + 3D0/2*log(Pij2/Qij2) + pi**2/6 - 
     &   Li2(Pij2/Qij2) - 2*log(Pij2/Qij2) - mi2/Pij2*log(mi2/Qij2)

	 Gij = log(mu_ir2/Pij2)**2/4 + 
     &   log(mu_ir2/Pij2)*log(mi2/Pij2)/2 + vv(0) - pi**2/3
	 gg(0) = log(mi2/mu_ir2)/2 - 2

	 Gij = Gij + gg(0) + log(mu_ir2/Pij2)*3/2 + 3D0/2 
     &  + 7D0/2 - pi**2/6
	 else
	 print*,'Error in softphoton_sub.F: special case FIM0_D not done'
	 stop
	 endif
	 return
	endif
***]for DimReg
	endif

	if ( lambda2 .gt. 0 ) then
***[for MasReg
!	Pij2 = -2*DBLE(Pair(k(i), k(j)))+mi2		
      Pij2 = oext(i,1)*oext(j,1)-oext(i,2)*oext(j,2)-oext(i,3)*oext(j,3)-oext(i,4)*oext(j,4)
      Pij2 = -2*Pij2 + mi2
	
        Gij = log(mi2/lambda2)*log(2-Pij2/mi2)
     &	+2*log(lambda*mi/(mi2-Pij2))-2*Li2(Pij2/(Pij2-2*mi2))
     &  +0.5D0*log(2-Pij2/mi2)**2+(Pij2-mi2)**2/(2*Pij2**2)
     &  *log(1-Pij2/mi2)-pi**2/6D0+1.5D0+mi2/(2*Pij2)
***]for MasReg
	endif

        return   
        endif 
***]final (i) - initial (j)
	
***[initial (i) - final (j)
! Note: Pij2_off .ne. Pij2_OnShell here, 
! must use the on-shell one to match the Dipole subtraction 
! in the (n+1) contribution

	if (i .le. 2 .and. j .ge. 3) then 
	if ( lambda2 .le. 0 ) then
***[for DimReg
!	Pij2 = 2*DBLE(Pair(k(i), k(j)))
      Pij2 = oext(i,1)*oext(j,1)-oext(i,2)*oext(j,2)-oext(i,3)*oext(j,3)-oext(i,4)*oext(j,4)
      Pij2 = 2*Pij2
	Qij2 = Pij2 + mi2 + mj2
	Qij = sqrt(Qij2)

* Ref3. Eq. (6.52)
	if ( lambda2 .eq. -2 ) then
	 if ( mj2 .gt. 0 ) then
         Gij = 1D0/2
	 else
	 print*,'Error in softphoton_sub.F: special case IFM0_D not done'
	 stop
	 endif
	 return
	elseif ( lambda2 .eq. -1 ) then
	 if ( mj2 .gt. 0 ) then
* Ref3. Eq. (6.20)
         vv(-1) = log(mj2/Pij2)/2
	 gg(-1) = 3D0/2
	 Gij = vv(-1) + log(mu_ir2/Pij2)/2 + gg(-1)

c	 print*,'IF =', i, j, Gij

	 else
	 print*,'Error in softphoton_sub.F: special case IFM0_D not done'
	 stop
	 endif
	 return
	else
	 if ( mj2 .gt. 0 ) then
* Ref3. Eq. (6.20)
	 vv(0) = -log(mj2/Pij2)**2/4 - pi**2/12 
     &   - log(mj2/Pij2)*log(Pij2/Qij2)/2 - 
     &   log(mj2/Qij2)*log(Pij2/Qij2)/2

* Ref3. Eq. (6.23), the non-singular part
	 vv(0) = vv(0) + 3D0/2*(log(Pij2/Qij2) - 2*log((Qij - mj)/Qij) 
     &   - 2*mj/(mj + Qij)) + pi**2/6 - Li2(Pij2/Qij2)

	 Gij = log(mu_ir2/Pij2)**2/4 + 
     &   log(mu_ir2/Pij2)*log(mj2/Pij2)/2 + vv(0) - pi**2/3
	 gg(0) = 0
	 Gij = Gij + gg(0) + log(mu_ir2/Pij2)*3/2 + 3D0/2 
     &  + 7D0/2 - pi**2/6
	 else
	 print*,'Error in softphoton_sub.F: special case IFM0_D not done'
	 stop
	 endif
	 return
	endif
***]for DimReg
	endif

	if ( lambda2 .gt. 0 ) then
***[for MasReg
!	Pij2 = -2*DBLE(Pair(k(i), k(j)))+mj2	
      Pij2 = oext(i,1)*oext(j,1)-oext(i,2)*oext(j,2)-oext(i,3)*oext(j,3)-oext(i,4)*oext(j,4)
      Pij2 = -2*Pij2 + mj2
	
        Gij = log(lambda2/mi2)*(log(mi2*(2*mj2-Pij2)/(mj2-Pij2)**2)+1)
     &  +2*Li2(Pij2/(2*mj2-Pij2))-2*Li2(mj2/(2*mj2-Pij2))
     &  +2*log(mj2*mi2/(mj2-Pij2)/(2*mj2-Pij2))
     &  *log((2*mj2-Pij2)/(mj2-Pij2))+0.5D0*log(mi2/(2*mj2-Pij2))**2
     &  +1.5D0*log(mi2/(mj2-Pij2))+mj2*(mj2-4*Pij2)/(2*Pij2**2)
     &  *log(1-Pij2/mj2)+pi**2/3D0-1D0+mj2/(2*Pij2) 	
***]for MasReg
	endif
	
	return 
        endif
***]initial (i) - final (j)
	
***[final - final
	if (i .ge. 3 .and. j .ge. 3) then
	if ( lambda2 .le. 0 ) then
#ifdef NINH_CATANI
***[for DimReg Catani
!	Pij2 = 2*DBLE(Pair(k(i), k(j)))
      Pij2 = oext(i,1)*oext(j,1)-oext(i,2)*oext(j,2)-oext(i,3)*oext(j,3)-oext(i,4)*oext(j,4)
      Pij2 = 2*Pij2
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
         Gij = 0
	 else
	 print*,'Error in softphoton_sub.F: special case FFM0_D not done'
	 stop
	 endif
	 return
	elseif ( lambda2 .eq. -1 ) then
	 if ( mi2 .gt. 0 .and. mj2 .gt. 0 ) then
* Ref3. Eq. (6.20)
         vv(-1) = log(rho)/vij
	 gg(-1) = 1D0
	 Gij = vv(-1) + gg(-1)

c	 print*,'FF =', i, j, Gij, rho, vij

	 else
	 print*,'Error in softphoton_sub.F: special case FFM0_D not done'
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

	 Gij = log(mu_ir2/Pij2)*log(rho)/vij + vv(0) - pi**2/3
	 gg(0) = log(mi2/mu_ir2)/2 - 2
	 Gij = Gij + gg(0) + log(mu_ir2/Pij2)*3/2 + 3D0/2 
     &  + 7D0/2 - pi**2/6
	 else
	 print*,'Error in softphoton_sub.F: special case FFM0_D not done'
	 stop
	 endif
	 return
	endif
***]for DimReg Catani
#endif

#ifdef NINH_DITTMAIER
***[for DimReg Dittmaier
!	Pij2 = 2*DBLE(Pair(k(i), k(j)))	
      Pij2 = oext(i,1)*oext(j,1)-oext(i,2)*oext(j,2)-oext(i,3)*oext(j,3)-oext(i,4)*oext(j,4)
      Pij2 = 2*Pij2
	lamij = Pij2**2 - 4*mi2*mj2
        rlamij = sqrt(lamij)
	
	a1 = (Pij2 + 2*mi2 - rlamij)/(Pij2 + 2*mi2 + rlamij)
	a2 = (Pij2 - rlamij)/(Pij2 + rlamij)	
	a3 = mi/(sqrt(Pij2+mi2+mj2) - mj)

	if ( lambda2 .eq. -2 ) then
	 if ( mi2 .gt. 0 .and. mj2 .gt. 0 ) then
         Gij = 0
	 else
	 print*,'Error in softphoton_sub.F: special case FFM0_D not done'
	 stop
	 endif
	 return
	elseif ( lambda2 .eq. -1 ) then
	 if ( mi2 .gt. 0 .and. mj2 .gt. 0 ) then
* Ref3. Eq. (4.10)
         vv(-1) = Pij2/rlamij*log(a1)
	 gg(-1) = 1D0
	 Gij = vv(-1) + gg(-1)

c	 print*,'FF =', i, j, Gij, rho, vij

	 else
	 print*,'Error in softphoton_sub.F: special case FFM0_D not done'
	 stop
	 endif
	 return
	else
	 if ( mi2 .gt. 0 .and. mj2 .gt. 0 ) then
	 Gij = 1.5D0+log(mu_ir2*a3**3/mi2)-2*log(1D0-a3**2)+a3**2/2D0
     &   +(Pij2/rlamij)*(log(a1)*log(mu_ir2*mj2/(lamij*a2))
     &   + 2*Li2(a1)+4*Li2(-sqrt(a2/a1))-4*Li2(-sqrt(a1*a2))
     &   + 0.5D0*log(a1)**2-pi**2/3D0)
	 else
	 print*,'Error in softphoton_sub.F: special case FFM0_D not done'
	 stop
	 endif
	return
	endif
***]for DimReg Dittmaier
#endif

	endif

	if ( lambda2 .gt. 0 ) then
***[for MasReg
!	Pij2 = 2*DBLE(Pair(k(i), k(j)))	
      Pij2 = oext(i,1)*oext(j,1)-oext(i,2)*oext(j,2)-oext(i,3)*oext(j,3)-oext(i,4)*oext(j,4)
      Pij2 = 2*Pij2
	lamij = Pij2**2 - 4*mi2*mj2
        rlamij = sqrt(lamij)
	
	a1 = (Pij2 + 2*mi2 - rlamij)/(Pij2 + 2*mi2 + rlamij)
	a2 = (Pij2 - rlamij)/(Pij2 + rlamij)	
	a3 = mi/(sqrt(Pij2+mi2+mj2) - mj)			
	 
	 Gij = 1.5D0+log(lambda2*a3**3/mi2)-2*log(1D0-a3**2)+a3**2/2D0
     &   +(Pij2/rlamij)*(log(a1)*log(lambda2*mj2/(lamij*a2))
     &   + 2*Li2(a1)+4*Li2(-sqrt(a2/a1))-4*Li2(-sqrt(a1*a2))
     &   + 0.5D0*log(a1)**2-pi**2/3D0)
***]for MasReg
	endif

        return
	endif
***]final - final

	end

