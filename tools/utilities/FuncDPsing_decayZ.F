	double precision function FuncDPsing_decayZ(fac)
	use global_vars
	implicit none

#include "model.h"

	double precision GdZ
	external GdZ

	integer l1, l2
!        double precision MVint2(2) ! the off-shell invariant mass of the gauge bosons	
!	double precision qV(4)
	double precision Jacb, fac
	
C We need the off-shell gauge-boson invariant masses, 
C to pass to the reduced amplitudes. 
C Here we take the simplest choice that those invariant 
C masses are independent of dipoles (= emitter + spectator).
C Note: pext(6,4) are the true qq -> 4l + gamma momenta
!	MVint2(1) = pext(3,1)*pext(4,1) - pext(3,2)*pext(4,2) 
!     &  - pext(3,3)*pext(4,3) - pext(3,4)*pext(4,4)
!        MVint2(1) = MVint2(1)*2
!        do i = 1, 4
!         qV(i) = pext(5,i) + pext(6,i) + pext(7,i)
!        enddo
!	MVint2(2) = qV(1)**2 - qV(2)**2 - qV(3)**2 - qV(4)**2
!        print*,"FuncDPsing_decayZ: MVint2 =",MVint2

	Jacb = -1D0*fac

      if (lNLO_QCD) then
       !Jacb = -16*pi*Jacb/3D0
       print*,"error: hadronic decay of Z not implemented"
       stop
      else
       Jacb = -4*pi*Jacb
      endif

! order of particles: p(1)p(2) > V(3)V(4) > photon(5) leptons(6,7,8,9)
      if (lDECAYV1_NLO_EW) then
       l1 = 6
       l2 = 7
      else
       l1 = 8
       l2 = 9
      endif

       FuncDPsing_decayZ =
     & + Jacb*GdZ(l1, l2, Jacb)
     & + Jacb*GdZ(l2, l1, Jacb)

	end

************************************************************************
	double precision function GdZ(i, j, Jacb)
        use global_vars
        use global_mom_xpdf
	implicit none
	integer i, j, i1, nlegs
        parameter (nlegs = 5) ! for parent q1 q2 > V3 V4 g5

! proc_nlegs = 7 
! proc_nlegs_os = 9
        double precision Jacb
	double precision temp, temp1
	double precision PiK, PjK, PiPj		
	double precision xij, yij, zij	

        double precision tp(2,4) ! Z > l8 + l9 + phot, phot = oext(5)
        double precision tp_off(2,4) ! off shell version of tp()
        double precision PiPj_off,PiK_off,PjK_off,yij_off
	integer i_off, j_off, shift
!        double precision ma2 ! squared mass of the decaying particle

!        ma2 = MZ2 ! Z decay

      if (lDECAYV1_NLO_EW) then
       shift = 5
      else
       shift = 7
      endif

*** convention: i = emitter, j = spectator ***

! oext: q(1)q(2) > W(3)Z(4)gamma(5) > e(6)nu_e(7)mu+(8)mu-(9)
! pext: q(1)q(2) > e(3)nu_e(4)mu+(5)mu-(6)gamma(7), proc_nlegs = 7
        i_off = i-3
        j_off = j-3

***[ final(emitter) - final(spectator)
        PiK = oext(i,1)*oext(nlegs,1)-oext(i,2)*oext(nlegs,2)
     &  -oext(i,3)*oext(nlegs,3)-oext(i,4)*oext(nlegs,4)

        PjK = oext(j,1)*oext(nlegs,1)-oext(j,2)*oext(nlegs,2)
     &  -oext(j,3)*oext(nlegs,3)-oext(j,4)*oext(nlegs,4)
     
        PiPj = oext(i,1)*oext(j,1)-oext(i,2)*oext(j,2)
     &  -oext(i,3)*oext(j,3)-oext(i,4)*oext(j,4)	
	
c 
        yij = PiK/(PiK+PjK+PiPj)
	zij = PiPj/(PiPj+PjK) 
	

***[the off-shell correspondence
      PiK_off=pext(i_off,1)*pext(proc_nlegs,1)-pext(i_off,2)*pext(proc_nlegs,2)
     &  -pext(i_off,3)*pext(proc_nlegs,3)-pext(i_off,4)*pext(proc_nlegs,4)
      PjK_off=pext(j_off,1)*pext(proc_nlegs,1)-pext(j_off,2)*pext(proc_nlegs,2)
     &  -pext(j_off,3)*pext(proc_nlegs,3)-pext(j_off,4)*pext(proc_nlegs,4)
      PiPj_off=pext(i_off,1)*pext(j_off,1)-pext(i_off,2)*pext(j_off,2)
     &  -pext(i_off,3)*pext(j_off,3)-pext(i_off,4)*pext(j_off,4)	

        yij_off = PiK_off/(PiK_off+PjK_off+PiPj_off)
***]the off-shell correspondence

!        print*,"FuncDPsing yij =",yij,yij_off

! Note: yij_off = yij and zij_off = zij if the OS mapping of 
! Denner ZZ 2021 paper is used.

***[get the on-shell momenta
	do i1 = 1,4
	tp(j-shift,i1)=oext(j,i1)/(1-yij)
        tp(i-shift,i1)=oext(i,i1)+oext(nlegs,i1)-oext(j,i1)*yij/(1-yij)
	enddo
***]get the on-shell momenta
	 
***[get the off-shell momenta
! oext: q(1)q(2) > W(3)Z(4)gamma(5) > e(6)nu_e(7)mu+(8)mu-(9)
! pext: q(1)q(2) > e(3)nu_e(4)mu+(5)mu-(6)gamma(7)
	do i1 = 1,4
	tp_off(j-shift,i1)=pext(j_off,i1)/(1-yij_off)
        tp_off(i-shift,i1)=pext(i_off,i1)
     & +pext(proc_nlegs,i1)-pext(j_off,i1)*yij_off/(1-yij_off)
	enddo
***]get the off-shell momenta

! Note for GdZ: this PiK kinematic factor must be on-shell 
! to match the singular propagator in the on-shell 
! (n+1) amplitude. Note: PiK = yij*MZ2/2, 
! where yij = yij_off is the integration variable 
! when integrating over the photon phase space [dk]

! Other choice: using PiK_off = yij*ma2_off/2 
! may be OK (convergence behavior ?), but we have to use Pij2 = ma2_off in 
! the Endpoint to match this subtraction term.

! Note: PiK = yij*MZ2/2 (default choice)

      GdZ = (2D0/(1-zij*(1-yij))-1-zij)/PiK/(1-yij)

      call DoTildeAmpHis_decayZ_OffS(tp,tp_off,1,GdZ,Jacb)

        return
***] final(emitter) - final(spectator)

	end
