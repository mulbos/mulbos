	double precision function FuncDPsing_fin_rad(fac)
	use global_vars
	implicit none

#include "model.h"

	double precision Gfr
	external Gfr

	integer i, j, nlegs
        parameter (nlegs = 5) ! for the on-shell 2to3 process
!        double precision MVint2(2) ! the off-shell invariant mass of the gauge bosons
	double precision charge(nlegs), r
	
	double precision Jacb, fac

      if (lNLO_EW) then	
       charge(1) = charge_qed(1)
       charge(2) = charge_qed(2)
       charge(3) = charge_qed(3) + charge_qed(4)
       charge(4) = charge_qed(5) + charge_qed(6)
       charge(5) = 0D0
      else
       charge(1) = charge_qcd(1)
       charge(2) = charge_qcd(2)
       charge(3) = 0D0
       charge(4) = 0D0
       charge(5) = 0D0
      endif
      
!      print*,"charge =",charge
	
	FuncDPsing_fin_rad = 0

C We need the off-shell gauge-boson invariant masses, 
C to pass to the reduced amplitudes. 
C Here we take the simplest choice that those invariant 
C masses are independent of dipoles (= emitter + spectator).
C Note: pext(6,4) are the true qq -> 4l + gamma momenta
!	MVint2(1) = pext(3,1)*pext(4,1) - pext(3,2)*pext(4,2) 
!     &  - pext(3,3)*pext(4,3) - pext(3,4)*pext(4,4)
!        MVint2(1) = MVint2(1)*2
!	MVint2(2) = pext(5,1)*pext(6,1) - pext(5,2)*pext(6,2) 
!     &  - pext(5,3)*pext(6,3) - pext(5,4)*pext(6,4)
!        MVint2(2) = MVint2(2)*2

	do i = 1, nlegs
	  do j = 1, nlegs
	    r = charge(i)*charge(j)
	    if( r .ne. 0 ) then
	    Jacb = sign(1, 2 - i)*sign(1, 2 - j)*r  

      if (lNLO_QCD) then
            Jacb = -Alfas*16*pi*Jacb*fac/3D0
      else if (lNLO_EW) then
            Jacb = -Alfa*4*pi*Jacb*fac
      endif

	    FuncDPsing_fin_rad = FuncDPsing_fin_rad + Jacb*Gfr(i,j,Jacb)/fac

c	   print*,'Gfr =',i,j,Gfr(i,j,Jacb),fac

            endif   
	  enddo
	enddo

c        stop
	end

************************************************************************
	double precision function Gfr(i,j,Jacb)
        use global_vars
        use global_mom_xpdf
	implicit none
	double precision Jacb
	integer i, j, i1, j1, j2, nlegs
        parameter (nlegs = 5) ! for parent q1 q2 > V3 V4 g5

	double precision mi, mj, mi2, mj2, Pij2, BPij2, lamij, rlamij
	double precision temp, temp1, temp2	
	double precision PiK, PjK, PiPj		
	double precision xij, yij, zij	

        double precision Pij(4), TPij(4)
        double precision tq(proc_nlegs-1,4) ! off-shell tilde momenta
!        double precision MVint2(2) ! the off-shell invariant mass of the gauge bosons (N+1 kin)
!        double precision MVint2_x(2) ! the off-shell invariant mass of the gauge bosons (reduced kin)
        double precision qW(4), tqW(4)
        double precision pW(3:4,4), tpW(3:4,4) ! this is for the WW process
      double precision PiPj_off, PiK_off, PjK_off, xij_off
      double precision Pij_off(4), Pij2_off
      double precision BPij2_off, lamij_off, rlamij_off, yij_off
      double precision p(4,2), q(4,2), MW_off, mi2_off, mj2_off
      double precision dot_p2

	if( i .eq. j ) then
          Gfr = 0 
	  return
	endif

* Ninh checking
c	if( i .le. 2 .or. j .gt. 2 ) then
c          Gfr = 0 
c          return
c	endif


*** convention: i = emitter, j = spectator ***

*** common kinematics *** 
	temp1 = sqrt(oext(i,2)**2 + oext(i,3)**2 + oext(i,4)**2)
        mi2 = (oext(i,1)-temp1)*(oext(i,1)+temp1)

        temp2 = sqrt(oext(j,2)**2 + oext(j,3)**2 + oext(j,4)**2)
        mj2 = (oext(j,1)-temp2)*(oext(j,1)+temp2)

	if (abs(mi2) .lt. 1D-8*temp1**2) mi2 = 0
	if (abs(mj2) .lt. 1D-8*temp2**2) mj2 = 0

        mi = sqrt(mi2)
	mj = sqrt(mj2)

c        print*,'Gfr mass hello =',i,j,mi,mj,mi2,mj2
#define NINH_CHECK_INI_INI
#define NINH_CHECK_FIN_INI
#define NINH_CHECK_INI_FIN
#define NINH_CHECK_FIN_FIN

***[ initial(emitter) - initial(spectator)
#ifdef NINH_CHECK_INI_INI
	if (i .le. 2 .and. j .le. 2) then 

        if (mi2 + mj2 .gt. 0) then
	 print*,'massive initial case is not coded'
	 print*,'mi, mj =',i,j,mi,mj
	call PrintPSP()
        call PrintPSP_OS()
	 stop
	endif 

***[for DPA approximation
***[for NLO QCD and EW
! mapping method to re-produce Denner WZ 2020 results. 
! idea is to use exact momenta to calculate tilde kinematics.
! Convergence is excellent.

! for the radiation process:
! Note: pext(7,4) are the qq -> 4l + j off-shell momenta

        PiK = pext(i,1)*pext(proc_nlegs,1)-pext(i,2)*pext(proc_nlegs,2)
     &  -pext(i,3)*pext(proc_nlegs,3)-pext(i,4)*pext(proc_nlegs,4)

        PjK = pext(j,1)*pext(proc_nlegs,1)-pext(j,2)*pext(proc_nlegs,2)
     &  -pext(j,3)*pext(proc_nlegs,3)-pext(j,4)*pext(proc_nlegs,4)
     
        PiPj = pext(i,1)*pext(j,1)-pext(i,2)*pext(j,2)
     &  -pext(i,3)*pext(j,3)-pext(i,4)*pext(j,4)

	xij = (PiPj-PiK-PjK)/PiPj

! on-shell = off-shell momenta here, 
! because on-shell projection does not change 
! initial and photon momenta

        Gfr = (2D0/(1-xij)-1-xij)/(PiK*xij)

        do i1 = 1,4
	Pij(i1)=pext(i,i1)+pext(j,i1)-pext(proc_nlegs,i1)  
	enddo
!	Pij2 = Pij(1)**2-Pij(2)**2-Pij(3)**2-Pij(4)**2
	Pij2 = dot_p2(Pij(1))

***[get the reduced momenta of the qq > 4l + g exact kinematics 
	do i1 = 1,4
        tq(i,i1)=pext(i,i1)*xij
        tq(j,i1)=pext(j,i1)	
	TPij(i1)=pext(j,i1)+tq(i,i1)
	enddo

	temp = Pij(1)*TPij(1)-Pij(2)*TPij(2)
     &	-Pij(3)*TPij(3)-Pij(4)*TPij(4)

	do j1 = 3, proc_nlegs-1 
	temp1 = Pij(1)*pext(j1,1)-Pij(2)*pext(j1,2)
     &	-Pij(3)*pext(j1,3)-Pij(4)*pext(j1,4) 	 	
     
	temp2 = TPij(1)*pext(j1,1)-TPij(2)*pext(j1,2)
     &	-TPij(3)*pext(j1,3)-TPij(4)*pext(j1,4)
     
        do i1 = 1,4
	tq(j1,i1)=pext(j1,i1)
     &	-(Pij(i1)+TPij(i1))*(temp1+temp2)/(Pij2+temp)
     &  +2*TPij(i1)*temp1/Pij2
	enddo  	 	     
     
        enddo
***]get the reduced momenta of the qq > 4l + g exact kinematics 

!!![Ninh checking
c        print*,"FuncDPsing.F: exact reduced momenta"
c        print*,"emitter-spectator:",i,j
c        do i1 = 1, proc_nlegs-1
c        print*,"particle :",i1,tq(i1,:)
c        enddo
!!!]Ninh checking

       call DoTildeAmpHis_OffS(tq,1,Gfr,Jacb)
***]for NLO QCD and EW
***]for DPA approximation
	return
	endif
#endif
***] initial(emitter) - initial(spectator)

***[ final (i-emitter) - initial (j-spectator)
#ifdef NINH_CHECK_FIN_INI
	if (i .ge. 3 .and. j .le. 2) then 

!!![Lorentz inv, use on-shell kinematics
! the photon: oext(nlegs,:) = on-shell = off-shell
        PiK = oext(i,1)*oext(nlegs,1)-oext(i,2)*oext(nlegs,2)
     &  -oext(i,3)*oext(nlegs,3)-oext(i,4)*oext(nlegs,4)

        PjK = oext(j,1)*oext(nlegs,1)-oext(j,2)*oext(nlegs,2)
     &  -oext(j,3)*oext(nlegs,3)-oext(j,4)*oext(nlegs,4)
     
        PiPj = oext(i,1)*oext(j,1)-oext(i,2)*oext(j,2)
     &  -oext(i,3)*oext(j,3)-oext(i,4)*oext(j,4)

	xij = (PiPj-PiK+PjK)/(PiPj+PjK)
	zij = PiPj/(PiPj+PjK)

! Must use on-shell momenta here to cancel 
! the IR divergences in the on-shell qq > WZ gamma amplitude. 
! This means that the on-shell momenta must be used for the 
! corresponding terms (mi^2 and Pia^2) in the collinear 
! and endpoint contributions.
	 
        Gfr = (2D0/(2-xij-zij)-1-zij-mi2/PiK)/(PiK*xij)
!!!]Lorentz inv, use on-shell kinematics

***[the off-shell correspondence
	do i1 = 1,4
        qW(i1) = pext(3,i1) + pext(4,i1) ! off-shell emitter W
        enddo

        PiK_off = qW(1)*pext(proc_nlegs,1)-qW(2)*pext(proc_nlegs,2)
     &  -qW(3)*pext(proc_nlegs,3)-qW(4)*pext(proc_nlegs,4)

        PjK_off = pext(j,1)*pext(proc_nlegs,1)-pext(j,2)*pext(proc_nlegs,2)
     &  -pext(j,3)*pext(proc_nlegs,3)-pext(j,4)*pext(proc_nlegs,4)

        PiPj_off = qW(1)*pext(j,1)-qW(2)*pext(j,2)
     &  -qW(3)*pext(j,3)-qW(4)*pext(j,4)

	xij_off = (PiPj_off-PiK_off+PjK_off)/(PiPj_off+PjK_off)
***]the off-shell correspondence

***[get the off-shell 2to4 momenta
	do i1 = 1,4
! eq. (3.12) of Dittmaier 1999
! must use xij_off for tqW^2 = qW^2 condition
        tqW(i1)=qW(i1)+pext(proc_nlegs,i1)-(1-xij_off)*pext(j,i1) ! combine the photon with the W
	tq(j,i1)=xij_off*pext(j,i1) ! spectator
	enddo

! the remaining initial-state particle
	do j1 = 1, 2
	if (j1 .ne. j) then
	 do i1 = 1,4
	 tq(j1,i1)=pext(j1,i1)
	 enddo
	endif
	enddo

! the W decay leptons
      do i1 = 1,4
       p(i1,1) = pext(3,i1)
       p(i1,2) = pext(4,i1)
      enddo
      MW_off = dot_p2(tqW(1))
      MW_off = sqrt(MW_off) 
      call OSMap_decayV_LO(p(1,1),MW_off,tqW(1),q(1,1))
      do i1 = 1,4
        tq(3,i1) = q(i1,1)
        tq(4,i1) = q(i1,2)
      enddo

! the Z decay leptons
	do j1 = 5, 6
	 do i1 = 1,4
	 tq(j1,i1)=pext(j1,i1)
	 enddo
	enddo
***]get the off-shell 2to4 momenta
      call DoTildeAmpHis_OffS(tq,1,Gfr,Jacb)
        return   
        endif
#endif
***] final (i-emitter) - initial (j-spectator) 


***[ initial (i-emitter) - final (j-spectator)
#ifdef NINH_CHECK_INI_FIN
	if (i .le. 2 .and. j .ge. 3) then 

!!![Lorentz inv, use on-shell kinematics
! the photon: oext(nlegs,:) = on-shell = off-shell
        PiK = oext(i,1)*oext(nlegs,1)-oext(i,2)*oext(nlegs,2)
     &  -oext(i,3)*oext(nlegs,3)-oext(i,4)*oext(nlegs,4)

        PjK = oext(j,1)*oext(nlegs,1)-oext(j,2)*oext(nlegs,2)
     &  -oext(j,3)*oext(nlegs,3)-oext(j,4)*oext(nlegs,4)
     
        PiPj = oext(i,1)*oext(j,1)-oext(i,2)*oext(j,2)
     &  -oext(i,3)*oext(j,3)-oext(i,4)*oext(j,4)

	xij = (PiPj-PjK+PiK)/(PiPj+PiK)
	zij = PiPj/(PiPj+PiK)
	 
! must use the on-shell momenta here (same as above)
        Gfr = (2D0/(2-xij-zij)-1-xij)/(PiK*xij)
!!!]Lorentz inv, use on-shell kinematics

***[the off-shell correspondence
	do i1 = 1,4
        qW(i1) = pext(3,i1)+pext(4,i1) ! off-shell spectator W
        enddo

        PiK_off = pext(i,1)*pext(proc_nlegs,1)-pext(i,2)*pext(proc_nlegs,2)
     &  -pext(i,3)*pext(proc_nlegs,3)-pext(i,4)*pext(proc_nlegs,4)

        PjK_off = qW(1)*pext(proc_nlegs,1)-qW(2)*pext(proc_nlegs,2)
     &  -qW(3)*pext(proc_nlegs,3)-qW(4)*pext(proc_nlegs,4)

        PiPj_off = pext(i,1)*qW(1)-pext(i,2)*qW(2)
     &  -pext(i,3)*qW(3)-pext(i,4)*qW(4)

	xij_off = (PiPj_off-PjK_off+PiK_off)/(PiPj_off+PiK_off)
***]the off-shell correspondence

***[get the off-shell 2to4 momenta
	do i1 = 1,4
! eq. (3.12) of Dittmaier 1999
! must use xij_off for tqW^2 = qW^2 condition
	tq(i,i1)=xij_off*pext(i,i1) ! emitter
        tqW(i1)=qW(i1)+pext(proc_nlegs,i1)-(1-xij_off)*pext(i,i1) ! spectator
	enddo

! the remaining initial-state particle
	do j1 = 1, 2
	if (j1 .ne. i) then
	 do i1 = 1,4
	 tq(j1,i1)=pext(j1,i1)
	 enddo
	endif
	enddo

! the W decay leptons
! Note: pext(7,4) are the qq -> 4l + j off-shell momenta
      do i1 = 1,4
       p(i1,1) = pext(3,i1)
       p(i1,2) = pext(4,i1)
      enddo
      MW_off = dot_p2(tqW(1))
      MW_off = sqrt(MW_off) 
      call OSMap_decayV_LO(p(1,1),MW_off,tqW(1),q(1,1))
      do i1 = 1, 4
       tq(3,i1) = q(i1,1)
       tq(4,i1) = q(i1,2)
      enddo

! the Z decay leptons
	do j1 = 5, 6
	 do i1 = 1,4
	 tq(j1,i1)=pext(j1,i1)
	 enddo
	enddo	
***]get the off-shell 2to4 momenta

      call DoTildeAmpHis_OffS(tq,1,Gfr,Jacb)
        return 
        endif
#endif
***] initial (i-emitter) - final (j-spectator)
	
***[ final(emitter) - final(spectator)
#ifdef NINH_CHECK_FIN_FIN
	if (i .ge. 3 .and. j .ge. 3) then
!!![needed for WW case
! on-shell momenta for Gfr (same as final-initial and initial-final)
! oext(nlegs_os_max,4): on-shell momenta in the partonic CMS
! the photon: oext(nlegs,:) = on-shell = off-shell, nlegs = 5
c        print*,"mi, mj =",i,j,sqrt(mi2),sqrt(mj2)

        do i1 = 1,4
	Pij(i1)=oext(i,i1)+oext(j,i1)+oext(nlegs,i1)  
	enddo
       Pij2 = dot_p2(Pij(1))
	BPij2 = Pij2 - mi2 -mj2
	lamij = Pij2**2+mi2**2+mj2**2-2*(Pij2*mi2+Pij2*mj2+mi2*mj2)
	rlamij = sqrt(lamij)
	
        PiK = oext(i,1)*oext(nlegs,1)-oext(i,2)*oext(nlegs,2)
     &  -oext(i,3)*oext(nlegs,3)-oext(i,4)*oext(nlegs,4)

        PjK = oext(j,1)*oext(nlegs,1)-oext(j,2)*oext(nlegs,2)
     &  -oext(j,3)*oext(nlegs,3)-oext(j,4)*oext(nlegs,4)
     
        PiPj = oext(i,1)*oext(j,1)-oext(i,2)*oext(j,2)
     &  -oext(i,3)*oext(j,3)-oext(i,4)*oext(j,4)	
	
c 
        yij = PiK/(PiK+PjK+PiPj)
	zij = PiPj/(PiPj+PjK) 
	
	if ( mi2+mj2 .eq. 0 ) then
	Gfr = (2D0/(1-zij*(1-yij))-1-zij)/PiK/(1-yij)
	else
	
        Gfr = sqrt((2*mj2+BPij2-BPij2*yij)**2-4*Pij2*mj2)/rlamij
	Gfr = (2D0/(1-zij*(1-yij))-1-zij-mi2/PiK)/(PiK*Gfr)
        endif
        
***[the off-shell correspondence
	do i1 = 1,4
         pW(3,i1) = pext(3,i1) + pext(4,i1) ! off-shell W1
         pW(4,i1) = pext(5,i1) + pext(6,i1) ! off-shell W2
        enddo
!!![off-shell masses
	temp1 = sqrt(pW(i,2)**2 + pW(i,3)**2 + pW(i,4)**2)
        mi2_off = (pW(i,1)-temp1)*(pW(i,1)+temp1)

        temp2 = sqrt(pW(j,2)**2 + pW(j,3)**2 + pW(j,4)**2)
        mj2_off = (pW(j,1)-temp2)*(pW(j,1)+temp2)

	if (abs(mi2_off) .lt. 1D-10*temp1**2) mi2_off = 0
	if (abs(mj2_off) .lt. 1D-10*temp2**2) mj2_off = 0
!!!]off-shell masses

        do i1 = 1,4
	Pij_off(i1)=pW(3,i1)+pW(4,i1)+pext(proc_nlegs,i1)  
	enddo
	Pij2_off = dot_p2(Pij_off(1))
	BPij2_off = Pij2_off - mi2_off - mj2_off
	lamij_off = Pij2_off**2+mi2_off**2+mj2_off**2
     & -2*(Pij2_off*mi2_off+Pij2_off*mj2_off+mi2_off*mj2_off)
	rlamij_off = sqrt(lamij_off)
	
       PiK_off = pW(i,1)*pext(proc_nlegs,1)-pW(i,2)*pext(proc_nlegs,2)
     &  -pW(i,3)*pext(proc_nlegs,3)-pW(i,4)*pext(proc_nlegs,4)
     
       PjK_off = pW(j,1)*pext(proc_nlegs,1)-pW(j,2)*pext(proc_nlegs,2)
     &  -pW(j,3)*pext(proc_nlegs,3)-pW(j,4)*pext(proc_nlegs,4)
     
       PiPj_off = pW(i,1)*pW(j,1)-pW(i,2)*pW(j,2)
     &  -pW(i,3)*pW(j,3)-pW(i,4)*pW(j,4)	

       yij_off = PiK_off/(PiK_off+PjK_off+PiPj_off)
	if ( mi2_off+mj2_off .eq. 0 ) then
	do i1 = 1,4
        tpW(j,i1)=pW(j,i1)/(1-yij_off)
        tpW(i,i1)=pW(i,i1)+pext(proc_nlegs,i1)-pW(j,i1)*yij_off/(1-yij_off)
	enddo
	
        else 
       temp = 2*PiK_off+mi2_off
       temp = temp**2+Pij2_off**2+mj2_off**2
     & -2*(temp*Pij2_off+temp*mj2_off+Pij2_off*mj2_off)
       temp = sqrt(temp)
       temp1 = Pij_off(1)*pW(j,1)-Pij_off(2)*pW(j,2)
     &	-Pij_off(3)*pW(j,3)-Pij_off(4)*pW(j,4) 
	
	do i1 = 1,4
	tpW(j,i1)=(Pij2_off+mj2_off-mi2_off)/(2*Pij2_off)*Pij_off(i1)
     &  + rlamij_off/temp*(pW(j,i1)-temp1/Pij2_off*Pij_off(i1))
        tpW(i,i1)=Pij_off(i1)-tpW(j,i1) 	
	enddo
	
	endif

! initial quark momenta 	
	do j1 = 1, 2
	 do i1 = 1,4
	 tq(j1,i1)=pext(j1,i1)
	 enddo
	enddo
	
c	do i1 = 1, 4
c	print*,"check EP =",i1,tq(1,i1)+tq(2,i1)-tpW(i,i1)-tpW(j,i1)	
c        enddo
c        stop

!!![W decays
! the W1 decay leptons
      do i1 = 1,4
       p(i1,1) = pext(3,i1)
       p(i1,2) = pext(4,i1)
      enddo
      tqW(:)=tpW(3,:)
      MW_off = dot_p2(tqW(1))
      MW_off = sqrt(MW_off)
      call OSMap_decayV_LO(p(1,1),MW_off,tqW(1),q(1,1))
      do i1 = 1,4
        tq(3,i1) = q(i1,1)
        tq(4,i1) = q(i1,2)
      enddo

! the W2 decay leptons
      do i1 = 1,4
       p(i1,1) = pext(5,i1)
       p(i1,2) = pext(6,i1)
      enddo
      tqW(:)=tpW(4,:)
      MW_off = dot_p2(tqW(1))
      MW_off = sqrt(MW_off)
      call OSMap_decayV_LO(p(1,1),MW_off,tqW(1),q(1,1))
      do i1 = 1,4
        tq(5,i1) = q(i1,1)
        tq(6,i1) = q(i1,2)
      enddo
!!!]W decays
***]the off-shell correspondence

c	do i1 = 1, 4
c	print*,"check EP2 =",i1,tq(1,i1)+tq(2,i1)-tq(3,i1)-tq(4,i1)-tq(5,i1)-tq(6,i1)	
c        enddo
c        stop

      call DoTildeAmpHis_OffS(tq,1,Gfr,Jacb)
!!!]needed for WW case

        return
	endif
#endif
***] final(emitter) - final(spectator)

	end
