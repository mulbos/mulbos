***[ Use Ninh code
	double precision function FuncDPsing_fin_ind(fac)
	use global_vars
	implicit none

#include "model.h"

	double precision Gfi
	external Gfi

	integer i, j, nlegs
        parameter (nlegs = 5) ! for the on-shell 2to3 process

!        double precision MVint2(2) ! the off-shell invariant mass of the gauge bosons
	double precision charge(nlegs), r
	
	double precision Jacb, fac

      charge = 0D0
      if (lNLO_EW) then	
       charge(3) = charge_qed(3) ! final jet at the 3rd position
      else
       charge(3) = charge_qcd(3) ! final jet at the 3rd position
      endif

C We need the off-shell gauge-boson invariant masses, 
C to pass to the reduced amplitudes. 
C Here we take the simplest choice that those invariant 
C masses are independent of dipoles (= emitter + spectator).
C Note: pext(7,4) are the true q g -> q + 4l momenta
!	MVint2(1) = pext(4,1)*pext(5,1) - pext(4,2)*pext(5,2) 
!     &  - pext(4,3)*pext(5,3) - pext(4,4)*pext(5,4)
!        MVint2(1) = MVint2(1)*2
!	MVint2(2) = pext(6,1)*pext(7,1) - pext(6,2)*pext(7,2) 
!     &  - pext(6,3)*pext(7,3) - pext(6,4)*pext(7,4)
!        MVint2(2) = MVint2(2)*2

! WARNING: this is process dependent
	r = charge(3)**2
	r = r*fac ! include a minus sign for being subtraction term
	
c	print*,"new r =", r
c	stop

! change here for Photon (1) or Gluon (4/3 = CF) 
      if (lNLO_QCD) then
       Jacb = r*4*pi*4/3D0
      else if (lNLO_EW) then
       Jacb = r*4*pi
      endif
	    
	    FuncDPsing_fin_ind = Jacb*Gfi(Jacb)

c        print*,'from FuncDPsing_ind.F =', FuncDPsing_fin_ind, Jacb, r, fac, Alfas

	end

************************************************************************
* Ref0. arXiv:hep-ph/0509084v1 (Diener, Dittmaier, Hollik)
* Ref1. arXiv: 0802.1405v1 (Stefan Dittmaier, Alois Kabelschacht, Tobias Kasprzik)

	double precision function Gfi(Jacb)
! Gfi includes Alfas and pdfs
        use global_vars
        use global_mom_xpdf
	implicit none
	integer i1, j1, i

        double precision Jacb
	double precision Pij2
	double precision temp, temp1, temp2
	double precision PiK, PjK, PiPj		
	double precision xij	
        double precision Pij(4), TPij(4)
        double precision tq(proc_nlegs-1,4) ! off-shell tilde momenta
        double precision dot_p2
       double precision tpp(4), Hfi

	Gfi = 0
	Hfi = 0
	
! oext(nlegs_os_max,4): jet at the 3rd position, V(4), V(5), leptons(6,7,8,9)

*** initial(emitter) - initial(spectator) ***
* For gamma -> ff splitting, 
* gamma is the 2nd particle.

C check particle masses ordering
c#define CHECK_MASS_ORDER
#ifdef CHECK_MASS_ORDER
	do i1 = 1,5
        temp = oext(i1,1)*oext(i1,1)-oext(i1,2)*oext(i1,2)
     &  -oext(i1,3)*oext(i1,3)-oext(i1,4)*oext(i1,4)

        print*,"mass2 =", i1, temp
        enddo
        stop
#endif

***[for NLO QCD and EW
! mapping method to re-produce Denner WZ 2020 results. 
! idea is to use exact momenta to calculate tilde kinematics.
! Convergence is excellent.

C Note: pext(7,4) are the qq -> j + 4l off-shell momenta

        PiK = pext(1,1)*pext(2,1)-pext(1,2)*pext(2,2)
     &  -pext(1,3)*pext(2,3)-pext(1,4)*pext(2,4)

        PjK = pext(3,1)*pext(2,1)-pext(3,2)*pext(2,2)
     &  -pext(3,3)*pext(2,3)-pext(3,4)*pext(2,4)

        PiPj = pext(1,1)*pext(3,1)-pext(1,2)*pext(3,2)
     &  -pext(1,3)*pext(3,3)-pext(1,4)*pext(3,4)

	xij = (PiK-PjK-PiPj)/PiK

        Gfi = (1+2*xij**2-2*xij)/(PjK*xij)

        do i1 = 1,4
	Pij(i1)=pext(1,i1)+pext(2,i1)-pext(3,i1)  
	enddo
	Pij2 = dot_p2(Pij(1))

***[get the q1q2 > l3 l4 l5 l6 momenta (tilde)
	do i1 = 1,4
        tq(2,i1)=pext(2,i1)*xij
        tq(1,i1)=pext(1,i1)		
	
	TPij(i1)=pext(1,i1)+pext(2,i1)*xij
	enddo

	temp = Pij(1)*TPij(1)-Pij(2)*TPij(2)
     &	-Pij(3)*TPij(3)-Pij(4)*TPij(4)
        
	do j1 = 3, proc_nlegs-1
	temp1 = Pij(1)*pext(j1+1,1)-Pij(2)*pext(j1+1,2)
     &	-Pij(3)*pext(j1+1,3)-Pij(4)*pext(j1+1,4) 	 	
     
	temp2 = TPij(1)*pext(j1+1,1)-TPij(2)*pext(j1+1,2)
     &	-TPij(3)*pext(j1+1,3)-TPij(4)*pext(j1+1,4)
     
        do i1 = 1,4
	tq(j1,i1)=pext(j1+1,i1)
     &	-(Pij(i1)+TPij(i1))*(temp1+temp2)/(Pij2+temp)
     &  +2*TPij(i1)*temp1/Pij2
	enddo  	 	     
     
        enddo
***]get the q1q2 > l3 l4 l5 l6 momenta (tilde)

      call DoTildeAmpHis_OffS(tq(1,1),1,Gfi,Jacb)
***]for NLO QCD and EW

***[ for f -> f' gam^*
      if (lWW.and.lNLO_EW) then
C to be done for W+ W-
c         print*,"from FuncDPsing_fin_ind.F: PHOTON PHOTON contribution"
c         stop

! f is the 1st particle
! Note: pext(7,4) are the qq -> j + 4l off-shell momenta
! tq(proc_nlegs-1,4): off-shell tilde momenta
***[ mapping
c        print*,'xij =', xij
	do i1 = 1,4
        tq(2,i1) = pext(2,i1) ! spectator photon
        tq(1,i1) = pext(1,i1)*xij ! photon from f
	TPij(i1) = tq(1,i1) + tq(2,i1)
* After Ref1. Eq. (5.11), exact momenta in the partonic CMS:
        tpp(i1) = pext(3,i1) - PiPj/PiK*pext(2,i1)
	enddo

	temp = Pij(1)*TPij(1)-Pij(2)*TPij(2)
     &	-Pij(3)*TPij(3)-Pij(4)*TPij(4)
        
	do j1 = 3, proc_nlegs-1
	temp1 = Pij(1)*pext(j1+1,1)-Pij(2)*pext(j1+1,2)
     &	-Pij(3)*pext(j1+1,3)-Pij(4)*pext(j1+1,4) 	 	
     
	temp2 = TPij(1)*pext(j1+1,1)-TPij(2)*pext(j1+1,2)
     &	-TPij(3)*pext(j1+1,3)-TPij(4)*pext(j1+1,4)
     
        do i1 = 1,4
	tq(j1,i1) = pext(j1+1,i1)
     &	-(Pij(i1)+TPij(i1))*(temp1+temp2)/(Pij2+temp)
     &  +2*TPij(i1)*temp1/Pij2
	enddo  	 	     
     
        enddo
***] mapping

* Ref1. Eq. (5.9)
* sum over the f' helicities (kappa_f) leads to a factor 2
	temp1 = 1D0/2/PiPj*2
	temp = tpp(1)**2 - tpp(2)**2 - tpp(3)**2 - tpp(4)**2
	temp2 = -4*(1 - xij)/xij**2/temp

c        print*,'temp1, temp2 =', temp1, temp2

      temp1 = temp1*3 ! adding color factor for quarks

      call DoTildeAmpHis_OffS_aa(tq(1,1),1,Hfi,Jacb,temp1,temp2,tpp)
      Gfi = Gfi + Hfi
c      stop
      endif
***] for f -> f' gam^*

	end
***] Use Ninh code

