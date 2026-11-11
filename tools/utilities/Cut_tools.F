* Author: Le Duc Ninh
* Note: some routines in this file are taken from VBFNLO
* (search "VBFNLO" for details)
* last modification 2016

***[interface
* For Tree, Virt, Dipole(Slicing)-Finite cases

	subroutine DoCut_parCM(cfail, fext)
	use global_vars
	use global_mom_xpdf
C input momenta in partonic center of mass system
	implicit none
	integer i, j, cfail

c#include "xsection.h"	

        double precision rgam, rbetgam
	double precision x1, x2
	double precision rext(proc_nlegs,4), fext(proc_nlegs,4)

c        double precision qext(4,proc_nlegs)
c	double precision pt(3:6), py(3:6), phi(3:6)
c        double precision pTW, pTZ, MinvW, MinvZ, cos_theta_W, cos_theta_Z
c        double precision pW(4), pZ(4), pTW2, pTZ2, ppW, ppZ, ppW2, ppZ2

	do j = 1,4 
	 do i = 1,proc_nlegs
	  rext(i,j) = fext(i,j)
	 enddo
	enddo

***{boost parameters
* mq = mP = 0 here
* boost for the subprocess: P1/PARTON1 + P2/PARTON2 -> XXX
* Rules: OZ(x1), -OZ(x2)
	x2 = tau_pdf**2
	x2 = x2 + (1-x2)*x_pdf

        x1 = tau_pdf**2/x2

        rgam = 0.5D0*(x1 + x2)/sqrt(x1*x2)        

        rbetgam=0.5D0*(x2 - x1)/sqrt(x1*x2)
	do i = 1, proc_nlegs
	call boost0(rgam, rbetgam, 
     &  rext(i,1), rext(i,2), rext(i,3), rext(i,4))
	enddo
***}boost parameters
	call DoCut(proc_nlegs, cfail, rext)
	end

* For NLO collinear (convolution) part.
* One has to boost since each initial particle can 
* radiate photon (gluon).

	subroutine DoCut_parCM_coll(cfail, fext, ie, ifail)
	use global_vars
	use global_mom_xpdf
	implicit none
	integer i, j, ie, ifail, cfail

c#include "xsection.h"	

        double precision qext(proc_nlegs,4)
	double precision rext(proc_nlegs,4), fext(proc_nlegs,4)
        double precision rgam,rbetgam,x1,x2

***{boost parameters 1
* ME = 0 here
* to be in the CMS of tau_pdf*ecm

        if (ie .gt. 2) then
	 print*,'DoCut_CMS() error: ie > 2'
	 print*,'ie =', ie
	 stop
	endif 

* Either
	if ( ifail .eq. 1 ) then
        x1 = y_coll
	else
        x1 = 1D0
	endif
* Or
* This sets x1 = 1 (Boost off) automatically 
* for the [x1=1]-term in the plus distribution. 
c        x1 = (sqrtS/(tau_pdf*ecm))**2

        rgam = 0.5D0*(x1+1)/sqrt(x1)        

        if (ie .eq. 1) then 
        rbetgam = 0.5D0*(1-x1)/sqrt(x1)
        else     
        rbetgam = -0.5D0*(1-x1)/sqrt(x1)
	endif
***}boost parameters 1

	do j = 1,4 
	 do i = 1,proc_nlegs
	  qext(i,j) = fext(i,j)
	 enddo
	enddo

        do i=1,proc_nlegs
        call boost0(rgam, rbetgam, 
     &  qext(i,1),qext(i,2),qext(i,3),qext(i,4))
	enddo	

	do j = 1,4 
	 do i = 1,proc_nlegs
	  rext(i,j) = qext(i,j)
	 enddo
	enddo

***{boost parameters 2a
* mq = mP = 0 here
* boost for the subprocess: P1/PARTON1 + P2/PARTON2 -> XXX
* Rules: OZ(x1), -OZ(x2)
	x2 = tau_pdf**2
	x2 = x2 + (1-x2)*x_pdf

        x1 = tau_pdf**2/x2

        rgam = 0.5D0*(x1 + x2)/sqrt(x1*x2)        

        rbetgam=0.5D0*(x2 - x1)/sqrt(x1*x2)
	do i = 1, proc_nlegs
	call boost0(rgam, rbetgam, 
     &  rext(i,1), rext(i,2), rext(i,3), rext(i,4))
	enddo
***}boost parameters 2a
	call DoCut(proc_nlegs, cfail, rext)
        return

        end
***]interface


***[technical cuts
	subroutine technical_cut_fin_rad(lcut_fail)
	use global_vars
	use global_mom_xpdf
	implicit none
C arguments
      logical lcut_fail
      double precision s12, s_IR, eps
      integer l1, l2, l3, l4, phot
      integer i
      
      lcut_fail = .false. ! initial value
        
      eps = 1D-8  ! cut threshold
      phot = 5 ! photon at the 5th position
      l1=6
      l2=7
      l3=8
      l4=9
        
      IF (lDECAYV2_NLO_EW) THEN
       s12 = mass_VV(2)**2
       if (lWZ.or.lZZ) then
       do i = l3, l4
        s_IR = oext(phot,1)*oext(i,1)-oext(phot,2)*oext(i,2)-oext(phot,3)*oext(i,3)-oext(phot,4)*oext(i,4)
        s_IR = abs(s_IR*2)
        if ( s_IR .le. s12*eps ) then
         lcut_fail = .true.
         return
        endif
       enddo
       else ! W case
       i=l3
        s_IR = oext(phot,1)*oext(i,1)-oext(phot,2)*oext(i,2)-oext(phot,3)*oext(i,3)-oext(phot,4)*oext(i,4)
        s_IR = abs(s_IR*2)
        if ( s_IR .le. s12*eps ) then
         lcut_fail = .true.
         return
        endif
       endif
      
      ELSE IF (lDECAYV1_NLO_EW) THEN
       s12 = mass_VV(1)**2
       if (lZZ) then
       do i = l1, l2
        s_IR = oext(phot,1)*oext(i,1)-oext(phot,2)*oext(i,2)-oext(phot,3)*oext(i,3)-oext(phot,4)*oext(i,4)
        s_IR = abs(s_IR*2)
        if ( s_IR .le. s12*eps ) then
         lcut_fail = .true.
         return
        endif
       enddo
       else ! W case
       i=l1
        s_IR = oext(phot,1)*oext(i,1)-oext(phot,2)*oext(i,2)-oext(phot,3)*oext(i,3)-oext(phot,4)*oext(i,4)
        s_IR = abs(s_IR*2)
        if ( s_IR .le. s12*eps ) then
         lcut_fail = .true.
         return
        endif
       endif       
       
      ELSE
       s12 = oext(1,1)*oext(2,1)-oext(1,2)*oext(2,2)-oext(1,3)*oext(2,3)-oext(1,4)*oext(2,4)
       s12 = abs(s12*2)
       do i = 1, 2 ! initial quarks
        s_IR = oext(phot,1)*oext(i,1)-oext(phot,2)*oext(i,2)-oext(phot,3)*oext(i,3)-oext(phot,4)*oext(i,4)
        s_IR = abs(s_IR*2)
        if ( s_IR .le. s12*eps ) then
         lcut_fail = .true.
         return
        endif
       enddo
      ENDIF 

        end subroutine
        
      	subroutine technical_cut_fin_ind(lcut_fail)
	use global_vars
	use global_mom_xpdf
	implicit none
C arguments
      logical lcut_fail
      double precision s12, s_IR, eps
      integer jet
      integer i
      
      lcut_fail = .false. ! initial value
        
      eps = 1D-8  ! cut threshold
      jet = 3 ! quark jet at the 3rd position
      
      s12 = oext(1,1)*oext(2,1)-oext(1,2)*oext(2,2)-oext(1,3)*oext(2,3)-oext(1,4)*oext(2,4)
      s12 = abs(s12*2)
      
       do i = 1, 2 ! initial quarks
        s_IR = oext(jet,1)*oext(i,1)-oext(jet,2)*oext(i,2)-oext(jet,3)*oext(i,3)-oext(jet,4)*oext(i,4)
        s_IR = abs(s_IR*2)
       if ( s_IR .le. s12*eps ) then
        lcut_fail = .true.
        return
       endif
       enddo

        end subroutine
***]technical cuts

***[Gram determinant cuts
	subroutine Check_GramDet_OS_1234(lcut_fail)
	use global_vars
	use global_mom_xpdf
	implicit none
C arguments
      logical lcut_fail
      double precision detG,p1(4),p2(4),p3(4)
      double precision G11, G22, G33, G12, G13, G23
      double precision dot, dot_p2      

      lcut_fail = .false. ! initial value

! oext(nlegs_os_max,4)
      p1(:) = oext(1,:)
      p2(:) = oext(2,:)
      p3(:) = oext(3,:)
      G11 = dot_p2(p1)
      G22 = dot_p2(p2)
      G33 = dot_p2(p3)
      G12 = dot(p1,p2)
      G13 = dot(p1,p3)
      G23 = dot(p2,p3)
      detG=G11*G22*G33+2*G12*G23*G13
     &  -G13*G22*G13-G11*G23*G23
     &  -G33*G12*G12
      detG = abs(detG/G12**3)
       if ( detG .le. eps_GramDet4 ) then
        lcut_fail = .true.
        return
       endif
        end subroutine
***]Gram determinant cuts
***********************************************************************
	subroutine DoCut(Next, cfail, fext)
	use global_vars
	implicit none
	integer Next, cfail
	double precision fext(Next,4)
	if (lDiBoson) then
	 call DoCut_DiBoson(Next, cfail, fext)
	else
	 call DoCut_TriBoson(Next, cfail, fext)
	endif
        end subroutine DoCut
***********************************************************************
	subroutine DoCut_TriBoson(Next, cfail, fext)
	use global_vars
	use global_mom_xpdf
	use global_cuts
	implicit none
#include "model.h"
        integer Next, cfail
        double precision fext(Next,4)
        double precision qext(4,Next)
        double precision q_clep(8,Nlep_charged) ! charged lep mom and particle id number
	double precision pt(1:Next), py(1:Next), phi(1:Next)
	double precision qV(4,3), qVV(4,3), temp
        double precision inv_mass, inv_mass_p, dot_p2
	integer i, j

!        print*,"DoCut_TriBoson hellos, cfail =",cfail

! the order: q_1 q_2 > l_3 + l_4 + ... + jet:
        do j = 1,4
         do i = 1, Next
          qext(j,i) = fext(i,j)
         enddo
        enddo

***[all final-sate particles
      do i = 3, Next
       call lego( qext(1,i), pt(i), py(i), phi(i) )
      enddo
***]all final-sate particles

      call Get_charged_leptons(Next,qext(1,1),pt(1),py(1),phi(1),q_clep(1,1))

C for individual leptons
      do i = 3, 8
        if ( pt(i) .le. pTmin_L(i-2) .or. abs(py(i)) .ge. pymax_L(i-2) ) then
         cfail = 0
         return       
        endif
      enddo

***[momentum of resonances
      qV(:,1) = qext(:,3) + qext(:,4)
      qV(:,2) = qext(:,5) + qext(:,6)
      qV(:,3) = qext(:,7) + qext(:,8)
      qVV(:,3) = qV(:,1) + qV(:,2) ! 4L_1234
      qVV(:,2) = qV(:,1) + qV(:,3) ! 4L_1256
      qVV(:,1) = qV(:,2) + qV(:,3) ! 4L_3456
***]momentum of resonances

***[4-lepton invariant mass cut
      do i = 1, 3
      temp = inv_mass_p(qVV(1,i))
      if ( temp .le. minv_4lep_min(i) .or. temp .ge. minv_4lep_max(i) ) then
         cfail = 0
         return       
      endif
      enddo
***]4-lepton invariant mass cut

***[DR separations
      do i=1,Nlep_charged-1
       do j=i+1,Nlep_charged
        call rsep_simple(q_clep(6,i),q_clep(7,i),q_clep(6,j),q_clep(7,j),temp)
        if ( temp .le. DRmin_LL_all ) then
         cfail = 0
         return
        endif
       enddo
      enddo
***]DR separations


        end subroutine DoCut_TriBoson
***********************************************************************
	subroutine DoCut_DiBoson(Next, cfail, fext)
	use global_vars
	use global_mom_xpdf
	use global_cuts
	implicit none
#include "model.h"

C Next = 6 for LO and 7 for real emission
	integer Next, cfail, itemp

	double precision fext(Next,4)
	double precision qext(4,Next), clep_bared(7,Nlep_charged), 
     & cLep_dressed(7,Nlep_charged), dot
        double precision qWZ(4), m2_WZ, m_WZ, qT_WZ

	double precision m2_Z, m_Z, mT_Z
	double precision m2_W, m_W, mT_W
	double precision DRLL_34, DRLL_35, DRLL_36, DRLL_45, DRLL_46, DRLL_56
	double precision pt(3:6), py(3:6), phi(3:6)
	double precision pt_phot, py_phot, phi_phot
	double precision phot(7), pt_jet, py_jet, phi_jet
        integer cut_i
        double precision p_W(4), p_Z(4), pT_W, pT_Z
        double precision p_miss(4), pT_miss, pt_lead, temp, p_temp(4)
        double precision pt_lead2

	integer i, j

        do j = 1,4
         do i = 1, Next
          qext(j,i) = fext(i,j)
         enddo
        enddo

***[all leptons
      do i = 3, 6
       call lego( qext(1,i), pt(i), py(i), phi(i) )
      enddo
***]all leptons

***[photon-chargedlepton recombination
c Final-state photons with rapidity |y| > 5 are
c considered as lost in the beam pipe and excluded from recombination. 
      IF (lNLO_FIN_RAD) THEN
      IF (lNLO_EW) THEN

      call lego( qext(1,Next), pt_phot, py_phot, phi_phot )
      If ( abs(py_phot) .le. py_phot_recomb ) then

      do j = 1,4
       if (lWZ) then
        clep_bared(j,1) = qext(j,3)
        clep_bared(j,2) = qext(j,5)
        clep_bared(j,3) = qext(j,6)
       else if (lWW) then
        clep_bared(j,1) = qext(j,3)
        clep_bared(j,2) = qext(j,5)
       else if (lZZ) then
        clep_bared(j,1) = qext(j,3)
        clep_bared(j,2) = qext(j,4)
        clep_bared(j,3) = qext(j,5)
        clep_bared(j,4) = qext(j,6)
       endif
       phot(j) = qext(j,Next)
      enddo

      phot(5) = pt_phot
      phot(6) = py_phot
      phot(7) = phi_phot

      clep_bared(5,1) = pt(3) ! electron
      clep_bared(6,1) = py(3)
      clep_bared(7,1) = phi(3)

      if (lWZ) then
      do i = 2, 3 ! muons
      j = i+3
      clep_bared(5,i) = pt(j)
      clep_bared(6,i) = py(j)
      clep_bared(7,i) = phi(j)
      enddo
      else if (lWW) then
      clep_bared(5,2) = pt(5)
      clep_bared(6,2) = py(5)
      clep_bared(7,2) = phi(5)
      else if (lZZ) then
      do i = 2, 4
      j = i+2
      clep_bared(5,i) = pt(j)
      clep_bared(6,i) = py(j)
      clep_bared(7,i) = phi(j)
      enddo
      endif

      call PhotLepComb(Nlep_charged,DRLA,phot(1),clep_bared(1,1),clep_dressed(1,1))

      do j = 1, 4
       if (lWZ) then
        qext(j,3) = clep_dressed(j,1)
        qext(j,5) = clep_dressed(j,2)
        qext(j,6) = clep_dressed(j,3)
       else if (lWW) then
        qext(j,3) = clep_dressed(j,1)
        qext(j,5) = clep_dressed(j,2)
       else if (lZZ) then
        qext(j,3) = clep_dressed(j,1)
        qext(j,4) = clep_dressed(j,2)
        qext(j,5) = clep_dressed(j,3)
        qext(j,6) = clep_dressed(j,4)
       endif
      enddo
! for W1 electron
       pt(3) = clep_dressed(5,1)
       py(3) = clep_dressed(6,1)
       phi(3) = clep_dressed(7,1)
       if (lWZ) then
        do i = 5, 6
         j = i-3
         pt(i) = clep_dressed(5,j)
         py(i) = clep_dressed(6,j)
         phi(i) = clep_dressed(7,j)
        enddo
       else if (lWW) then
! for W2 muon
        pt(5) = clep_dressed(5,2)
        py(5) = clep_dressed(6,2)
        phi(5) = clep_dressed(7,2)
       else if (lZZ) then
        do i = 4, 6
         j = i-2
         pt(i) = clep_dressed(5,j)
         py(i) = clep_dressed(6,j)
         phi(i) = clep_dressed(7,j)
        enddo
       endif

      Endif
      ENDIF
      ENDIF
***]photon-chargedlepton recombination

!!![leading-pt lepton
      if (lWW) then
       pt_lead = max(pt(3),pt(5))
       pt_lead2 = min(pt(3),pt(5))
      else if (lWZ) then
       pt_lead = max(pt(3),pt(5),pt(6))
       pt_lead2 = min(pt(3),pt(5),pt(6))
       pt_lead2 = pt(3)+pt(5)+pt(6)-pt_lead-pt_lead2
      else if (lZZ) then
       ! 3456
       pt_lead = max(pt(3),pt(4),pt(5))
       temp = min(pt(3),pt(4),pt(5))
       pt_lead2 = pt(3)+pt(4)+pt(5)-pt_lead-temp
       if (pt_lead .lt. pt(6)) then
        pt_lead2 = pt_lead
        pt_lead = pt(6)
       else if (pt_lead2 .lt. pt(6)) then
        pt_lead2 = pt(6)
       endif
      endif
!!!]leading-pt lepton

C need qext(4,leptons), pt(leptons), py(leptons), phi(leptons)

***[jet veto
! apply for NLO QCD RAD OR IND
! and for NLO EW IND (not for photon radiation)
c Final-state jets with rapidity |y| > 4.5 are
c considered as lost in the beam pipe and excluded from jet veto. 
      IF ((lNLO_FIN.and.lNLO_QCD).or.lNLO_FIN_IND) THEN
       call lego( qext(1,7), pt_jet, py_jet, phi_jet )
       if (py_jet.le.py_jet_veto_max.and.pt_jet.ge.pt_jet_veto_min) then
        cfail = 0
        return
       endif
      ENDIF
***]jet veto

***[4-lepton invariant mass cut
      do i = 1, 4
       qWZ(i) = qext(i,3) + qext(i,4) + qext(i,5) + qext(i,6)
      enddo
      m2_WZ = dot(qWZ(1),qWZ(1))
      m_WZ = sqrt(m2_WZ)

      if ( m_WZ .le. minv_4l_min .or. m_WZ .ge. minv_4l_max ) then
         cfail = 0
         return       
      endif
***]4-lepton invariant mass cut

***[invariant mass of 2 leptons
      temp = 2*dot(qext(1,3),qext(1,5))
      if ( temp .le. minv_LL13_min**2 .or. temp .ge. minv_LL13_max**2) then
        cfail = 0
        return
      endif

      m2_W = 2*dot(qext(1,3),qext(1,4)) ! V1
      if ( m2_W .le. minv_LL12_min**2 .or. m2_W .ge. minv_LL12_max**2) then
        cfail = 0
        return
      endif

      m2_Z = 2*dot(qext(1,5),qext(1,6)) ! V2
      if ( m2_Z .le. minv_LL34_min**2 .or. m2_Z .ge. minv_LL34_max**2 ) then
         cfail = 0
         return       
      endif
      
C V1, V2 invariant mass shift
      m_W = sqrt(m2_W)
      if ( m_W .le. (mass_VV(1)-m_V1_shift) .or. m_W .ge. (mass_VV(1)+m_V1_shift) ) then
         cfail = 0
         return
      endif
      
      m_Z = sqrt(m2_Z)
      if ( m_Z .le. (mass_VV(2)-m_V2_shift) .or. m_Z .ge. (mass_VV(2)+m_V2_shift) ) then
         cfail = 0
         return       
      endif
      
c Ref: https://arxiv.org/abs/1606.04017
      mT_W = cos(phi(3) - phi(4)) ! V1
      mT_W = sqrt( 2*pt(3)*pt(4)*(1 - mT_W) )
      if ( mT_W .le. mT_V1_min ) then
         cfail = 0
         return       
      endif
      
      mT_Z = cos(phi(5) - phi(6)) ! V2
      mT_Z = sqrt( 2*pt(5)*pt(6)*(1 - mT_Z) )
      if ( mT_Z .le. mT_V2_min ) then
         cfail = 0
         return
      endif
***]invariant mass of 2 leptons

C for individual leptons
      do i = 3,6
        if ( pt(i) .le. pTmin_L(i-2) .or. abs(py(i)) .ge. pymax_L(i-2) ) then
         cfail = 0
         return       
        endif
      enddo
      
! leading-pt lepton cut
      if (pt_lead .le. pTmin_Llead) then
        cfail = 0
        return
      endif
      if (pt_lead2 .le. pTmin_Llead2) then
        cfail = 0
        return
      endif

C DR separations
cc      DRLL_Z = sqrt((py(5) - py(6))**2 + (phi(5) - phi(6))**2)
      call rsep_simple(py(3),phi(3),py(4),phi(4),DRLL_34)
      call rsep_simple(py(3),phi(3),py(5),phi(5),DRLL_35)
      call rsep_simple(py(3),phi(3),py(6),phi(6),DRLL_36)
      call rsep_simple(py(4),phi(4),py(5),phi(5),DRLL_45)
      call rsep_simple(py(4),phi(4),py(6),phi(6),DRLL_46)
      call rsep_simple(py(5),phi(5),py(6),phi(6),DRLL_56)
      if ( DRLL_34 .le. DRmin_LL_12 ) then
         cfail = 0
         return       
      endif
      if ( DRLL_35 .le. DRmin_LL_13 ) then
         cfail = 0
         return       
      endif
      if ( DRLL_36 .le. DRmin_LL_14 ) then
         cfail = 0
         return       
      endif
      if ( DRLL_45 .le. DRmin_LL_23 ) then
         cfail = 0
         return       
      endif
      if ( DRLL_46 .le. DRmin_LL_24 ) then
         cfail = 0
         return       
      endif
      if ( DRLL_56 .le. DRmin_LL_34 ) then
         cfail = 0
         return       
      endif
      
***[pTWZ cut
      qT_WZ = sqrt(qWZ(2)**2 + qWZ(3)**2)
      if ( qT_WZ .lt. pT_4L_min .or. qT_WZ .ge. pT_4L_max) then
         cfail = 0
         return
      endif      
***]pTWZ cut

***[pTV cut
      do i = 1, 4
      p_Z(i) = qext(i,5) + qext(i,6)
      p_W(i) = qext(i,3) + qext(i,4)
      enddo
      pT_Z = sqrt(p_Z(2)**2 + p_Z(3)**2)
      pT_W = sqrt(p_W(2)**2 + p_W(3)**2)
      if (pT_Z .lt. pT_LL34_min .or. pT_Z .ge. pT_LL34_max) then
         cfail = 0
         return
      endif      
      if (pT_W .lt. pT_LL12_min .or. pT_W .ge. pT_LL12_max) then
         cfail = 0
         return
      endif
      p_temp(:) = qext(:,3) + qext(:,5)
      temp = sqrt(p_temp(2)**2 + p_temp(3)**2)
      if (temp .lt. pT_LL13_min .or. temp .ge. pT_LL13_max) then
         cfail = 0
         return
      endif
***]pTV cut

***[pTmiss
      if (lWW) then
       p_miss(:) = qext(:,4) + qext(:,6)
      else if (lWZ) then
       p_miss(:) = qext(:,4)
      endif
      if (lWW.or.lWZ) then
       pT_miss = sqrt(p_miss(2)**2 + p_miss(3)**2)
       if ( pT_miss .lt. pT_miss_min .or. pT_miss .ge. pT_miss_max) then
        cfail = 0
        return
       endif      
      endif
***]pTmiss

      return
      end subroutine DoCut_DiBoson
***********************************************************************
	subroutine Get_charged_leptons(Next,qext,pt,py,phi,q_clep)
! input: qext(4,Next) with the order: q_1 q_2 > l_3 + l_4 + ... + jet.
! output: q_clep(:,Nlep_charged)
	use global_vars
	implicit none
	integer Next, i, j, mu, ic
	double precision qext(4,Next),q_clep(8,Nlep_charged)
	double precision pt(Next), py(Next), phi(Next)
	q_clep(:,:) = 0d0
	ic = 0
	do i = 3, Next
	IF (.not.lNLO_FIN_IND) THEN
	 if ((abs(leg(i))-11)*(abs(leg(i))-13)*(abs(leg(i))-15).eq.0) then
	  ic = ic + 1
	  q_clep(8,ic) = leg(i) ! particle identification number
	  do mu = 1, 4
	   q_clep(mu,ic) = qext(mu,i)
	  enddo
	   q_clep(5,ic) = pt(i)
	   q_clep(6,ic) = py(i)
	   q_clep(7,ic) = phi(i)
	 endif
	ELSE
	! the order in leg(:): q1 q2 > j3 l4 l5 ...
	 if ((abs(leg(i))-11)*(abs(leg(i))-13)*(abs(leg(i))-15).eq.0) then
	  j = i - 1 ! move the j3 to the end
	  ic = ic + 1
	  q_clep(8,ic) = leg(i)
	  do mu = 1, 4
	  q_clep(mu,ic) = qext(mu,j)
	  enddo
	  q_clep(5,ic) = pt(j)
	  q_clep(6,ic) = py(j)
	  q_clep(7,ic) = phi(j)
	 endif
	ENDIF
	enddo
        end subroutine Get_charged_leptons
***********************************************************************
	subroutine PhotLepComb(Nlep,DRcut,phot,clep_bared,clep_dressed)
	use global_cuts
	implicit none
        integer Nlep, i, j, ilep_nearest
        double precision DRcut
        double precision phot(7),clep_bared(7,Nlep),clep_dressed(7,Nlep)
        double precision DR(Nlep), DS(Nlep), DSmin
        double precision pt, py, phi, temp, temp1

c        print*,"from PhotLepComb(): Nlep, DRcut =", Nlep, DRcut

        do j = 1, Nlep
         do i = 1, 7
          clep_dressed(i,j) = clep_bared(i,j) ! initial values
         enddo
        enddo

        do i = 1, Nlep
        call rsep_simple(phot(6),phot(7),clep_bared(6,i),clep_bared(7,i),DR(i))
        enddo
        
        temp1 = phot(5)**p_recomb
        do i = 1, Nlep
         temp = min(temp1,clep_bared(5,i)**p_recomb)
         DS(i) = DR(i)*temp
        enddo

        DSmin = DS(1)
        ilep_nearest = 1
        do i = 2, Nlep
         if (DS(i) .lt. DSmin) then
          DSmin = DS(i)
          ilep_nearest = i
         endif
        enddo

! recombine
        if (DR(ilep_nearest) .le. DRcut) then
         do i = 1, 4
          clep_dressed(i,ilep_nearest) = clep_bared(i,ilep_nearest) + phot(i)
         enddo
         call lego( clep_dressed(1,ilep_nearest), pt, py, phi )
         clep_dressed(5,ilep_nearest) = pt
         clep_dressed(6,ilep_nearest) = py
         clep_dressed(7,ilep_nearest) = phi
        endif

        end
***********************************************************************
	subroutine DefineJets(jetp, Npar, qpar, Rcut, Njet, qjet)
C       kT-scheme as defined in
C       M.Seymour, Nucl.Phys.B513:269-300,1998 [hep-ph/9707338]
C       anti-KT as in Cacciari, Salam and Soyez: arXiv:0802.1189

C Input: jetp = 1 for kt, -1d0 for anti-kt,
C Input: Npar (number of partons), qpar(4, Npar) (momentum of partons), Rcut(jet-algorithm parameter).
C Output: Njet (number of jets), qjet(4, Npar) (inverse-pT-order momenta of jets). 
	implicit none
c        integer max_jets
c        parameter (max_jets = 3) ! default = max Npar
	integer Npar, Njet
        double precision qpar(4, Npar), qjet(4, Npar), qjet0(4, Npar), Rcut
        double precision jetp
C local variables
	integer i

***[ 1 parton
	IF( Npar .eq. 1 ) then
          Njet = 1
          do i=1,4
          qjet(i, 1) = qpar(i, 1)
          enddo
          return
	ENDIF
***] 1 parton

***[ 2 partons
	IF( Npar .eq. 2 ) then
         call DefineJets2(jetp, qpar(1,1), Rcut, Njet, qjet0(1,1))
C inverse-pT-order momenta of jets
         call OrderParticlesPt(Njet,qjet0(1, 1),qjet(1,1))
         return
	ENDIF
***] 2 partons

***[ 3 partons
	IF( Npar .eq. 3 ) then
         call DefineJets3(jetp, qpar(1,1), Rcut, Njet, qjet0(1,1))
C inverse-pT-order momenta of jets
         call OrderParticlesPt(Njet,qjet0(1, 1),qjet(1,1))
         return
	ENDIF
***] 3 partons

	IF( Npar .gt. 3 ) then
         print*,'jet algorithm todo: N > 3'
	 stop
	ENDIF

	end
***********************************************************************
	subroutine DefineJets2(jetp, qpar, Rcut, Njet, qjet)
C       kT-scheme as defined in
C       M.Seymour, Nucl.Phys.B513:269-300,1998 [hep-ph/9707338]
C       anti-KT as in Cacciari, Salam and Soyez: arXiv:0802.1189

C For 2 partons, which can be combined partons (from 3-parton events), hence can be massive.

C Input: jetp(1d0 for kt or -1d0 anti-kt), qpar(4, Npar) (momentum of partons), Rcut(jet-algorithm parameter).
C Output: Njet (number of jets), qjet(4, Npar). 
	implicit none
	integer Npar, Njet
        parameter (Npar = 2)
        double precision qpar(4, Npar), qjet(4, Npar), Rcut, jetp
C local variables
	integer i
        double precision DR

***[ 2 partons
         call GetDisPRapPhi(qpar(1,1), qpar(1,2), DR)
         IF ( DR .lt. Rcut ) then 
          Njet = 1
          do i=1,4
          qjet(i, 1) = qpar(i, 1) + qpar(i, 2)
          enddo
	 ELSE
          Njet = 2
          do i=1,4
           qjet(i, 1) = qpar(i, 1)
           qjet(i, 2) = qpar(i, 2)
          enddo
	 ENDIF
***] 2 partons
        end
***********************************************************************
	subroutine DefineJets3(jetp, qpar, Rcut, Njet, qjet)
C       kT-scheme as defined in
C       M.Seymour, Nucl.Phys.B513:269-300,1998 [hep-ph/9707338]
C       anti-KT as in Cacciari, Salam and Soyez: arXiv:0802.1189

C For 3 partons, which can be combined partons (from 4-parton events), hence can be massive.
C Input: jetp(1d0 for kt or -1d0 for anti-kt), qpar(4, Npar) (momentum of partons), Rcut(jet-algorithm parameter).
C Output: Njet (number of jets), qjet(4, Npar). 

C agree with VBFNLO jetdefinition.F !
	implicit none
	integer Npar, Njet
        parameter (Npar = 3)

        double precision qpar(4, Npar), qjet(4, Npar), Rcut, jetp
C local variables
	integer i
        double precision qpar1(4, Npar)
        double precision DR(Npar), pT2(Npar), dd(Npar), pT2min, pT2ddmin

***[ 3 partons
         call GetDisPRapPhi(qpar(1,1), qpar(1,2), DR(1))
         call GetDisPRapPhi(qpar(1,1), qpar(1,3), DR(2))
         call GetDisPRapPhi(qpar(1,2), qpar(1,3), DR(3))
         do i = 1, Npar
	 pT2(i) = qpar(2,i)**2 + qpar(3,i)**2
         pT2(i) = pT2(i)**jetp ! generalized kt scheme
         enddo

         dd(1) = min(pT2(1),pT2(2))*DR(1)**2/Rcut**2
         dd(2) = min(pT2(1),pT2(3))*DR(2)**2/Rcut**2
         dd(3) = min(pT2(2),pT2(3))*DR(3)**2/Rcut**2
         pT2min = min(pT2(1),pT2(2),pT2(3))
         pT2ddmin = min(dd(1),dd(2),dd(3))

         IF (pT2ddmin .lt. pT2min) THEN
          if (pT2ddmin .eq. dd(1)) then
           do i=1,4
           qpar1(i, 1) = qpar(i, 1) + qpar(i, 2)
           qpar1(i, 2) = qpar(i, 3)
           enddo
           call DefineJets2(jetp, qpar1(1,1), Rcut, Njet, qjet(1,1))
          elseif (pT2ddmin .eq. dd(2)) then
           do i=1,4
           qpar1(i, 1) = qpar(i, 1) + qpar(i, 3)
           qpar1(i, 2) = qpar(i, 2)
           enddo
           call DefineJets2(jetp, qpar1(1,1), Rcut, Njet, qjet(1,1))
          else
           do i=1,4
           qpar1(i, 1) = qpar(i, 2) + qpar(i, 3)
           qpar1(i, 2) = qpar(i, 1)
           enddo
           call DefineJets2(jetp, qpar1(1,1), Rcut, Njet, qjet(1,1))
          endif

         ELSE
          if (pT2min .eq. pT2(1)) then
           do i=1,4
           qpar1(i, 1) = qpar(i, 2)
           qpar1(i, 2) = qpar(i, 3)
           enddo
           call DefineJets2(jetp, qpar1(1,1), Rcut, Njet, qjet(1,1))
           Njet = Njet + 1
           do i=1,4
           qjet(i, Njet) = qpar(i, 1)
           enddo
          elseif (pT2min .eq. pT2(2)) then
           do i=1,4
           qpar1(i, 1) = qpar(i, 1)
           qpar1(i, 2) = qpar(i, 3)
           enddo
           call DefineJets2(jetp, qpar1(1,1), Rcut, Njet, qjet(1,1))
           Njet = Njet + 1
           do i=1,4
           qjet(i, Njet) = qpar(i, 2)
           enddo
          else
           do i=1,4
           qpar1(i, 1) = qpar(i, 1)
           qpar1(i, 2) = qpar(i, 2)
           enddo
           call DefineJets2(jetp, qpar1(1,1), Rcut, Njet, qjet(1,1))
           Njet = Njet + 1
           do i=1,4
           qjet(i, Njet) = qpar(i, 3)
           enddo
          endif
         ENDIF
***] 3 partons
	end
***********************************************************************
	subroutine Leptonc(Nlep, qlep, Rcut, qlep1, fail)
C For 4 leptons at most.
C Input: Nlep (number of leptons), qlep(4, Nlep) (momentum of leptons), Rcut(lepton-separation parameter).
C Output: qlep(4, Nlep) (inverse-pT-order momenta of leptons), fail (a flag). 
	implicit none
	integer Nlep, i, j, fail

        double precision qlep(4, Nlep), qlep1(4, Nlep), Rcut
C local variables
        double precision DR(15) ! it is Nlep*(Nlep-1)/2 = 15 for 6 leptons
        double precision pT2(Nlep)
	integer ii(Nlep), jj(Nlep)

C default value
	fail = 1

          do j=1,Nlep
          do i=1,4
          qlep1(i, j) = qlep(i, j)
          enddo
          enddo

***[ 1 lepton
	IF( Nlep .eq. 1 ) return
***] 1 lepton

***[ 2 leptons
	IF( Nlep .eq. 2 ) THEN
         call GetDisRapPhi(qlep(1,1), qlep(1,2), DR(1))
         if ( DR(1) .lt. Rcut ) then 
          fail = 0
	 else
C inverse-pT-order momenta of leptons
         call OrderParticlesPt(Nlep,qlep(1, 1),qlep1(1,1))
	 endif
	 return
	ENDIF
***] 2 leptons

***[ 3 leptons
	IF( Nlep .eq. 3 ) THEN
         call GetDisRapPhi(qlep(1,1), qlep(1,2), DR(1))
         call GetDisRapPhi(qlep(1,1), qlep(1,3), DR(2))
         call GetDisRapPhi(qlep(1,2), qlep(1,3), DR(3))

         if (DR(1) .ge. Rcut .and. DR(2) .ge. Rcut 
     & .and. DR(3) .ge. Rcut) then
C inverse-pT-order momenta of leptons
         call OrderParticlesPt(Nlep,qlep(1, 1),qlep1(1,1))
         else
          fail = 0
         endif
         return
	ENDIF
***] 3 leptons

***[ 4 leptons
	IF( Nlep .eq. 4 ) THEN
         call GetDisRapPhi(qlep(1,1), qlep(1,2), DR(1))
         call GetDisRapPhi(qlep(1,1), qlep(1,3), DR(2))
         call GetDisRapPhi(qlep(1,1), qlep(1,4), DR(3))
         call GetDisRapPhi(qlep(1,2), qlep(1,3), DR(4))
         call GetDisRapPhi(qlep(1,2), qlep(1,4), DR(5))
         call GetDisRapPhi(qlep(1,3), qlep(1,4), DR(6))

         if (DR(1) .ge. Rcut .and. DR(2) .ge. Rcut 
     & .and. DR(3) .ge. Rcut .and. DR(4) .ge. Rcut 
     & .and. DR(5) .ge. Rcut .and. DR(6) .ge. Rcut ) then
C inverse-pT-order momenta of leptons
         call OrderParticlesPt(Nlep,qlep(1, 1),qlep1(1,1))
         else
          fail = 0
         endif
         return
	ENDIF
***] 4 leptons
	end
***********************************************************************
	subroutine RJetJet(Njet, qjet, Rcut, fail)
	implicit none
	integer Njet, i, j, fail

        double precision qjet(4, Njet), Rcut
C local variables
        double precision DR

C default value
	fail = 1

        do i = 2, Njet
         do j = 1, i - 1
         call GetDisPRapPhi(qjet(1,j), qjet(1,i), DR)
         if ( DR .lt. Rcut ) then
	  fail = 0
          return
         endif
         enddo
        enddo

	end
***********************************************************************
	subroutine JetLep(Njet, Nlep, qjet, qlep, Rcut, fail)
	implicit none
	integer Njet, Nlep, i, j, fail

        double precision qjet(4, Njet), qlep(4, Nlep), Rcut
C local variables
        double precision DR

C default value
	fail = 1

        do i = 1, Nlep
         do j = 1, Njet
         call GetDisPRapPhi(qjet(1,j), qlep(1,i), DR)
         if ( DR .lt. Rcut ) then
	  fail = 0
          return
         endif
         enddo
        enddo

	end
***********************************************************************
	subroutine CutpTmiss(qq, pTmin, fail)
C Input: qq(4) (a momentum), pTmin(pT cut).
C Output: fail (a flag: 1 OK, 0 Not OK). 
	implicit none
	integer fail

        double precision qq(4), pTmin
C local variables
        double precision xx

	xx = sqrt(qq(2)**2 + qq(3)**2)

        if ( xx .ge. pTmin ) then
	 fail = 1
        else
	 fail = 0
        endif

	end
***********************************************************************
	subroutine CutMLL(kk, qq, Mmin, fail)
C Input: qq(4) (a momentum), Mmin(invariant mass cut).
C Output: fail (a flag: 1 OK, 0 Not OK). 
	implicit none
	integer fail, i

        double precision kk(4), qq(4), pp(4), Mmin
C local variables
        double precision xx

        do i = 1, 4
         pp(i) = kk(i) + qq(i)
        enddo
        xx = sqrt(pp(1)**2 - pp(2)**2 - pp(3)**2 - pp(4)**2)

        if ( xx .lt. Mmin ) then 
         fail = 0
        else
         fail = 1
        endif

	end
***********************************************************************
	subroutine OrderParticlesPt(Njet, qjet0, qjet)
	implicit none
	integer Njet
        double precision qjet(4, Njet), qjet0(4, Njet)
C local variables
        double precision pT2(Njet)
	integer i, j
        integer ii(Njet), jj(Njet)

          do i = 1, Njet
	  pT2(i) = qjet0(2,i)**2 + qjet0(3,i)**2
          ii(i) = i
          enddo
          call OrderParticles(Njet,pT2(1),ii(1),jj(1))
          do j=1,Njet
          do i=1,4
          qjet(i, j) = qjet0(i, jj(j))
          enddo
          enddo
        End
***********************************************************************
	subroutine OrderParticles(NN,pT2,ii,jj)
C Input: NN (number of particles), pT2(NN) (a momentum parameter), ii(particle indices).
C Output: jj (inverse-pT2-order particle indices).
	implicit none
	integer NN, i, j, ii(NN), jj(NN), temp
        double precision pT2(NN)

        do i = 1, NN
         jj(i) = ii(i)
        enddo

        IF ( NN .eq. 1 ) THEN
         return
        ELSE
        do i = 1, NN - 1
         do j = i + 1, NN 
         if ( pT2(jj(i)) .lt. pT2(jj(j)) ) then
          temp = jj(i)
          jj(i) = jj(j)
          jj(j) = temp
         endif
         enddo
        enddo
        ENDIF
        end
***********************************************************************
	subroutine GetDisRapPhi(p1, p2, DR)
C Note: Delta_R is a measure of closeness of 2 particles. 
C This means Delta_Phi must be in [-pi, pi]. 
	implicit none
        double precision p1(4), p2(4), DR, temp
	double precision pi
	parameter (pi = 3.1415926535897932384626433832795029D0)

       DR = (p1(1)+p1(4))*((p2(1)-p2(4)))/
     & ((p1(1)-p1(4))*(p2(1)+p2(4)))
       DR = log(abs(DR))/2

       temp = abs(atan2( p1(3), p1(2) ) - atan2( p2(3), p2(2) ))

       if ( temp .gt. pi ) temp = 2*pi - temp
       DR = sqrt(DR**2 + temp**2)
       end
***********************************************************************
	subroutine GetDisPRapPhi(p1, p2, DR)
C Note: Delta_R is a measure of closeness of 2 particles. 
C This means Delta_Phi must be in [-pi, pi]. 
	implicit none
        double precision p1(4), p2(4), DR, temp, pp1, pp2
	double precision pi
	parameter (pi = 3.1415926535897932384626433832795029D0)

       pp1 = sqrt(p1(2)**2 + p1(3)**2 + p1(4)**2)
       pp2 = sqrt(p2(2)**2 + p2(3)**2 + p2(4)**2)

       DR = (pp1+p1(4))*((pp2-p2(4)))/
     & ((pp1-p1(4))*(pp2+p2(4)))
       DR = log(abs(DR))/2

       temp = abs(atan2( p1(3), p1(2) ) - atan2( p2(3), p2(2) ))

       if ( temp .gt. pi ) temp = 2*pi - temp
       DR = sqrt(DR**2 + temp**2)
       end
***********************************************************************
	subroutine rsep_simple(x1,phi1,x2,phi2,rsep)
C Note: Delta_R is a measure of closeness of 2 particles. 
C This means Delta_Phi must be in [-pi, pi]. 
	implicit none
        double precision x1, x2 ! can be rapidity or pseudo-rapidity
        double precision phi1, phi2, rsep
	double precision pi
	parameter (pi = 3.1415926535897932384626433832795029D0)
        rsep = abs(phi1 - phi2)
        if ( rsep .gt. pi ) rsep = 2*pi - rsep
        rsep = sqrt(rsep**2 + (x1 - x2)**2)
        end
***********************************************************************
	subroutine CutInc(qq, pTmin, pymax, fail)
C Input: qq(4) (a momentum), pTmin(pT cut), pymax(pseudorapidity cut).
C Output: fail (a flag: 1 OK, 0 Not OK). 
	implicit none
	integer fail

        double precision qq(4), pTmin, pymax
C local variables
        double precision xx(2)

	xx(1)=qq(2)**2 + qq(3)**2
        xx(2)=xx(1)+qq(4)**2
        xx(1)=sqrt(xx(1))
	xx(2)=sqrt(xx(2))
	xx(2)=log((xx(2)+qq(4))/(xx(2)-qq(4)))/2D0

        if ( xx(1) .ge. pTmin .and. abs(xx(2)) .le. pymax ) then
	 fail = 1
        else
	 fail = 0
        endif

	end
***********************************************************************
***[codes taken from VBFNLO-2.7.1
c------------------------------------------------------------------------------
c This code is taken from VBFNLO-2.7.1
c	This subroutine converts a four-momentum into
c	pT, pseudorapidity and phi variables.
c
      subroutine lego(
     &              p,		!in:  particle four-momentum
     &              pt,		!out:  corresponding pT (transverse momentum)
     &              y,		!out:  corresponding eta (pseudorapidity)
     &              phi		!out:  corresponding phi (azimuthal angle)
     &               )
      implicit none
c
c declare input/output variables
c
      real*8 p(0:3), pt, y, phi
c
c declare local variables
c
      real*8 r_zero, r_half, r_pi, r_large
      parameter( r_zero=0.0d0, r_half=0.5d0 )
      parameter( r_pi=3.141592653589793238d0 )
      parameter( r_large=1.0d+3 )
c
      real*8 pt2, pabs, arg
c
      pt2 = p(1)**2 + p(2)**2
      if ( pt2 .gt. r_zero ) then
         pt = sqrt( pt2 )
         pabs = sqrt(pt2 + p(3)**2)
c
         if ( p(3) .ge. r_zero ) then
            arg = (pabs + p(3)) / pt
         else
            arg = pt / (pabs - p(3))
         end if
         y = log( arg )
         phi = atan2( p(2), p(1) )
c
      else
         if ( p(3) .gt. r_zero ) then
            pt = r_zero
            y = r_large
            phi = r_zero
         else if ( p(3) .lt. r_zero ) then
            pt = r_zero
            y = -r_large
            phi = r_zero
         else
            pt = r_zero
            y = r_zero
            phi = r_zero
         end if
      end if
c
      return
      end
c----------------------------------------------------------------------------
***]codes taken from VBFNLO-2.7.1


