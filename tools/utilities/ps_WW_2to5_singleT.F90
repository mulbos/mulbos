! order of integrations: search for "Vegas random number order"

	subroutine ps_WW_2to5_singleT(ndim, x, fac, iflag)
        use global_vars
        use global_mom_xpdf
!c Generate the phase space for 2 -> N jets + Massive gauge bosons
!c Input: 
!c	ndim: number of Vegas integration variables
!c	x(ndim): random-number array in (0,1)
!c	fac: some initial weight (e.g. color-spin average factor)
!c Output:
!c	fac: updated weight with flux-2Pi-UnitConversion-JacPs included
!c	X1, X2: partonic momentum fractions calculated, 
!c	all external momenta calculated

	implicit none

	integer ndim, iflag
	double precision x(ndim), RN, fac, tau1, range, sqrtSY
	integer Nj, i, j
        double precision W, X1, X2, Y, m_4l
!c For 1 final jets and 2 gauge bosons (4 fermions)
        double precision rext(0:3,proc_nlegs_os), Vext(0:3,proc_nlegs)
        double precision rgam, rbetgam, mext(proc_nlegs_os), tmp1, tmp2

      double precision RM2(2), MVs(2), rmg(2), rm2min(2), rm2max(2)
      double precision RM2_t, rmg_t, rm2min_t, rm2max_t, MTs
      double precision MVs_rem ! invariant mass of the W NOT from the top quark
      double precision MVs_top ! invariant mass of the W from the top decay
      double precision Tint(0:3) ! momentum of the top quark
      double precision Vint_rem(0:3) ! momentum of the W NOT from the top quark
      double precision Vint_top(0:3) ! momentum of the W from the top decay
      double precision Lext_top(0:3,2) ! leptons from top decay
      double precision Lext_rem(0:3,2) ! leptons NOT from top decay
      double precision dot
      logical NonResonance, RESONANCE, TWOBODYDECAY, TwoToJetsPlusX
      external NonResonance, RESONANCE, TWOBODYDECAY, TwoToJetsPlusX

!      double precision Vext_Lab(0:3,proc_nlegs)

!***[kinematic constants
        mext(1) = 0D0
        mext(2) = 0D0
        mext(3) = 0D0 ! final jet
        mext(4) = mass_VV(1) ! input
        mext(5) = mass_VV(2) ! input
        mext(6) = 0D0 ! lepton 1
        mext(7) = 0D0 ! lepton 2
        mext(8) = 0D0 ! lepton 3
        mext(9) = 0D0 ! lepton 4
        rm2(1) = mass_VV(1)**2        ! mass squared input
        rm2(2) = mass_VV(2)**2        ! mass squared input
        rmg(1) = mass_VV(1)*width_VV(1)        ! mass * decay width input
        rmg(2) = mass_VV(2)*width_VV(2)        ! mass * decay width input
        rm2_t = mass_Top**2
        rmg_t = mass_Top*width_Top
!***]kinematic constants

! number of final jets
	Nj = 1 ! no effect at LO
	RN = 0D0

        if (iflag .eq. 1) then
         Y = x(ndim)
        else
         Y = 1D0
        endif
	
! the weight: fac = fac*JacPs*flux*Tofb/(2*pi)**(3*Nout - 4)
        W = 1D0

        sqrtSY = sqrt(Y)*ecm
!***[generate masses of V
      rm2min(1) = 1d-3        ! minimal/maximal value of mass squared
      rm2max(1) = sqrtSY**2   ! for tan mapping

      rm2min(2) = 1d-3        ! input, minimal/maximal value of mass squared
      rm2max(2) = rm2max(1)   ! input, for tan mapping

! Vegas random number order: x(1), x(2)
          if (.not. Resonance(rm2(1), rmg(1), rm2min(1), rm2max(1),&
     &                      x(1), W, MVs(1))) then
        fac = 0
        return
      endif
          if (.not. Resonance(rm2(2), rmg(2), rm2min(2), rm2max(2),&
     &                      x(2), W, MVs(2))) then
        fac = 0
        return
      endif
!***]generate masses of V

!!![generate inv mass of top
      if (leg(1).eq.-5) then ! b_bar gluon process
       rm2min_t = MVs(2)
       MVs_top = MVs(2)
       MVs_rem = MVs(1)
      else ! b gluon process
       rm2min_t = MVs(1)
       MVs_top = MVs(1)
       MVs_rem = MVs(2)
      endif
!      rm2min_t = 1d-3 ! Ninh testing (slightly worse)
      rm2max_t = rm2max(1)
      ! use random number x(3):
      if (.not. Resonance(rm2_t, rmg_t, rm2min_t, rm2max_t,&
     &                      x(3), W, MTs)) then
        fac = 0
        return
      endif
!!!]generate inv mass of top

! q g > top + W 
! Vegas random number order: x(4), x(5), x(6)
        call TwoToTwo(x(4),RN,sqrtSY,MTs,MVs_rem,Vext(0,1),Vext(0,2),&
     & X1,X2,Tint(0),Vint_rem(0),W)
         if ( w.eq.0 ) then
          fac = 0
          return
         endif
         
!c decay of Top > b + W (m_b = 0)
!c Vegas random number order: rd(7), ..., rd(8)
      if( .not. TwoBodyDecay(x(7), x(8), Tint(0), MTs, 0D0, MVs_top,&
     &        Vext(0,3), Vint_top(0), W)) then
        fac = 0
        return
      endif

!***[decays of V -> 2l
!c Vegas random number order: x(9), ..., x(10)
          call TwoBodyDecay0(x(9),x(10),Vint_top(0),MVs_top,&
     &                  Lext_top(0,1),Lext_top(0,2),W)
         if ( w.eq.0 ) then
          fac = 0
          return
         endif

!c Vegas random number order: x(11), ..., x(12)
          call TwoBodyDecay0(x(11),x(12),Vint_rem(0),MVs_rem,&
     &                  Lext_rem(0,1),Lext_rem(0,2),W)
         if ( w.eq.0 ) then
          fac = 0
          return
         endif

! re-ordering:
      if (leg(1).eq.-5) then ! b_bar gluon process
       Vext(:,4) = Lext_rem(:,1)
       Vext(:,5) = Lext_rem(:,2)
       Vext(:,6) = Lext_top(:,1) ! t_bar > W- = V(2) = V_top
       Vext(:,7) = Lext_top(:,2)
      else ! b gluon process
       Vext(:,4) = Lext_top(:,1) ! t > W+ = V(1) = V_top
       Vext(:,5) = Lext_top(:,2)
       Vext(:,6) = Lext_rem(:,1)
       Vext(:,7) = Lext_rem(:,2)
      endif
!***]decays of V -> 2l

!!! output: X1, X2, Vext(0,i); i <= 7.

!C Note: here Vext(0,7) are the true qq -> j + 4l momenta

	if( X1 .lt. 1D-10 .or. X1 .gt. (1 - 1D-10) ) then
	  fac = 0
	endif

	if( X2 .lt. 1D-10 .or. X2 .gt. (1 - 1D-10) ) then
	  fac = 0
	endif

!***[ change from (x1, x2) to (tau, x)
        tau1 = X1*X2
        range = 1 - tau1
        tau_pdf = sqrt(tau1)
        x_pdf = (X2 - tau1)/range
        fac = fac*X2/(range*2*tau_pdf) ! undo the Jacobian in lumi_hadron.F
!***] change from (x1, x2) to (tau, x)
      x1_pdf = x1
      x2_pdf = x2

        IF (lDPA) THEN
        tmp1 = mext(4) + mext(5)
        m_4l = (Vext(0,4) + Vext(0,5) + Vext(0,6) + Vext(0,7))**2 &
     & - (Vext(1,4) + Vext(1,5) + Vext(1,6) + Vext(1,7))**2 &
     & - (Vext(2,4) + Vext(2,5) + Vext(2,6) + Vext(2,7))**2 &
     & - (Vext(3,4) + Vext(3,5) + Vext(3,6) + Vext(3,7))**2 
        m_4l = sqrt(m_4l)
        if (m_4l .le. tmp1) then 
          fac = 0
          return
        endif

        IF (lNLO_FIN_IND_TW) THEN
        ! for on-shell TW contribution
        tmp1 = sqrt(rm2_t) + mext(4)
        m_4l = (Vext(0,3)+Vext(0,4) + Vext(0,5) + Vext(0,6) + Vext(0,7))**2 &
     & - (Vext(1,3)+Vext(1,4) + Vext(1,5) + Vext(1,6) + Vext(1,7))**2 &
     & - (Vext(2,3)+Vext(2,4) + Vext(2,5) + Vext(2,6) + Vext(2,7))**2 &
     & - (Vext(3,3)+Vext(3,4) + Vext(3,5) + Vext(3,6) + Vext(3,7))**2 
        m_4l = sqrt(m_4l)
        if (m_4l .le. tmp1) then
         fac = 0
         return
        endif
        ENDIF
        ENDIF
        
        fac = fac*W/1D3 ! to picobarn

!C Note: Vext_Lab() are the exact momenta in the LAB frame
!      Vext_Lab = Vext

!***[ boost to the partonic CMS
        rgam = 0.5D0*(x1 + x2)/sqrt(x1*x2)

        rbetgam=-0.5D0*(x2 - x1)/sqrt(x1*x2) ! boost from Proton -> parton CMS

! to store in xsection.h
        rgam_boost = rgam
        rbetgam_boost = rbetgam

        do i = 1, proc_nlegs 
        call boost0(rgam, rbetgam,&
     &  Vext(0,i), Vext(1,i), Vext(2,i), Vext(3,i))
        enddo
! Ninh checking
!        print*,"2to5_VBF.F: after partonic CMS boost"
!        call PrintPSPMG(Vext, proc_nlegs)
!***] boost to the partonic CMS

!C Note: pext(proc_nlegs,4) are the exact momenta in the partonic CMS
      call translate_V2P(Vext(0,1),proc_nlegs) ! calculate pext(:,:)
      call Cal_mom_ID()
      
!!![input momenta
!#define FIX_MOM
#ifdef FIX_MOM
! pext(nlegs_max,4)
! from PrintPSP: true momenta (pext):
! number of legs =           7
   pext(1,1)=91.1612621576649d0
   pext(1,2)=0d0
   pext(1,3)=0d0
   pext(1,4)=pext(1,1)
   pext(2,1)=91.1612621576648d0
   pext(2,2)=0d0
   pext(2,3)=0d0
   pext(2,4)=-pext(2,1)  
   pext(3,1)=37.0517586059594d0
   pext(3,2)=-17.8739876528907d0
   pext(3,3)=-31.3459503903795d0
   pext(3,4)=8.41336884383793d0

   pext(4,1)=42.3982372368540d0
   pext(4,2)=29.3556519406349d0
   pext(4,3)=27.7858261024788d0
   pext(4,4)=12.7985971002945d0

   pext(5,1)=50.1720868372695d0
   pext(5,2)=-45.2002091034675d0
   pext(5,3)=-6.91371548801667d0
   pext(5,4)=-20.6489692905017d0

   pext(6,1)=43.7203030741737d0
   pext(6,2)=42.3341562064514d0
   pext(6,3)=10.4738397759175d0
   pext(6,4)=-3.09560971922716d0

   pext(7,1)=8.98013856107320d0
   pext(7,2)=-8.61561139072809d0
   pext(7,3)=-1.055108091235263d-015
   pext(7,4)=2.53261306559651d0

        x1 = 5.797778927368320d-002
        x2 = 3.392597160647159d-003
        x1_pdf = x1
        x2_pdf = x2
    
        rgam = 0.5D0*(x1 + x2)/sqrt(x1*x2)
        rbetgam=-0.5D0*(x2 - x1)/sqrt(x1*x2) ! boost from Proton -> parton CMS

! to store in xsection.h
        rgam_boost = rgam
        rbetgam_boost = rbetgam
#endif
!!!]input momenta

        return
	end subroutine
