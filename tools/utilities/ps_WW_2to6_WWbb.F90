! order of integrations: search for "Vegas random number order"

	subroutine ps_WW_2to6_WWbb(ndim, x, fac, iflag, MB)
        use global_vars
        use global_mom_xpdf
!c Generate the phase space for 2 -> 4l + bb, with TTbar resonances
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
	integer i, j
        double precision W, X1, X2, Y, m_4l
!c For 1 final jets and 2 gauge bosons (4 fermions)
        double precision rext(0:3,proc_nlegs_os), Vext(0:3,proc_nlegs)
        double precision rgam, rbetgam, mext(proc_nlegs), tmp1, tmp2

      double precision RM2(2), rmg(2), rm2min(2), rm2max(2)
      double precision MVs(2) ! invariant masses of the W bosons
      double precision RM2_t(2), rmg_t(2), rm2min_t(2), rm2max_t(2)
      double precision MTs(2) ! invariant masses of the top quarks
      double precision Tint(0:3,2) ! momentum of the top quarks
      double precision Vint(0:3,2) ! momentum of the intermediate Ws
      double precision Lext(0:3,4) ! momenta of leptons
      double precision dot, MB
      logical NonResonance, RESONANCE, TWOBODYDECAY, TwoToJetsPlusX
      external NonResonance, RESONANCE, TWOBODYDECAY, TwoToJetsPlusX

!      double precision Vext_Lab(0:3,proc_nlegs)

!***[kinematic constants
        mext(1) = 0D0
        mext(2) = 0D0
        mext(3) = 0D0 ! lepton 1
        mext(4) = 0D0 ! lepton 2
        mext(5) = 0D0 ! lepton 3
        mext(6) = 0D0 ! lepton 4
        mext(7) = MB ! final b_bar
        mext(8) = MB ! final b

! top resonances:
        rm2_t(:) = mass_Top**2
        rmg_t(:) = mass_Top*width_Top
! W resonances:
        rm2(1) = mass_VV(1)**2        ! mass squared input
        rm2(2) = mass_VV(2)**2        ! mass squared input
        rmg(1) = mass_VV(1)*width_VV(1)        ! mass * decay width input
        rmg(2) = mass_VV(2)*width_VV(2)        ! mass * decay width input
!***]kinematic constants

! an uniform random number (for phi: 2*pi):
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

      !print*,"ps_WW, MVs =",MVs(:),fac,W

!!![generate inv masses of tops
       rm2min_t(1) = MVs(1)
       rm2min_t(2) = MVs(2)
       rm2max_t(1) = rm2max(1)
       rm2max_t(2) = rm2max_t(1)
      ! use random number x(3):
      !print*,"resonance top1:",rm2_t(1), rmg_t(1)
      !print*,"bounds =",rm2min_t(1), rm2max_t(1)
      !print*,"x3 =",x(3)
      if (.not. Resonance(rm2_t(1), rmg_t(1), rm2min_t(1), rm2max_t(1),&
     &  x(3), W, MTs(1))) then
        fac = 0
        return
      endif
      !print*,"ps_WW, MTs1 =",MTs(1),fac,W
      if (.not. Resonance(rm2_t(2), rmg_t(2), rm2min_t(2), rm2max_t(2),&
     &  x(4), W, MTs(2))) then
        fac = 0
        return
      endif
!!!]generate inv masses of tops
      !print*,"ps_WW, MTs2 =",MTs(2),fac,W

! q q > top + top
! Vegas random number order: x(5), x(6), x(7)
        call TwoToTwo(x(5),RN,sqrtSY,MTs(1),MTs(2),Vext(0,1),Vext(0,2),&
     & X1,X2,Tint(0,1),Tint(0,2),W)
         
       !print*,"ps_WW, Tint, w =",w
       !print*,"sqrtY =", sqrtSY
       !print*,"Vext1:",Vext(:,1)
       !print*,"Vext2:",Vext(:,2)
       !print*,"Tint1:",Tint(:,1)
       !print*,"Tint2:",Tint(:,2)
       !print*,"check:",Vext(:,1)+Vext(:,2)-Tint(:,1)-Tint(:,2)

         if ( w.eq.0 ) then
          fac = 0
          return
         endif

         
!c decay of Top > b + W (m_b > 0): bottom-1 at position 3
!c Vegas random number order: rd(8), ..., rd(9)
      if( .not. TwoBodyDecay(x(8), x(9), Tint(0,1), MTs(1), mext(7)**2, MVs(1),&
     &        Vext(0,3), Vint(0,1), W)) then
        fac = 0
        return
      endif
!c Vegas random number order: rd(10), ..., rd(11)
      if ( .not. TwoBodyDecay(x(10), x(11), Tint(0,2), MTs(2), mext(8)**2, MVs(2),&
     &  Vext(0,4), Vint(0,2), W)) then
        fac = 0
        return
      endif

!***[decays of V -> 2l
!c Vegas random number order: x(12), ..., x(13)
          call TwoBodyDecay0(x(12),x(13),Vint(0,1),MVs(1),&
     &          Lext(0,1),Lext(0,2),W)
         if ( w.eq.0 ) then
          fac = 0
          return
         endif

!c Vegas random number order: x(14), ..., x(15)
          call TwoBodyDecay0(x(14),x(15),Vint(0,2),MVs(2),&
     &          Lext(0,3),Lext(0,4),W)
         if ( w.eq.0 ) then
          fac = 0
          return
         endif

! re-ordering:
       Vext(:,7) = Vext(:,3) ! b_bar
       Vext(:,8) = Vext(:,4) ! b
       Vext(:,3) = Lext(:,1)
       Vext(:,4) = Lext(:,2)
       Vext(:,5) = Lext(:,3)
       Vext(:,6) = Lext(:,4)
!***]decays of V -> 2l

!C Note: here Vext(0,8) are the true qq -> 4l + bb momenta

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
        tmp1 = mass_VV(1) + mass_VV(2)
        m_4l = (Vext(0,4) + Vext(0,5) + Vext(0,6) + Vext(0,3))**2 &
     & - (Vext(1,4) + Vext(1,5) + Vext(1,6) + Vext(1,3))**2 &
     & - (Vext(2,4) + Vext(2,5) + Vext(2,6) + Vext(2,3))**2 &
     & - (Vext(3,4) + Vext(3,5) + Vext(3,6) + Vext(3,3))**2 
        m_4l = sqrt(m_4l)
        if (m_4l .le. tmp1) then 
          fac = 0
          return
        endif

        IF (lTTbar_OS_4FS) THEN
        ! for on-shell TTbar contribution
        tmp1 = sqrt(rm2_t(1)) + sqrt(rm2_t(2))
        m_4l = (Vext(0,3)+Vext(0,4)+Vext(0,5)+Vext(0,6)+Vext(0,7)+Vext(0,8))**2 &
     & - (Vext(1,3)+Vext(1,4)+Vext(1,5)+Vext(1,6)+Vext(1,7)+Vext(1,8))**2 &
     & - (Vext(2,3)+Vext(2,4)+Vext(2,5)+Vext(2,6)+Vext(2,7)+Vext(2,8))**2 &
     & - (Vext(3,3)+Vext(3,4)+Vext(3,5)+Vext(3,6)+Vext(3,7)+Vext(3,8))**2 
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

!store the momenta to: pext(proc_nlegs,4), exact momenta in the partonic CMS, 
! with same particle order as processes() in global_vars.
       pext = 0d0 ! initialize ALL to zero
      do i = 1, proc_nlegs
       pext(i,:) = Vext(:,i)
      enddo
      call Cal_mom_ID()

        return
	end subroutine
