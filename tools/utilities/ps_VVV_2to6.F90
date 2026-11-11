!*************************************************************************  
! This file helps to calculate global momenta.
! Coded by Duc Ninh LE (2026)
!*************************************************************************  
	subroutine ps_VVV_2to6(ndim, x, fac, iflag)
        use global_vars
        use global_mom_xpdf
!************************************************************************************
! Input: ndim (number of integration variables), x(ndim) (random numbers in [0,1]), 
! iflag (whether a collinear term ?)
! Output: fac (updated), momenta.
!************************************************************************************
! Generate the phase space for 2 -> N jets + Massive gauge bosons
! Input: 
!	ndim: number of Vegas integration variables
!	x(ndim): random-number array in (0,1)
!	fac: some initial weight (e.g. color-spin average factor)
! Output:
!	fac: updated weight with flux-2Pi-UnitConversion-JacPs included
!	X1, X2: partonic momentum fractions calculated, see common in xsection.h
!	all external momenta calculated, see common in xsection.h
	implicit none

	integer ndim, iflag
        logical flag
	double precision x(ndim), RN, fac, tau1, range, sqrtSY 
	integer Nj, i, j
        double precision W, X1, X2, Y, m_VVV
! For resonace gauge bosons:
        double precision Vint(0:3,nlegs_res), Vext(0:3,proc_nlegs)
        double precision rgam, rbetgam, mext(proc_nlegs_os), threshold_os

      double precision RM2(nlegs_res), MVs(nlegs_res), rmg(nlegs_res)
      double precision rm2min(nlegs_res), rm2max(nlegs_res), dot
      double precision rm2min_X, rm2max_X, MXs, Vint_X(0:4) ! intermediate X > VVV
      logical NonResonance, RESONANCE, TWOBODYDECAY, ThreeBodyDecay
      external NonResonance, RESONANCE, TWOBODYDECAY, ThreeBodyDecay

!***[kinematic constants
        mext(:) = 0D0
        mext(3) = mass_VV(1) ! input
        mext(4) = mass_VV(2) ! input
        mext(5) = mass_VV(3) ! input

        rm2(1) = mass_VV(1)**2 ! mass squared input
        rm2(2) = mass_VV(2)**2 ! mass squared input
        rm2(3) = mass_VV(3)**2 ! mass squared input
        rmg(1) = mass_VV(1)*width_VV(1) ! mass * decay width input
        rmg(2) = mass_VV(2)*width_VV(2) ! mass * decay width input
        rmg(3) = mass_VV(3)*width_VV(3) ! mass * decay width input
!***]kinematic constants

!      print*,"hello from ps_VVV_2to6",ndim, iflag
!      print*,"mext(:) =",mext(:)

! number of final jets
	Nj = 1 ! no effect at LO
	RN = 0D0 ! a trivial random number for phi integration (=2pi)

        if (iflag .eq. 1) then
         Y = x(ndim)
        else
         Y = 1D0
        endif
	
! the weight: fac = fac*JacPs*flux*Tofb/(2*pi)**(3*Nout - 4)
!	flux = 1D0/(sqrtSX*sqrtSX*2)
        W = 1D0

        sqrtSY = sqrt(Y)*ecm

!***[generate masses of V
      rm2min(1) = 1d-3        ! minimal/maximal value of mass squared
      rm2max(1) = ecm**2   ! for tan mapping

      rm2min(2) = 1d-3        ! input, minimal/maximal value of mass squared
      rm2max(2) = ecm**2   ! input, for tan mapping

      rm2min(3) = 1d-3
      rm2max(3) = ecm**2

! Vegas random number order: x(1), x(2), x(3)
! Resonance() uses 1 random number:
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
          if (.not. Resonance(rm2(3), rmg(3), rm2min(3), rm2max(3),&
     &                      x(3), W, MVs(3))) then
        fac = 0
        return
      endif
!***]generate masses of V

!!![generate X = VVV
      rm2min_X = (sqrt(MVs(1)) + sqrt(MVs(2)) + sqrt(MVs(3)))**2
      rm2max_X = ecm**2
! NonResonance() uses 1 random number: x(4)
      if (.not. NonResonance(rm2min_X, rm2max_X, x(4), W, MXs)) then
       fac = 0d0
       return
      endif
!!!]generate X = VVV

! Vegas random number order: x(5)
! TwoToOne() uses 1 random number.
        call TwoToOne(x(5),sqrtSY,MXs,Vext(0,1),Vext(0,2),X1,X2,Vint_X(0),W)
         if ( w.eq.0 ) then
          fac = 0
          return
         endif

	if( X1 .lt. 1D-10 .or. X1 .gt. (1 - 1D-10) ) then
	  fac = 0
          return
	endif

	if( X2 .lt. 1D-10 .or. X2 .gt. (1 - 1D-10) ) then
	  fac = 0
          return
	endif

! invariant mass of the VVV:
      m_VVV = sqrt(MXs)

        IF (lDPA) THEN
        threshold_os = mext(3) + mext(4) + mext(5)
        if (m_VVV .le. threshold_os) then 
          fac = 0
          return
        endif
        ENDIF

!***[ change from (x1, x2) to (tau, x)
        tau1 = X1*X2
        range = 1 - tau1
        tau_pdf = sqrt(tau1)
        x_pdf = (X2 - tau1)/range
        fac = fac*X2/(range*2*tau_pdf) ! undo the Jacobian in lumi_hadron.F
        if (lNLO_COLL_IND.or.lNLO_COLL_RAD) then
         y_coll = x(ndim)
        endif
!***] change from (x1, x2) to (tau, x)
      x1_pdf = x1
      x2_pdf = x2

!        sqrtSX = sqrt(X1*X2)*sqrtSY

!!![decay of X -> VVV
!      LOGICAL FUNCTION ThreeBodyDecay(R, RN, P_in, m1sq, m2sq, m3sq, 
!     &                  ql1, ql2, ql3, w)
! ThreeBodyDecay() uses 4 random numbers [x6,x7,x8,x9] + 1 trivial random number (RN)
! NOTE: Vint_X(:) has 5 entries, including the squared invariant mass.
      Vint_X(4) = MXs ! squared invariant mass
      if (.not.ThreeBodyDecay(x(6), RN, Vint_X(0), MVs(1), MVs(2), MVs(3), &
     & Vint(0,1), Vint(0,2), Vint(0,3), W)) then
       fac = 0
       return
      endif

!      print*,"check ThreeBodyDecay:"
!      print*,"Vint(1) =",Vint(:,1)
!      print*,"Vint(2) =",Vint(:,2)
!      print*,"Vint(3) =",Vint(:,3)
!      print*,"Vint_X(:)     =", Vint_X(:)
!      print*,"sum 123 check =", Vint(:,1)+Vint(:,2)+Vint(:,3)
!      print*,"Vint_X2   =",Vint_X(0)**2 - Vint_X(1)**2 - Vint_X(2)**2 - Vint_X(3)**2
!      print*,"MXs check =",MXs,sqrt(MXs)
!      print*,"MVs(1) =",MVs(1),sqrt(MVs(1))
!      print*,"MVs(2) =",MVs(2),sqrt(MVs(2))
!      print*,"MVs(3) =",MVs(3),sqrt(MVs(3))

!!!]decay of X -> VVV

!***[decays of V -> 2l
! Vegas random number order: x10, x11
          call TwoBodyDecay0(x(10),x(11),Vint(0,1),MVs(1),&
     &                  Vext(0,3),Vext(0,4),W)
         if ( w.eq.0 ) then
          fac = 0
          return
         endif

! Vegas random number order: x12, x13
          call TwoBodyDecay0(x(12),x(13),Vint(0,2),MVs(2),&
     &                  Vext(0,5),Vext(0,6),W)

!      print*,"energy Vext(0,6) =", Vext(0,6)
!      print*,"mom V2 =", Vint(:,2)
!      print*,"MVs(2) =",MVs(2),sqrt(MVs(2))

         if ( w.eq.0 ) then
          fac = 0
          return
         endif

! Vegas random number order: x14, x15
          call TwoBodyDecay0(x(14),x(15),Vint(0,3),MVs(3),&
     &                  Vext(0,7),Vext(0,8),W)
         if (w.eq.0) then
          fac = 0d0
          return
         endif
!***]decays of V -> 2l

        fac = fac*W/1D3 ! to picobarn

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
        do i = 1, 2
        call boost0(rgam, rbetgam,&
     &  Vint(0,i), Vint(1,i), Vint(2,i), Vint(3,i))
        enddo

!***] boost to the partonic CMS

! Note: pext(6,4) are the true qq -> 4l momenta in the partonic CMS
       pext = 0d0 ! initialize ALL to zero
       do i = 1,proc_nlegs
         do j = 1,4
          pext(i,j) = Vext(j-1,i)
         enddo
       enddo
      call Cal_mom_ID()

!       do i = 1, proc_nlegs
!        print*,"mom =",i,pext(i,:)
!       enddo

! Ninh fix momenta:
!   pext(2,:)=(/ 2910.5053438602813003d0,    0.0000000000000000d0,    0.0000000000000000d0, 2910.5053438602813003d0/)
!   pext(1,:)=(/ 2910.5053438602813003d0,    0.0000000000000000d0,    0.0000000000000000d0,-2910.5053438602813003d0/)
!   pext(3,:)=(/ 1662.8663575180635235d0, -912.0810400620719065d0,  719.9727628705633151d0, 1189.4838880998804598d0/)
!   pext(4,:)=(/  539.6718475496548990d0, -257.4265203153563562d0,  245.8201149591465082d0,  405.6473354497199466d0/)
!   pext(5,:)=(/  257.6495081408880878d0, -182.8545595600546392d0,  -68.8545833706323833d0, -167.9479843350215731d0/)
!   pext(6,:)=(/ 1044.8359137905538319d0, -753.0477590851933201d0, -398.8742464911039178d0, -604.5663692014787785d0/)
!   pext(7,:)=(/ 1873.7435849741400489d0, 1694.3434852823083929d0, -396.8524365019629840d0, -694.7109613037598592d0/)
!   pext(8,:)=(/  442.2434757472627211d0,  411.0663937403678574d0, -101.2116114660105382d0, -127.9059087093406220d0/)


        return
      end subroutine
!**********************************************************************
