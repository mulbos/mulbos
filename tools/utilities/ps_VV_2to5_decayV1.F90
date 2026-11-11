! the kinematics for a 2 -> 5 process, for W > 2l + photon decay

! order of integrations: search for "Vegas random number order"

!***********************************************************************
	subroutine ps_VV_2to5_decayV1(ndim, x, fac, iflag)
        use global_vars
        use global_mom_xpdf
!c Generate the phase space for 2 -> N jets + Massive gauge bosons
!c Input: 
!c	ndim: number of Vegas integration variables
!c	x(ndim): random-number array in (0,1)
!c	fac: some initial weight (e.g. color-spin average factor)
!c Output:
!c	fac: updated weight with flux-2Pi-UnitConversion-JacPs included
!c	X1, X2: partonic momentum fractions calculated, see common in xsection.h
!c	all external momenta calculated, see common in xsection.h

	implicit none

	integer ndim, iflag
        logical flag
	double precision x(ndim), RN, fac, tau1, range, sqrtSY
	integer Nj, i, j
        double precision W, X1, X2, Y, m_4l, ex, ey, ez
!c For 1 final jets and 2 gauge bosons (4 fermions)
        double precision Vint(0:3,3), rext(0:3,proc_nlegs_os), sext(0:3,proc_nlegs_os),&
     &  Vext(0:3,proc_nlegs), pX(0:3)
        double precision rgam, rbetgam, mext(proc_nlegs_os), tmp1, tmp2

      double precision RM2(2), MVs(3), rmg(2), rm2min(2), rm2max(2)
      double precision MVV, dot, rm2minx, rm2maxx, ptjmin
      logical NonResonance, RESONANCE, TWOBODYDECAY, TwoToJetsPlusX, &
     & NonResonance_Log, NonResonance_Pow
      external NonResonance, RESONANCE, TWOBODYDECAY, TwoToJetsPlusX,&
     & NonResonance_Log, NonResonance_Pow

! proc_nlegs = 7 
! proc_nlegs_os = 9

!!![kinematic constants
        mext(1) = 0D0
        mext(2) = 0D0
        mext(3) = mass_VV(1) ! input
        mext(4) = mass_VV(2) ! input

        mext(5) = 0D0 ! final photon
        mext(6) = 0D0 ! lepton 1
        mext(7) = 0D0 ! lepton 2
        mext(8) = 0D0 ! lepton 3
        mext(9) = 0D0 ! lepton 4
        rm2(1) = mass_VV(1)**2 ! mass squared input
        rm2(2) = mass_VV(2)**2 ! mass squared input
        rmg(1) = mass_VV(1)*width_VV(1) ! mass * decay width input
        rmg(2) = mass_VV(2)*width_VV(2) ! mass * decay width input
!!!]kinematic constants

! number of final jets
	Nj = 1 ! no effect at LO
	RN = 0D0

        if (iflag .eq. 1) then
         Y = x(ndim)
        else
         Y = 1D0
        endif
	
! the weight: fac = fac*JacPs*flux*Tofb/(2*pi)**(3*Nout - 4)
        W = 1D0 ! for the exact phase-space weight

        sqrtSY = sqrt(Y)*ecm

!!![generate masses of V
      rm2min(1) = 1d-2        ! minimal/maximal value of mass squared
      rm2max(1) = ecm**2   ! for tan mapping

      rm2min(2) = 1d-2        ! input, minimal/maximal value of mass squared
      rm2max(2) = ecm**2   ! input, for tan mapping

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

!!!]generate masses of V

! Vegas random number order: x(3), x(4), x(5)
        call TwoToTwo(x(3),RN,sqrtSY,MVs(1),MVs(2),Vext(0,1),Vext(0,2),&
     & X1,X2,Vint(0,1),Vint(0,2),W)
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

      IF (lDPA) THEN
        tmp1 = mext(3) + mext(4)
        m_4l = (Vint(0,1) + Vint(0,2))**2 &
     & - (Vint(1,1) + Vint(1,2))**2 &
     & - (Vint(2,1) + Vint(2,2))**2 &
     & - (Vint(3,1) + Vint(3,2))**2 
        m_4l = sqrt(m_4l)
        if (m_4l .le. tmp1) then 
          fac = 0
          return
        endif
      ENDIF

!!![ change from (x1, x2) to (tau, x)
        tau1 = X1*X2
        range = 1 - tau1
        tau_pdf = sqrt(tau1)
        x_pdf = (X2 - tau1)/range
        fac = fac*X2/(range*2*tau_pdf) ! undo the Jacobian in lumi_hadron.F
!!!] change from (x1, x2) to (tau, x)
      x1_pdf = x1
      x2_pdf = x2

!!![decays of V -> 2l
! Vegas random number order: x(6), x(9), ..., x(10)
          call ThreeBodyDecay0(x(6),Vint(0,1),MVs(1),&
     &                  Vext(0,3),Vext(0,4),Vext(0,7),W)
         if ( w.eq.0 ) then
          fac = 0
          return
         endif

! Vegas random number order: x(11), ..., x(12)
          call TwoBodyDecay0(x(11),x(12),Vint(0,2),MVs(2),&
     &                  Vext(0,5),Vext(0,6),W)
         if ( w.eq.0 ) then
          fac = 0
          return
         endif
!!!]decays of V -> 2l
 
!!![ boost to the partonic CMS
        rgam = 0.5D0*(x1 + x2)/sqrt(x1*x2)

        rbetgam=-0.5D0*(x2 - x1)/sqrt(x1*x2) ! boost from Proton -> parton CMS

! to store in xsection.h
        rgam_boost = rgam
        rbetgam_boost = rbetgam

        do i = 1, proc_nlegs 
        call boost0(rgam, rbetgam,&
     &  Vext(0,i), Vext(1,i), Vext(2,i), Vext(3,i))
        enddo
!!!] boost to the partonic CMS

! Note: pext(6,4) are the true qq -> 4l + j momenta
       pext = 0d0 ! initialize ALL to zero
       do i = 1,7 
         do j = 1,4
          pext(i,j) = Vext(j-1,i)
         enddo
        enddo
      call Cal_mom_ID()
        
       fac = fac*W/1D3 ! to picobarn
       return
	end subroutine
