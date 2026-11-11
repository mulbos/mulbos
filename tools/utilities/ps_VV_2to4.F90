!*************************************************************************  
! This file helps to calculate global momenta.
!*************************************************************************  
	subroutine ps_VV_2to4(ndim, x, fac, iflag)
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
        double precision W, X1, X2, Y, m_4l
! For 1 final jets and 2 gauge bosons (4 fermions)
        double precision pjet(0:3,1),Vint(0:3,2), rext(0:3,proc_nlegs_os), Vext(0:3,proc_nlegs)
        double precision rgam, rbetgam, mext(proc_nlegs_os), tmp1, tmp2

      double precision RM2(2), MVs(2), rmg(2), rm2min(2), rm2max(2), dot
      logical NonResonance, RESONANCE, TWOBODYDECAY, TwoToJetsPlusX
      external NonResonance, RESONANCE, TWOBODYDECAY, TwoToJetsPlusX

!***[kinematic constants
        mext(1) = 0D0
        mext(2) = 0D0
        mext(3) = mass_VV(1) ! input
        mext(4) = mass_VV(2) ! input
        mext(5) = 0D0
        mext(6) = 0D0
        mext(7) = 0D0
        mext(8) = 0D0

        rm2(1) = mass_VV(1)**2        ! mass squared input
        rm2(2) = mass_VV(2)**2        ! mass squared input
        rmg(1) = mass_VV(1)*width_VV(1)        ! mass * decay width input
        rmg(2) = mass_VV(2)*width_VV(2)        ! mass * decay width input
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
!	flux = 1D0/(sqrtSX*sqrtSX*2)
        W = 1D0

        sqrtSY = sqrt(Y)*ecm

!***[generate masses of V
      rm2min(1) = 1d-3        ! minimal/maximal value of mass squared
      rm2max(1) = ecm**2   ! for tan mapping

      rm2min(2) = 1d-3        ! input, minimal/maximal value of mass squared
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

!***]generate masses of V

! Vegas random number order: x(3), x(4), x(5)
        call TwoToTwo(x(3),RN,sqrtSY,MVs(1),MVs(2),Vext(0,1),Vext(0,2),&
     & X1,X2,Vint(0,1),Vint(0,2),W)
         if ( w.eq.0 ) then
          fac = 0
          return
         endif

!        if (w .ne. 0) then
!         print*,"hello from 2to4_VBF.F, w=",w
!         stop
!        endif

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

!***[decays of V -> 2l
! Vegas random number order: x(6), ..., x(7)
          call TwoBodyDecay0(x(6),x(7),Vint(0,1),MVs(1),&
     &                  Vext(0,3),Vext(0,4),W)
         if ( w.eq.0 ) then
          fac = 0
          return
         endif

! Vegas random number order: x(8), ..., x(9)
          call TwoBodyDecay0(x(8),x(9),Vint(0,2),MVs(2),&
     &                  Vext(0,5),Vext(0,6),W)
         if ( w.eq.0 ) then
          fac = 0
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

! Ninh fix momenta:
!      pext(1,:) = (/249.741055728578d0, 0.0d0, 0.0d0, -249.741055728578d0/)     
!      pext(2,:) = (/249.741055728577d0, 0.0d0, 0.0d0, 249.741055728577d0/)     
!      pext(3,:) = (/149.584674890381d0, -44.2028897471014d0, -4.72191407084658d0, -142.826408719025d0/)     
!      pext(4,:) = (/97.7112666786286d0, 35.4229525338694d0, 4.72191407084658d0, -90.9417923577101d0/)     
!      pext(5,:) = (/41.9100852915986d0, -26.9671961771702d0, -16.9371288030402d0, 27.2462703392653d0/)     
!      pext(6,:) = (/210.276084596546d0, 35.7471333904022d0, 16.9371288030402d0, 206.521930737470d0/)     
        return
      end subroutine
!**********************************************************************
