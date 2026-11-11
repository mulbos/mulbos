! order of integrations: search for "Vegas random number order"

	subroutine ps_VV_2to5(ndim, x, fac, iflag)
	use, intrinsic :: ieee_arithmetic
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
        double precision Vint(0:3,2), rext(0:3,proc_nlegs_os), Vext(0:3,proc_nlegs), pX(0:3)
        double precision rgam, rbetgam, mext(proc_nlegs_os), tmp1, tmp2

      double precision RM2(2), MVs(2), rmg(2), rm2min(2), rm2max(2)
      double precision MVV, dot, rm2minx, rm2maxx, ptjmin
      logical lcheck_mom
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
!***]kinematic constants

!      print*,"ps_2to5, input:",mass_VV(1), mass_VV(2), width_VV(1), width_VV(2), ecm

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
      rm2max(2) = sqrtSY**2   ! input, for tan mapping

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

      rm2minx = (mext(4) + mext(5))**2
      rm2maxx = sqrtSY**2                 ! for intermediate particle X -> V1 V2
! cutoff for numerical purposes
      ptjmin = 20d0

!c generate mass of VV cluster
!c Vegas random number order: rd(1)
!c      if(.not. Resonance(xxm**2,xxwidth*xxm,rm2minx,rm2maxx,rd(3),W,pX(4))) return
       if(.not. NonResonance(rm2minx,rm2maxx,x(3),W,MVV)) then
        fac = 0
        return
      endif

!c generate X(=VV) + N Jets, pjet = Vext(0,3)
!c Vegas random number order: rd(4), ..., rd(6)
      if( .not. TwoToJetsPlusX(Nj,x(4),RN,sqrtSY,ptjmin,MVV,&
     &    Vext(0,1),Vext(0,2),X1,X2,pX(0),Vext(0,3),W)) then
        fac = 0
        return
      endif

!      do i = 1, 3
!      print*,"pi =",i,Vext(:,i),Vext(0,i)**2-Vext(1,i)**2-Vext(2,i)**2-Vext(3,i)**2
!      enddo
!      print*,"pX =",pX(:),pX(0)**2-pX(1)**2-pX(2)**2-pX(3)**2
!      print*,"EP check:",Vext(:,1)+Vext(:,2)-Vext(:,3)-pX(:)

!c decay of X to V + V
!c Vegas random number order: rd(7), ..., rd(8)
      if( .not. TwoBodyDecay(x(7), x(8), pX(0), MVV, MVs(1), MVs(2),&
     &        Vint(0,1), Vint(0,2), W)) then
        fac = 0
        return
      endif

!***[decays of V -> 2l
!c Vegas random number order: x(9), ..., x(10)
          call TwoBodyDecay0(x(9),x(10),Vint(0,1),MVs(1),&
     &                  Vext(0,4),Vext(0,5),W)
         if ( w.eq.0 ) then
          fac = 0
          return
         endif

!c Vegas random number order: x(11), ..., x(12)
          call TwoBodyDecay0(x(11),x(12),Vint(0,2),MVs(2),&
     &                  Vext(0,6),Vext(0,7),W)
         if ( w.eq.0 ) then
          fac = 0
          return
         endif
!***]decays of V -> 2l

!      do i = 1, 7
!      print*,"pi all =",i,Vext(:,i),Vext(0,i)**2-Vext(1,i)**2-Vext(2,i)**2-Vext(3,i)**2
!      enddo

      if (ieee_is_nan(w) .or. .not. ieee_is_finite(w)) then
       fac = 0
       return
      endif

      call check_mom(Vext(:,1),proc_nlegs,lcheck_mom)
      if (.not.lcheck_mom) then
       fac = 0
       return
      endif

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
!        print*,"ps_2to5: after partonic CMS boost"
!      do i = 1, 7
!      print*,"pi partonic =",i,Vext(:,i),Vext(0,i)**2-Vext(1,i)**2-Vext(2,i)**2-Vext(3,i)**2
!      enddo
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
