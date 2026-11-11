***********************************************************************
** DifferentialCS computes the differential cross-section at x.
** For all integration variables (those with zero step) it factors in
** the Jacobian, too.

	subroutine DifferentialCS(ndim, x, ncomp, result)
        use, intrinsic :: ieee_arithmetic
        use global_vars
        use global_xsection
        use global_mom_xpdf
	implicit none
#include "model.h"
	integer ndim, ncomp, i, j, cfail
	double precision x(ndim), temp(ncomp), result(ncomp)
        double precision kext, temp_endp
        double precision amp(2)
        double precision s12, s_IR

	double precision FuncDPsing_fin_rad, FuncDPsing_decayW, FuncDPsing_decayZ
	external FuncDPsing_fin_rad, FuncDPsing_decayW, FuncDPsing_decayZ

	double precision FuncDPsing_coll, A0sq(10)
	external FuncDPsing_coll

	double precision fac, range, flux, fac1 
	double precision facx1, facx2, facx3
       integer c

	double precision FuncDPsing_fin_ind
	external FuncDPsing_fin_ind

	double precision flambda
	external flambda

c	double precision a, b, y

	double precision qext(proc_nlegs,4), qext_his(proc_nlegs,4)
	integer flagreset
	logical lcut_fail
	
!	print*,"hello from Integrand"

	fac = avgfac
	fac_pdf = 1D0
	lpdf_fixed = .false.
	flagreset = 0

	do c = 1, ncomp
	  result(c) = 0
	  temp(c) = 0	  
	enddo

        cfail = 1 ! initialize phase space cut
        amp = 0

***[ for the induced finite part 
      if (lNLO_FIN_IND) then

c        print*,"hello lNLO_FIN_IND =",lNLO_FIN_IND

        fac = avgfac
        if (lWW.and.abs(leg(1)).eq.5) then ! bg induced process
         call ps_WW_2to5_singleT(ndim, x, fac, 0)
        else
         call ps_VV_2to5(ndim, x, fac, 0)
        endif
        if( fac .eq. 0 ) return
!!![on-shell map and vecset
        IF (lNLO_FIN_IND_TW) THEN
         call OSmap_TW_2to5(fac)
         if( fac .eq. 0 ) return
         call Do_VecSet_DPA_TW()
        ELSE
         call OSmap_VV_2to5(fac)
         if( fac .eq. 0 ) return
         call Do_VecSet_DPA()
        ENDIF
!!!]on-shell map and vecset

	if( fac .eq. 0 ) return

        do i = 1,proc_nlegs
         do j = 1,4
          qext(i,j) = pext(i,j)
         enddo
        enddo

***[ get the whole 2to(N+1) amplitude

C for cuts and distributions (mom must be in the partonic CMS)
	call GetMomenta_NLO_FIN2(qext_his)
	! qext_his(proc_nlegs,4)

***[technical cuts
      call technical_cut_fin_ind(lcut_fail)
      if (lcut_fail) then
       badcount = badcount + 1
       result(1) = 0
       return
      endif
***]technical cuts

***[impose cuts
	call DoCut_parCM(cfail, qext_his(1,1))
	IF (lNLO_FIN_IND_TW) THEN
	 if ( cfail .eq. 0 ) return
	ELSE
	 if ( cfail .eq. 0 ) goto 13
	ENDIF
***]impose cuts

       if (lscaleFR_dyn) then
        call scaleFR(.false.) ! choose false for non-subtraction terms
        if (lNLO_QCD) then
         call GetAlfas_dyn()
        endif
       endif

        fac_pdf = 1D0
	call MultiplyPdfs(fac_pdf)
	if (lscaleFR_fixed) then
	 lpdf_fixed = .true.
	endif
	if( fac_pdf .eq. 0 ) then
	 goto 13
	endif
	call SquaredMES(result, qext, flags)
        amp(1) = result(1)

c	print*,'squared amp summed and averaged over heli and color'
c	print*,'Amp(N+1) =',amp(1)*avgfac
c	stop
***] get the whole 2to(N+1) amplitude

C check for NaN
        if (ieee_is_nan(result(1)) .or. .not. ieee_is_finite(result(1))) then
         nancount = nancount + 1

c	call PrintPSP()
c	call PrintPSP_OS()
c	print*,'squared amp summed and averaged over heli and color'
c	print*,'Amp(N+1) =',amp(1)*avgfac
c        stop

         result(1) = 0
         return
        endif

        if (lHist_Fill) then
c OS momenta here, to match with the dipole tilde kinematics
c important mpmenta: q1 q2 > l1 l2 l3 l4
        call FillHIS(result(1), fac*fac_pdf, qext_his)
        endif 

 13     continue
      if (lNLO_FIN_IND_TW) then
       result(1) = amp(1)*fac*fac_pdf
       goto 16
      else
       result(1) = amp(1)*fac*fac_pdf
      endif

! Distributions are filled in FuncDPsing_fin_ind().
      amp(2) = FuncDPsing_fin_ind(-fac)

c	print*,'Subtract(N+1) =',amp(2)
c	if (abs(amp(2)).gt.0) then
c	 stop
c	endif

	result(1) = result(1) + amp(2)

C check for NaN
        if (ieee_is_nan(result(1)) .or. .not. ieee_is_finite(result(1))) then
         nancount = nancount + 1

c	call PrintPSP()
c	call PrintPSP_OS()
c	print*,'squared amp summed and averaged over heli and color'
c	print*,'Amp(N+1) =',amp(1)*avgfac
c	print*,'Subtract(N+1) =',amp(2)*avgfac
c        stop
         result(1) = 0
         return
        endif

	goto 16
      endif
***] for the induced finite part 

!        print*,"hello lNLO_COLL 3 check=",lNLO_COLL

***[ for the dipole coll part 
      if (lNLO_COLL) then
        fac = avgfac
        
!        print*,"hello lNLO_COLL 3 case angfac=",avgfac
***[for triboson case
      if (lTriBoson) then
        call ps_VVV_2to6(ndim, x, fac, 1)
!        print*,"hello lNLO_COLL 31 case fac=",fac
        if( fac .eq. 0 ) goto 14
        call OSmap_VVV_2to6(fac)
      endif
***]for triboson case

***[for diboson case
      if (lDiBoson) then
        call ps_VV_2to4(ndim, x, fac, 1)       
        if( fac .eq. 0 ) goto 14
        call OSmap_VV_2to4(fac)
      endif
***]for diboson case

        if( fac .eq. 0 ) goto 14
        call Do_VecSet_DPA()

       do i = 1,proc_nlegs
         do j = 1,4
          qext(i,j) = pext(i,j)
         enddo
        enddo

C check for NaN
        if (ieee_is_nan(fac) .or. .not. ieee_is_finite(fac)) then
         nancount = nancount + 1

c	call PrintPSP()
c	call PrintPSP_OS()
c	print*,'from xsection.F coll 1'
c	print*,'fac =', fac
c        stop

         result(1) = 0
         return
        endif

       if (lscaleFR_dyn) then
        call scaleFR(.false.) ! choose false for non-subtraction terms
        ! no need to re-calculate Alfas as the born-amps are independent of it
       endif

	call SquaredMES(result, qext, flags)
	A0sq(1) = result(1)

      if (lWW.and.lNLO_EW) then ! photon-photon induced contribution
	call SquaredMES1(result, qext, flags)
	A0sq(2) = result(1)
      endif

c Distributions and cuts are done in FuncDPsing_coll()
c         print*,"amp coll 1 =", A0sq(1) ! NinhXXX

	call MultiplyPdfs_coll(facx1, facx2, facx3)
	if (lNLO_QCD) then
	 if (lscaleR_dyn) then
          call GetAlfas_dyn()
         endif
	endif
	result(1) = FuncDPsing_coll(A0sq, fac, facx1, facx2, facx3, 1)
c         print*,"amp coll 2 =", result(1) ! NinhXXX

C check for NaN
        if (ieee_is_nan(result(1)) .or. .not. ieee_is_finite(result(1))) then
         nancount = nancount + 1

c	call PrintPSP()
c	call PrintPSP_OS()
c	print*,'from xsection.F coll 1'
c	print*,'result(1) =',result(1), A0sq(1), fac, facx1, facx2, facx3
c        stop

         result(1) = 0
         return
        endif

 14     continue

cc#ifdef PHOTON_RAD
        fac1 = avgfac

***[for triboson case
      if (lTriBoson) then
        call ps_VVV_2to6(ndim, x, fac1, 0)
        if (fac1 .eq. 0) goto 16
        call OSmap_VVV_2to6(fac1)
      endif
***]for triboson case

***[for diboson case
      if (lDiBoson) then
        call ps_VV_2to4(ndim, x, fac1, 0)
        if (fac1 .eq. 0) goto 16
        call OSmap_VV_2to4(fac1)
      endif
***]for diboson case

        if( fac1 .eq. 0 ) goto 16
        call Do_VecSet_DPA()

       do i = 1,proc_nlegs
         do j = 1,4
          qext(i,j) = pext(i,j)
         enddo
        enddo

C check for NaN
        if (ieee_is_nan(fac1) .or. .not. ieee_is_finite(fac1)) then
         nancount = nancount + 1

c	call PrintPSP()
c	call PrintPSP_OS()
c	print*,'from xsection.F coll 0'
c	print*,'fac1 =', fac1
c        stop

         result(1) = 0
         return
        endif

       if (lscaleFR_dyn) then
        call scaleFR(.false.) ! choose false for non-subtraction terms
       endif

	temp(1) = result(1)

	call SquaredMES(result, qext, flags)
	A0sq(1) = result(1)
      if (lWW.and.lNLO_EW) then ! photon-photon induced contribution
	call SquaredMES1(result, qext, flags)
	A0sq(2) = result(1)
      endif
c Distributions and cuts are done in FuncDPsing_coll()
	call MultiplyPdfs_coll(facx1, facx2, facx3)
        if (lNLO_QCD) then
	 if (lscaleR_dyn) then
          call GetAlfas_dyn()
         endif
        endif
	result(1)=FuncDPsing_coll(A0sq, -fac1, facx1, facx2, facx3, 0)

C check for NaN
        if (ieee_is_nan(result(1)) .or. .not. ieee_is_finite(result(1))) then
         nancount = nancount + 1

c	call PrintPSP()
c	call PrintPSP_OS()
c	print*,'from xsection.F coll 0'
c	print*,'result(1) =',result(1), A0sq(1), fac1, facx1, facx2, facx3
c        stop

         result(1) = 0
         return
        endif

	do c = 1, ncomp
	  result(c) = temp(c) + result(c)	  
	enddo
cc#endif

	goto 16

      endif
***] for the dipole coll part 

        fac = avgfac

      if (lDiBoson) then
***[for diboson case
        IF (lNLO_FIN_RAD) then
         if (lDECAYV1_NLO_EW) then
          call ps_VV_2to5_decayV1(ndim, x, fac, 0)
          if( fac .eq. 0 ) return
          call OSmap_VV_2to5_decayV1(fac)
         else if (lDECAYV2_NLO_EW) then
          call ps_VV_2to5_decayV2(ndim, x, fac, 0)
          if( fac .eq. 0 ) return
          call OSmap_VV_2to5_decayV2(fac)
         else
          call ps_VV_2to5(ndim, x, fac, 0)
          if( fac .eq. 0 ) return
          call OSmap_VV_2to5(fac)
         endif
        ELSE
         if (proc_nlegs.eq.6) then
          call ps_VV_2to4(ndim, x, fac, 0)
         else if (proc_nlegs.eq.8.and.lBproc4FS) then
          call ps_WW_2to6_WWbb(ndim, x, fac, 0, MB)
         else
          print*,"VV phase space not implemented, proc_nlegs:",proc_nlegs
          stop
         endif

         if( fac .eq. 0 ) return
         if (lDPA) then
          if (proc_nlegs.eq.6) then
           call OSmap_VV_2to4(fac)
          else if (proc_nlegs.eq.8.and.lBproc4FS) then
           call OSmap_VV_2to6(fac)
          else
           print*,"OS mapping not implemented, proc_nlegs:",proc_nlegs
           stop
          endif
         endif
        ENDIF
***]for diboson case
      endif

***[for triboson case
      if (lTriBoson) then
! WWW:
!       print*,"hello TriBoson case"
       call ps_VVV_2to6(ndim, x, fac, 0)
       if( fac .eq. 0 ) return

!!![on-shell mapping
         if (lDPA) then
           call OSmap_VVV_2to6(fac)
           if( fac .eq. 0 ) return
         endif
!!!]on-shell mapping

!      print*,"passed OSmap step !"

      endif
***]for triboson case

!      print*,"moving to VecSet step:"
       
      	if( fac .eq. 0 ) return
      	if (lDPA) then
      	 if (lFC9) then
      	  !call Do_VecSet_DPA_FC9()
      	  print*,"FormCalc9 is not initialized"
      	  stop
      	 else
          if (proc_nlegs.eq.8.and.lBproc4FS) then
           call Do_VecSet_DPA_WWbb_4FS(MB)
          else
!       print*,"Do_VecSet_DPA ..."
           call Do_VecSet_DPA()
          endif
         endif
        else
         call Do_VecSet_FULL()
        endif

!      print*,"passed VecSet ... fac1 =", fac

       do i = 1,proc_nlegs
         do j = 1,4
          qext(i,j) = pext(i,j)
         enddo
        enddo
        
	if( fac .eq. 0 ) return

!      print*,"passed VecSet ... fac2 =", fac

***[ for the tree-level part 
      if (lTREE) then

!      print*,"enter Tree case, lTREE =", lTREE

      if (lDPA) then
       call GetMomenta_BORN(qext_his)
      endif

c      do i = 1, proc_nlegs
c	print*,"before cut qext =",i,qext_his(i,:),
c     & sqrt(abs(qext_his(i,1)**2-qext_his(i,2)**2-qext_his(i,3)**2-qext_his(i,4)**2))
c      enddo
c      print*,"check sum EP =",qext_his(1,:)+qext_his(2,:)-qext_his(3,:)-qext_his(4,:)-qext_his(5,:)-qext_his(6,:)-qext_his(7,:)
c     & -qext_his(8,:)

***[impose cuts
C physical cuts
      if (lDPA) then
	call DoCut_parCM(cfail, qext_his(1,1))
      else
	call DoCut_parCM(cfail, qext(1,1))
      endif
!      print*,"after cut cfail =",cfail
      if ( cfail .eq. 0 ) return
***]impose cuts

!!![Gram determinant cut for gluon-gluon process
      if (lGLUGLU) then
       call Check_GramDet_OS_1234(lcut_fail)
       if ( lcut_fail ) then
         gramcount4 = gramcount4 + 1
         return
       endif
      endif
!!!]Gram determinant cut for gluon-gluon process

c      print*,"passed cut: =", fac

       if (lscaleFR_dyn) then
        call scaleFR(.false.) ! choose false for non-subtraction terms
       endif

      if (lGLUGLU) then
       call clearcache ! free memory for loop integrals
       if (lscaleR_dyn) then
        call GetAlfas_dyn()
       endif
      endif

	call MultiplyPdfs(fac)
	call SquaredMES(result, qext, flags) 
c	print*,"amp2 result, fac=",result(1), fac
c	print*,"averaged amp2 =",result(1)/4D0/9D0
c	print*,"qext =",qext
c	stop

      IF (lHist_Fill) THEN
       if (lDPA) then
        call FillHIS(result(icomp), fac, qext_his) 
       else
        call FillHIS(result(icomp), fac, qext) 
       endif
      ENDIF

!      stop

	do c = 1, ncomp
	  result(c) = result(c)*fac
	enddo

C check for NaN
        if (ieee_is_nan(result(1)) .or. .not. ieee_is_finite(result(1))) then
         nancount = nancount + 1
         result(1) = 0
         return
        endif

	flagreset = 1
	goto 16

      endif
***] for the tree-level part 

***[ for the Virtual part 
      IF (lVIRT) THEN
C for cut and distributions
	call GetMomenta_BORN(qext_his)

***[impose cuts
C the momenta here has to be the OS-projected kinematics,
c to match with the subtraction term
	call DoCut_parCM(cfail, qext_his(1,1))
	if ( cfail .eq. 0 ) return
***]impose cuts

      call clearcache ! free memory for loop integrals

       if (lscaleFR_dyn) then
        call scaleFR(.false.) ! choose false for non-subtraction terms
       endif

!!![recalculate strong coupling
      if( lscaleR_dyn ) then
       if (lNLO_QCD) then
        call GetAlfas_dyn()
       endif
      endif
!!!]recalculate strong coupling

!!![recalculate renormalization constants if muR changes
! only needed for MS(bar) scheme
! no need for OS scheme
      if( lscaleR_dyn ) then
       if (lNLO_QCD) then
        call CalcRenConst_QCD
       endif
      endif
!!!]recalculate renormalization constants if muR changes

	call MultiplyPdfs(fac)
	call SquaredMES(result, qext, flags)
c	print*,'virt amp 1 =',result(2)*avgfac
c#define YES_CHECK_UV_IR

***[include endpoint
       if (lDECAYV1_NLO_EW) then
        if (lWZ.or.lWW) then
         call SoftFactor_decayW(result(1),temp_endp)
        else
         call SoftFactor_decayZ(result(1),temp_endp)
        endif
       else if (lDECAYV2_NLO_EW) then
        if (lWZ.or.lZZ) then
         call SoftFactor_decayZ(result(1),temp_endp)
        else
         call SoftFactor_decayW(result(1),temp_endp)
        endif
       else
	call SoftFactor(result(1),temp_endp)
       endif        
       result(2) = result(2) + temp_endp 
#ifdef YES_CHECK_UV_IR
	print*,'squared amp summed and averaged over heli and color'
	print*,'tree amp 1 =',result(1)*avgfac
	print*,'endp amp 1 =',temp_endp*avgfac
c	stop
#endif
***]include endpoint

        if (lHist_Fill) then
C the momenta here has to be the OS-projected kinematics,
c to match with the subtraction term
        call FillHIS(result(icomp), fac, qext_his) 
        endif

#ifdef YES_CHECK_UV_IR
	print*,'squared amp summed and averaged over heli and color'
	print*,'virt + endp amp =',result(2)*avgfac
	stop
#endif

	do c = 1, ncomp
	  result(c) = result(c)*fac
	enddo

C check for NaN
        if (ieee_is_nan(result(2)) .or. .not. ieee_is_finite(result(2))) then
         nancount = nancount + 1
         result(2) = 0
         return
        endif

	flagreset = 1
	goto 16

      ENDIF
***] for the Virtual part 

***[ for the radiation finite dipole
      if (lNLO_FIN_RAD) then

***[technical cuts
C collinear and soft cuts
      call technical_cut_fin_rad(lcut_fail)
       if ( lcut_fail ) then
         badcount = badcount + 1
         result(1) = 0
         return
        endif
***]technical cuts

C for cut and distributions
	call GetMomenta_NLO_FIN(qext_his)

***[impose cuts
C the momenta here has to be the OS-projected kinematics,
c to match with the subtraction term
	call DoCut_parCM(cfail, qext_his(1,1))
	if ( cfail .eq. 0 ) goto 15
***]impose cuts

!        print*,"Integrand.F: passed cuts for lNLO_FIN_RAD"

       if (lscaleFR_dyn) then
        call scaleFR(.false.) ! choose false for non-subtraction terms
        if (lNLO_QCD) then
         call GetAlfas_dyn()
        endif
       endif

        fac_pdf = 1D0
	call MultiplyPdfs(fac_pdf)
	if (lscaleFR_fixed) then
	 lpdf_fixed = .true.
	endif
	if( fac_pdf .eq. 0 ) then
	 goto 15
	endif

	call SquaredMES(result, qext, flags)
        amp(1) = result(1)*fac*fac_pdf

C check for NaN
        if (ieee_is_nan(result(1)) .or. .not. ieee_is_finite(result(1))) then
         nancount = nancount + 1
         result(1) = 0
         return
        endif
        
        if (lHist_Fill) then
        call FillHIS(result(1), fac*fac_pdf, qext_his)
        endif

!      print*,"End of first term"

 15     continue

!      print*,"Begin of subtraction term"

! Distributions are filled in FuncDPsing_xxx()
       if (lDECAYV1_NLO_EW) then
        if (lWZ.or.lWW) then
         amp(2) = FuncDPsing_decayW(-fac)
        else
         amp(2) = FuncDPsing_decayZ(-fac)
        endif
       else if (lDECAYV2_NLO_EW) then
        if (lWZ.or.lZZ) then
         amp(2) = FuncDPsing_decayZ(-fac)
        else
         amp(2) = FuncDPsing_decayW(-fac)        
        endif
       else
	amp(2) = FuncDPsing_fin_rad(-fac)
       endif

       result(1) = amp(1) + amp(2) ! minus sign included in the sub. term

#ifdef YES_CHECK_UV_IR
	print*,'full Integrand with all Jacobian (no PDF)'
	print*,'result(1) =',result(1)
	print*,'fac, fac_pdf =',fac, fac_pdf
#endif

!      print*,"end of subtraction term"

C check for NaN
        if (ieee_is_nan(result(1)) .or. .not. ieee_is_finite(result(1))) then
         nancount = nancount + 1
         result(1) = 0
         return
        endif
	goto 16
      endif
***] for the radiation finite dipole

 16     continue

        afcutcount = afcutcount + 1 ! all cuts passed, amp sucessfully calculated.

	if ( flagreset .eq. 1 ) then 
	flags = ibclr(flags, 0)
	flags1 = ibclr(flags1, 0)
	endif

#ifdef YES_CHECK_UV_IR
      if (lVIRT.or.lGLUGLU) then
	call ffexi
      endif
      stop
#endif
	end

************************************************************************
        subroutine SquaredMES(result, qext, flg)
        use global_vars
        use global_xsection
        implicit none
        double precision result(*), temp(2), qext(proc_nlegs,4)
        double precision p, m, ex, ey, ez
        integer flg, i

        call SquaredME(result, qext, 1)

#ifdef YES_CHECK_UV_IR
        call PrintPSP()
	call PrintPSP_OS()
	print*,'squared amp summed and averaged over heli and color'
	print*,'tree amp 0 =',result(1)*avgfac
	print*,'virt amp 0 =',result(2)*avgfac
#endif
	end

        subroutine SquaredMES1(result, qext, flg)
        use global_vars
        implicit none
        double precision result(*), temp(2), qext(proc_nlegs,4)
        double precision p, m, ex, ey, ez
        integer flg, i
        call P1SquaredME(result, qext, 1)
	end
***********************************************************************
