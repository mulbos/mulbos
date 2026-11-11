module global_virt
!double precision Divergence
double precision xrat, xrat1
double precision lambda2, lambda, Delta, mu_uv2, mu_ir2
contains
      subroutine set_virt_params()
      use global_vars
      implicit none
      double precision stab0, stab3, stab4, scale_ecm
      double precision sta03, sta04 ! quasi-collinear approximations
      integer rank3_max, rank4_max, rank5_max ! optimize runtime by choosing good max bounds

      Delta = 0D0
!      Delta = 1D3 ! check UV divergence
      if (lMassReg) then
       lambda2 = 1D0
!       lambda2 = 1D-6 ! check IR divergence
      else
       lambda2 = 0D0 ! finite part for dimensional regularization
!       lambda2 = -2D0 ! 1/eps^2 pole
!       lambda2 = -1D0 ! 1/eps pole       
      endif
      
      	if ( lambda2 .gt. 0) then
	 lambda = sqrt(lambda2) ! photon mass
	else
	 lambda = lambda2
	endif

      ! must set mu_uv = scaleren for MS(bar) renormalization scheme
      ! for OS renormalization: mu_uv can be anything positive
      ! mu_uv cancellation between loop corrections and counter terms
      !mu_uv2 = scaleren**2 ! for MS(bar) renormalization
      mu_uv2 = 1D0 ! default (OK for VVV production)
      ! mu_ir2 can be anything positive
      ! mu_ir cancellation between virtual corrections and endpoint (I-operator)
      mu_ir2 = 1D0 ! default (OK for VVV production)

      call setmudim_uv(mu_uv2) ! mudim^2 for UV
      call setmudim_ir(mu_ir2) ! mudim^2 for IR
      call setdelta(Delta)
      call setlambda(lambda2)
!      Divergence = Delta      

     ! xrat is of UV origin, occuring in LoopInts and FormCalc.
     ! xrat1: rational term in the fermion wave-function renormalization constants, 
     ! occuring in FeynArts. xrat1 is UV origin (hence = xrat) for mass regularization.
     ! xrat1 is of both UV and IR origin for dimensional regularization. 
	if (lambda2 .gt. 0) then
	 xrat = 1D0 ! UV origin
	 xrat1 = 1D0 ! UV origin
	elseif (lambda2 .eq. 0) then
	 xrat = 1D0 ! UV origin
	 xrat1 = 0D0 ! UV and IR cancellation
	else
	 xrat = 0D0 ! for pole extraction
	 xrat1 = 0D0 ! for pole extraction
	endif

      scale_ecm = (13d3/ecm)**2
      stab0 = 1D-6*scale_ecm
      sta03 = stab0
      sta04 = sta03
      stab3 = 1D-4
      if (lGLUGLU) then
       stab4 = 1D-2
      else
       stab4 = 1D-3
      endif
!!      stab5 1D-14
      call setstab0(stab0)
      call setsta03(sta03)
      call setsta04(sta04)
      call setstab3(stab3)      
      call setstab4(stab4)      

!!![which tensor ranks are needed
! rank3: 3,4,5 (possible values)
! rank4: 3,4,5 (possible values)
! rank5: fixed to 4 at the moment in LoopInts, can go to 5

      if (lDiBoson) then
       if (lGLUGLU) then
        rank4_max = 4 ! fermion loop
       else ! qq > VV cases
        rank4_max = 3 ! no fermion loop
       endif
      else ! VVV case
       if (lGLUGLU) then
        rank5_max = 5
        print*,"make sure to include 5pt-rank5 in Ecoeffs.F"
        stop
       else ! the qq > VVV cases
        rank5_max = 4
       endif
! rank4_max >= rank5_max:
       rank4_max = max(rank5_max,4) ! lowest value is 4
      endif

      rank3_max = max(rank4_max - 1, 3) ! lowest value is 3

!      call setrankmax(100, 100, 100) ! checking purpose
      call setrankmax(rank3_max, rank4_max, rank5_max) ! optimal runtime

!!!]which tensor ranks are needed

!       call setdebugkey(DEBUGKEY)	
!       call setversionkey(1*KeyBget)
!       call setversionkey(1*KeyC0)
!       call setversionkey(1*KeyD0)
!       call setversionkey(1*KeyEget)
!       call setversionkey(1*KeyFget)
!       call setversionkey(1*KeyAll)

      end subroutine
end module
