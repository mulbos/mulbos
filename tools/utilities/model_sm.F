* model_sm.F
* this file is adapted from the one of FormCalc

	subroutine ModelDefaults
	use global_vars
	implicit none
	double precision Mqscl, Mlscl
        double complex Alfa_complex

#include "model.h"

        MZ2 = MZ**2
        MW2 = MW**2
        MT2 = MT**2

c the complex masses squared
        CMWsq = MW2 - cI*MW*GamW
        CMZsq = MZ2 - cI*MZ*GamZ
        CMTsq = MT2 - cI*MT*GamTq

       if (lcomplex_mass_scheme) then
        CW2 = CMWsq/CMZsq ! complex
       else
        CW2 = MW2/MZ2
       endif
	SW2 = 1 - CW2
	CW = sqrt(CW2)
	SW = sqrt(SW2)

	Mqscl = 0
	Mlscl = Mqscl

       IF (lMassReg) THEN
       if (lVIRT.or.lGLUGLU) then
	Mqscl = 1D0
	Mlscl = Mqscl
       endif

       if (lNLO_COLL_IND.or.lNLO_COLL_RAD) then
	Mqscl = 0
	Mlscl = Mqscl
       endif
       ENDIF

      if (input_scheme.eq.3) then ! Alpha(0) scheme
***[lepton masses
	ME = .510998910D-3*Mlscl
	MM = 105.658367D-3*Mlscl
	ML = 1776.84D-3*Mlscl
***]lepton masses

***[quark masses
	MU = 66D-3*Mqscl
	MD = 66D-3*Mqscl
	MC = 1.20D0*Mqscl
	MS = 150D-3*Mqscl
	MB = 4.3D0*Mqscl
***]quark masses
      else
c Check light-fermion mass dependence: MB=1D-4,1D-3,1D-2
	MB = 1D-2

	ME = MB*Mlscl
	MM = ME
	ML = ME

	MU = MB*Mqscl
	MD = MU
	MC = MU
	MS = MU
	MB = MU
      endif
      
      if (lGLUGLU.or.lBproc4FS) then
       MB = MBot
      endif

	ME2 = ME**2
	MM2 = MM**2
	ML2 = ML**2

	MU2 = MU**2
	MD2 = MD**2
	MC2 = MC**2
	MS2 = MS**2

	MB2 = MB**2
      
      if (input_scheme.eq.1) then ! Alpha_Gmu scheme
        if (lcomplex_mass_scheme) then
	!Alfa_complex = sqrt(2D0)/pi*GF*CMWsq*(1-CMWsq/CMZsq)
	Alfa_complex = sqrt(2D0)/pi*GF*MW2*(1-MW2/MZ2)
        Alfa = abs(Alfa_complex)
        else
	Alfa = sqrt(2D0)/pi*GF*MW2*(1-MW2/MZ2) ! real coupling
        endif
      else if (input_scheme.eq.2) then ! Alpha(MZ) scheme
! Ref. hep-ph/0611353v2 (U, D, C, S), hep-ph/0105283
! The below value agrees with 128.907+-0.054. 
	Alfa = 1/128.926D0
      else if (input_scheme.eq.3) then ! Alpha(0) scheme
	Alfa = 1/137.035999679D0
      else
       print*,"invalid input scheme: input_scheme=",input_scheme
       stop
      endif
      
      	Alfa2 = Alfa**2

	EL = sqrt(4*pi*Alfa)

	end
	
************************************************************************

	subroutine ModelConstIni(fail)
	use global_vars
	implicit none
	integer fail

#include "model.h"

	double precision alphasPDF
	external alphasPDF

	fail = 0

	AlfasMZ = alphasPDF(MZ)

* Can fix at the center scale (since no counter-term for Alfas)
	Alfas = alphasPDF(scaleren)

	Alfas2 = Alfas**2
	GS = sqrt(4*pi*Alfas)

	MH2 = MH**2

	Mf(1,1) = 0
	Mf(2,1) = ME
	Mf(3,1) = MU
	Mf(4,1) = MD
	Mf(1,2) = 0
	Mf(2,2) = MM
	Mf(3,2) = MC
	Mf(4,2) = MS
	Mf(1,3) = 0
	Mf(2,3) = ML
	Mf(3,3) = MT
	Mf(4,3) = MB

	Mf2(1,1) = 0
	Mf2(2,1) = ME2
	Mf2(3,1) = MU2
	Mf2(4,1) = MD2
	Mf2(1,2) = 0
	Mf2(2,2) = MM2
	Mf2(3,2) = MC2
	Mf2(4,2) = MS2
	Mf2(1,3) = 0
	Mf2(2,3) = ML2
	Mf2(3,3) = MT2
	Mf2(4,3) = MB2	

	end
************************************************************************
	subroutine GetAlfas_dyn()
	use global_vars
	implicit none

#include "model.h"

	double precision alphasPDF
	external alphasPDF
	integer i

      if (lscaleR_dyn) then
	Alfas = alphasPDF(scaleren)
	Alfas2 = Alfas**2
	GS = sqrt(4*pi*Alfas)
      endif
      
	end
***********************************************************************





