! print out important run information
	subroutine ModelDigest
	use global_vars
	implicit none

#include "model.h"
	print*, "============================"
	print*, "Reference information:"
	print*, "----------------------------"
	print*, "This_program: MulBos, Duc Ninh Le, et al, to be published"
	print*, "----------------------------"
	print*, "Reference for external modules:"
	print*, "Monte_Carlo_Integrator: BASES, Comput.Phys.Commun. 88 (1995) 309"
	print*, "Phase_space_helper: VBFNLO, Eur.Phys.J.C 84 (2024) 1003"
	print*, "Amplitude_generator_1: FeynArts, Comput.Phys.Commun. 140 (2001) 418"
	print*, "Amplitude_generator_2: FormCalc, Comput.Phys.Commun. 118 (1999) 153"
	if (lVIRT.or.lGLUGLU) then
	print*, "Loop_integral_main: LoopInts, Duc Ninh Le, et al, to be published"
	print*, "Loop_integral_helper: LoopTools, Comput.Phys.Commun. 118 (1999) 153"
	endif
	print*, "============================"
	print*, "Process information:"
	print*, "----------------------------"
	print*, "Process_ID:", proc_id
	print*, "Contribution_ID:", proc_cont
	if (proc_nlegs.eq.6) then
	 write(*,101) 'Process:',leg(1),leg(2),leg(3),leg(4),leg(5),leg(6)
	else if (proc_nlegs.eq.7) then
	 write(*,102) 'Process:',leg(1),leg(2),leg(3),leg(4),leg(5),leg(6),leg(7)
	else if (proc_nlegs.eq.8) then
	 write(*,103) 'Process:',leg(1),leg(2),leg(3),leg(4),leg(5),leg(6),leg(7),leg(8)
	else
	 print*,"invalid value of proc_nlegs:",proc_nlegs
	endif
101     format (A,I5,2X,I5,2X,I5,2X,I5,2X,I5,2X,I5)
102     format (A,I5,2X,I5,2X,I5,2X,I5,2X,I5,2X,I5,2X,I5)
103     format (A,I5,2X,I5,2X,I5,2X,I5,2X,I5,2X,I5,2X,I5,2X,I5)
	print*, "Number_of_legs:", proc_nlegs
	print*, "Energy_Unit: GeV"
	print*, "Collision_Energy:", ecm
	print*, "============================"
	print*, "Reference frame information:"
       if (frame_polar.eq.1) then
	print*, "Reference frame: V...V cms frame"
       else if (frame_polar.eq.2) then
	print*, "Reference frame: Lab frame"
       else if (frame_polar.eq.3) then
	print*, "Reference frame: Partonic frame"
       else
	print*, "Unidentified frame"
       endif
        print*, "Polarization_ID:", polar_case
        print*, "OS_mapping_ID:", osmap_case
	print*, "============================"
	print*, "SM coupling parameters:"
	print*, "----------------------------"
	print*, "1/Alfa:", 1/Alfa
	print*, "Alfa:", Alfa
	print*, "Alfas(scaleR):", Alfas
	print*, "Alfas(MZ):", AlfasMZ
	print*, "GF:", GF
	print*, "SW2:", SW2
	print*, "----------------------------"
	print*, "Gauge-boson masses and widths:"
	print*, "MW:", MW
	print*, "MZ:", MZ
	print*, "GamW:", GamW
	print*, "GamZ:", GamZ
	print*, "----------------------------"
	print*, "Fermion masses and widths:"
	print*, "ME:", ME
	print*, "MM:", MM
	print*, "ML:", ML
	print*, "MU:", MU
	print*, "MC:", MC
	print*, "MT:", MT
	print*, "MD:", MD
	print*, "MS:", MS
	print*, "MB:", MB
	print*, "GamTq:", GamTq
	print*, "----------------------------"
	print*, "Higgs mass and width:"
	print*, "MH:", MH
	print*, "GamH:", GamH
	print*, "============================"
       if (lDimReg) then
	print*, "Use dimensional regularization for IR"
       else
	print*, "Use mass regularization for IR"
       endif
	print*, "============================"
	print*, "Scale info:"
	if (lscaleF_dyn) then
	 print*, "Factorization scale is dynamical"
	else
	 print*, "Factorization scale is fixed"
	endif
	if (lscaleR_dyn) then
	 print*, "Renormalization scale is dynamical"
	else
	 print*, "Renormalization scale is fixed"
	endif
	print*, "scale_muf_case:", scale_muf_case
	print*, "scale_mur_case:", scale_mur_case
	print*, "scale_facF:", scale_facF
	print*, "scale_facR:", scale_facR	
	print*, "factorization_scale:", scale
	print*, "renormalization_scale:", scaleren
	print*, "============================"
	if (lGLUGLU) then
	 print*, "Technical cuts:"
	 print*, "eps_GramDet4:", eps_GramDet4
	 print*, "============================"
	endif
	print*, "PDF info:"
	print*, "PDF_set:", trim(pdf_set)
	print*, "PDF_member:", pdf_mem
	if (lNLO_EW) then
	 if (lPDF_SCHEME_DIS) then
	  print*,"use DIS scheme for PDF counter term (NLO EW)"
	 else
	  print*,"use MSbar scheme for PDF counter term (NLO EW)"
	 endif
	endif
	if (lNLO_QCD) then
	 if (lPDF_SCHEME_DIS) then
	  print*,"use DIS scheme for PDF counter term (NLO QCD)"
	 else
	  print*,"use MSbar scheme for PDF counter term (NLO QCD)"
	 endif
	endif
	IF (lZZ.and.Ndown.eq.2.and.Nbottom.eq.0) then
	 if (lFLV5) then
	  print*, "b-induced process included in this contribution"
	 else
	  print*, "b-induced process NOT included in this contribution"
	 endif
	ENDIF
	if (lBproc4FS) then
	print*, "This is a bottom-induced process in 4FS"
	endif
	print*, "Monte_Carlo_random_seed:", seed
	print*, "============================"
	end
	
	
