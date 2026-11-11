* model.h
* common blocks for the model parameters
* this file is adapted from the one in FormCalc

	double complex bogus, cI

	parameter (bogus = (-1D123, -2D123))
*	  some weird number likely to noticeably distort the final result;
*	  used for initializing arrays to check that all components
*	  have been calculated

	parameter (cI = (0D0, 1D0))
	double precision pi, hbar_c2
	parameter (pi = 3.1415926535897932384626433832795029D0)
	parameter (hbar_c2 = 3.8937966D8) ! to picobarn

* SM parameters
	double precision ME, ME2, MM, MM2, ML, ML2
	double precision MU, MU2, MC, MC2
	double precision MD, MD2, MS, MS2
	double precision Mf(4,3), Mf2(4,3)

        common /sm_para/ ME, ME2, MM, MM2, ML, ML2
        common /sm_para/ MU, MU2, MC, MC2
        common /sm_para/ MD, MD2, MS, MS2
	common /sm_para/ Mf, Mf2

	double precision MZ2, MW2
	double complex CMWsq, CMZsq, CW, CW2, SW, SW2
	double complex CMTsq
	double precision EL, Alfa, Alfa2
	double precision MB, MB2, MT2
	double precision MH2
	double precision GS, Alfas, Alfas2, AlfasMZ
	common /sm_para/ MZ2, MW2, CW, CW2, SW, SW2
	common /sm_para/ CMWsq, CMZsq, CMTsq
	common /sm_para/ EL, Alfa, Alfa2
	common /sm_para/ MT2
	common /sm_para/ MB, MB2 
	common /sm_para/ MH2
	common /sm_para/ GS, Alfas, Alfas2, AlfasMZ

	double precision GF
	double precision MW, MZ, MT, MH, MBot
	double precision GamW, GamZ
	double precision GamTq, GamH
	common /sm_para_ext/ GF, MZ, MW, MT, MH, GamZ, GamW, GamTq, GamH, MBot
