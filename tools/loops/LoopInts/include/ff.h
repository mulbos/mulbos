	RVAR pi
	CVAR cI,czero,cone,cipi,nan

	parameter(
     &          cI = (0D0, 1D0),
     &		czero = (0D0, 0D0),
     &		cone = (1D0, 0D0),
     &		pi  = 3.141592653589793238462643383279502884197D0,
     &		cipi = (0D0, pi),
     &          nan = (1D123, 1D123)
     &           )

* common
	RVAR prec,precd,mudim,delta,lambda,minmass,
     &	stab0,stab2,stab3,stab4,stab5,stab6,xrat
	RVAR maxdev
	integer iers
	integer debugkey,rankmax3,rankmax4,rankmax5,rankmax6
	integer versionkey,runmode
	common /YPPREGUL/ prec,precd,mudim,delta,lambda,minmass,
     &	stab0,stab2,stab3,stab4,stab5,stab6,xrat,maxdev,
     &  iers,debugkey,rankmax3,rankmax4,rankmax5,rankmax6,
     &  versionkey,runmode
