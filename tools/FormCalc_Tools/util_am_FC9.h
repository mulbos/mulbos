* prototypes for the functions in util.a
* this file is part of FormCalc
* last modified 23 Sep 08 th

	double precision ThreeMom
	double precision SInvariant, TInvariant
	double complex Pair, Eps
	double complex SxS, SeS
	integer VxS, VeS, BxS, BeS

	external ThreeMom
	external SInvariant, TInvariant
	external Pair, Eps
	external SxS, SeS
	external VxS, VeS, BxS, BeS
	
#define k(i) (1+nvec*(i-1))
#define s(i) (3+nvec*(i-1))
#define e(i) (3+nvec*(i-1)+Hel(i))
#define ec(i) (3+nvec*(i-1)-Hel(i))
#define Spinor(i,af,d) (af*2+d+7+nvec*(i-1)+Hel(i))

c Can get the momenta, polarization vectors, spinors from VecSet_FC9():
c these parameters have to match with those in VecSet_FC9()
	integer nvec, nmom_max, minvec
	parameter (nvec = 12)
	parameter (nmom_max = 20) ! max number of particles
	parameter (minvec = -3) ! FC9 convention

	double complex Vec(2,2,minvec:nvec*nmom_max)
	common /vec0_FC9/ Vec

