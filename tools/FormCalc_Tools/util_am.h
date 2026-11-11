* prototypes for the functions in util.a
* this file is part of FormCalc
* last modified 23 Sep 08 th

	double precision SInvariant, TInvariant
	external SInvariant, TInvariant

#include "vecset.h"	

****[FormCalc notation
*#define k(i) (8*i+1)
*#define s(i) (8*i+3)
*#define e(i) (8*i+3+Hel(i))
*#define ec(i) (8*i+3-Hel(i))
*#define Spinor(i,s,om) (s*2*Hel(i)+16*i+om+5)
* s = sign of Helilcity, om = 6/7 (chirality), i = particle
*#define DottedSpinor(i,s,om) (s*2*Hel(i)+16*i+om+7)

#define k(i) vec_K(1,1,i)
#define e(i) vec_E(1,1,Hel(i),i)
#define ec(i) vec_EC(1,1,Hel(i),i)
#define Spinor(i,s,om) spin(1,s*Hel(i),om,i)
#define DottedSpinor(i,s,om) dspin(1,s*Hel(i),om,i)

****]FormCalc notation
c Can get the polarization vectors from util:
c the parameter nlegs_max = 20 has to match VecSet()
*	double complex vec(2,2,8,0:20)
*	common /vectors/ vec

*	double complex kcomp(32,0:20)
*	equivalence (vec, kcomp)







