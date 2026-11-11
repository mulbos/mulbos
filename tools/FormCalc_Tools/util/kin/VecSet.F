* VecSet.F
* construct the CMS representation of a momentum and the
* three corresponding polarization vectors
* this file is part of FormCalc-6.0
* last modified by Le Duc Ninh 28 Oct 2025

* i: the index of the momentum
* m, p: mass and three-momentum of the particle
* ex,ey,ez: the unit three-vector of the momentum

* Ref0: Dittmaier, hep-ph/9805445v1
c#define QCD_GAUGE_CHECK

	subroutine VecSet(i, m, p, ex,ey,ez)
* input: i, m, p, ex, ey, ez
* output: spinors, polarization vectors, momentum matrix K
	implicit none
	integer i
	double precision m, p, ex,ey,ez
#include "vecset.h"

#define K(x,y)     vec_K(x,y,i)
#define EpsTR(x,y) vec_E(x,y, 1, i)
#define EpsL(x,y)  vec_E(x,y, 0, i)
#define EpsTL(x,y) vec_E(x,y, -1, i)

#define Spinor(hel,om,x) spin(x,hel,om,i)
#define DottedSpinor(hel,om,x) dspin(x,hel,om,i)

* i.e., spi(x,  9, i) = left-handed,  undotted, upper
*       spi(x, 10, i) = left-handed,  undotted, lower
*       spi(x, 11, i) = left-handed,  dotted,   upper
*       spi(x, 12, i) = left-handed,  dotted,   lower
*       spi(x, 13, i) = right-handed, undotted, upper
*       spi(x, 14, i) = right-handed, undotted, lower
*       spi(x, 15, i) = right-handed, dotted,   upper
*       spi(x, 16, i) = right-handed, dotted,   lower

	double precision momspec(16,0:20)
	common /momenta/ momspec

        integer icall
        data icall/0/
        save icall

* the following must match the definitions in util.h:
#define SPEC_M 1
#define SPEC_K 2
#define SPEC_E 3
#define SPEC_KT 4
#define SPEC_ET 5
#define SPEC_PRAP 6
#define SPEC_RAP 7
#define SPEC_DELTAK 8
#define SPEC_PHI 9
#define SPEC_EX 10
#define SPEC_EY 11
#define SPEC_EZ 12

	double precision m2, p2, p0, deltap, sump
	double precision sinth, onePez, oneMez
	double complex expIphi

	double precision sqrt2
	parameter (sqrt2 = 1.4142135623730950488016887242096981D0)

	double complex sn1(2), sn2(2)
	integer j1, j2

#ifdef QCD_GAUGE_CHECK
	double complex sgp(2), sgm(2), ctemp1, ctemp2
#endif

c        icall = icall + 1
c        print*,"from VecSet(), icall =", icall

	momspec(SPEC_M,i) = m
	momspec(SPEC_K,i) = p

	m2 = m**2
	if( m2 .lt. 1D-14 ) then
	  p0 = p
	  deltap = 0
	else
	  p2 = p**2
	  p0 = sqrt(p2 + m2)
	  deltap = m2/(p0 + p)
	endif

	momspec(SPEC_E,i) = p0
	momspec(SPEC_DELTAK,i) = deltap
	momspec(SPEC_EX,i) = ex
	momspec(SPEC_EY,i) = ey
	momspec(SPEC_EZ,i) = ez

	sinth = ex**2 + ey**2
	onePez = 1 + ez
	if( onePez .lt. .5D0 ) onePez = sinth/(1 - ez)
	oneMez = 1 - ez
	if( oneMez .lt. .5D0 ) oneMez = sinth/(1 + ez)

	if( sinth .lt. 1D-14 ) then
* phi is irrelevant when theta = 0
	  expIphi = 1
	  momspec(SPEC_PHI,i) = 0
	else
	  sinth = sqrt(sinth)
	  expIphi = DCMPLX(ex, ey)/sinth
	  momspec(SPEC_PHI,i) = atan2(ey, ex)
	endif

	K(1,1) = p0*onePez - deltap*ez
	K(2,2) = p0*oneMez + deltap*ez
	K(2,1) = p*DCMPLX(ex, ey)
	K(1,2) = DCONJG(K(2,1))

#define safelog(x) log(max(x, 1D-300))

	momspec(SPEC_KT,i) = p*sinth
	momspec(SPEC_ET,i) = sqrt(momspec(SPEC_KT,i)**2 + m2)
	momspec(SPEC_PRAP,i) = safelog(onePez) -
     &    safelog(sinth)
	momspec(SPEC_RAP,i) = safelog(DBLE(K(1,1))) -
     &    safelog(momspec(SPEC_ET,i))

* This macro is defined in the makefile
c#ifdef YES_MADGRAPH
c	return
c#endif

	if( m .ne. 0 ) then
	EpsL(1,1) = (p*onePez + deltap*ez)/m
	EpsL(2,2) = (p*oneMez - deltap*ez)/m
	EpsL(2,1) = p0/m*DCMPLX(ex, ey)
	EpsL(1,2) = DCONJG(EpsL(2,1))
	endif

* this is also OK for massless particles
	EpsTR(2,2) = 1/sqrt2*sinth
	EpsTL(2,2) = EpsTR(2,2)
	EpsTR(1,1) = -EpsTR(2,2)
	EpsTL(1,1) = EpsTR(1,1)
	EpsTR(2,1) = -1/sqrt2*oneMez*expIphi
	EpsTL(1,2) = DCONJG(EpsTR(2,1))
	EpsTL(2,1) = 1/sqrt2*onePez*expIphi
	EpsTR(1,2) = DCONJG(EpsTL(2,1))

* this is E^(I phi) cos(th/2) = 1/sqrt2*sqrt(1 + ez)*expIphi:
	expIphi = 1/sqrt2*sqrt(onePez)*expIphi
* this is sin(th/2):
	sinth = 1/sqrt2*sqrt(oneMez)
	sump = sqrt(p0 + p)
	deltap = sqrt(deltap)

***[ for initially massless particles
c#define Ninh_CHECK
#ifdef Ninh_CHECK
c	if( m .eq. 0 .and. i .eq. 4 ) then
	if( m .eq. 0) then
* for incoming particles
* Ref0. Eq. (2.15):
	sn1(1) = DCONJG(expIphi)
	sn1(2) = sinth

	sn2(1) = sinth
	sn2(2) = -expIphi

#ifdef QCD_GAUGE_CHECK
***[ QCD gauge check
* chose arbitrary gauge spinors spg(2) and sgm(2)
* sgm = sgp is allowed. 
* the default choice:
        sgp(1) = sn2(1)
        sgp(2) = sn2(2)

        sgm(1) = sgp(1)
        sgm(2) = sgp(2)

* another choice e.g. :
c#define GAUGE_CHOICE2
#ifdef GAUGE_CHOICE2
        sgp(1) = 1
        sgp(2) = 2

        sgm(1) = 1
        sgm(2) = 2
#endif

	ctemp1 = sgp(1)*sn1(2)-sgp(2)*sn1(1)
	ctemp2 = sgm(1)*sn1(2)-sgm(2)*sn1(1)

	print*,'particle:', i
        print*,'Reference spinors'
        print*,'sgp(1) =',sgp(1)
        print*,'sgp(2) =',sgp(2)
        print*,'sgm(1) =',sgm(1)
        print*,'sgm(2) =',sgm(2)
        print*,'ctemp1 =',ctemp1
        print*,'ctemp2 =',ctemp2

        if (abs(ctemp1) .eq. 0) then
         print*,'Vecset.F: gauge ref. vector + ZERO'
         print*,'<sgp,sn1> =', ctemp1
         print*,'Choose another sgp(i)'
         stop   
        endif    

        if (abs(ctemp2) .eq. 0) then
         print*,'Vecset.F: gauge ref. vector - ZERO'
         print*,'<sgm,sn1> =', ctemp2
         print*,'Choose another sgm(i)'
         stop
        endif
        do j1 = 1, 2
         sgp(j1) = sgp(j1)/ctemp1
         sgm(j1) = sgm(j1)/ctemp2
        enddo
***] QCD gauge check
#endif

	do j1 = 1, 2
	do j2 = 1, 2
	print*,'Components of Eps:', j1, j2
	print*,'Gauge FC R, L:', EpsTR(j1,j2), EpsTL(j1,j2)

* Ref0. Eq. (3.16): convention Eps(A, B) = Eps(Adot, B)
	if( i .le. 2 ) then ! incoming particle
	EpsTR(j1,j2) = sqrt2*DCONJG(sn2(j1))*sn1(j2)
	EpsTL(j1,j2) = sqrt2*DCONJG(sn1(j1))*sn2(j2)
        else ! out-going particle
	EpsTR(j1,j2) = sqrt2*sn1(j1)*DCONJG(sn2(j2))
	EpsTL(j1,j2) = sqrt2*sn2(j1)*DCONJG(sn1(j2))
ccc
c	EpsTL(j1,j2) = sqrt2*DCONJG(sn1(j1))*sn2(j2)
c	EpsTR(j1,j2) = sqrt2*DCONJG(sn2(j1))*sn1(j2)
        endif

#ifdef QCD_GAUGE_CHECK
***[ QCD gauge check
c	print*,'Components of Eps:', j1, j2
	print*,'Gauge 0 R, L:', EpsTR(j1,j2), EpsTL(j1,j2)
* Ref0. Eq. (3.20)
	if( i .le. 2 ) then ! incoming particle
	EpsTR(j1,j2) = sqrt2*DCONJG(sgp(j1))*sn1(j2)
	EpsTL(j1,j2) = sqrt2*DCONJG(sn1(j1))*sgm(j2)
        else ! out-going particle
	EpsTR(j1,j2) = sqrt2*sn1(j1)*DCONJG(sgp(j2))
	EpsTL(j1,j2) = sqrt2*sgm(j1)*DCONJG(sn1(j2))
c        EpsTR(j1,j2) = DCONJG(EpsTL(j1,j2))
c        EpsTL(j1,j2) = DCONJG(EpsTR(j1,j2))
c        EpsTR(j1,j2) = EpsTL(j1,j2)
c        EpsTL(j1,j2) = EpsTR(j1,j2)
        endif
***] QCD gauge check
	print*,'Gauge 1 R, L:', EpsTR(j1,j2), EpsTL(j1,j2)
#endif

	enddo
	enddo

	endif
#endif
***] for initially massless particles

* vec_EC polarization vectors (helicity flip):
      vec_EC(:,:,1,i) = vec_E(:,:,-1,i)
      vec_EC(:,:,-1,i) = vec_E(:,:,1,i)
      vec_EC(:,:,0,i) = vec_E(:,:,0,i)
      
* calculate spinors:
       spin(:,0,:,i) = 0d0
       dspin(:,0,:,i) = 0d0

	Spinor(-1, 6, 1) = deltap*sinth
	DottedSpinor(-1, 6, 1) = Spinor(-1, 6, 1)
	Spinor(-1, 6, 2) = -deltap*expIphi
	DottedSpinor(-1, 6, 2) = DCONJG(Spinor(-1, 6, 2))

	Spinor(-1, 7, 1) = sump*sinth
	DottedSpinor(-1, 7, 1) = Spinor(-1, 7, 1)
	Spinor(-1, 7, 2) = -sump*expIphi
	DottedSpinor(-1, 7, 2) = DCONJG(Spinor(-1, 7, 2))

	DottedSpinor(1, 6, 1) = sump*expIphi
	Spinor(1, 6, 1) = DCONJG(DottedSpinor(1, 6, 1))
	DottedSpinor(1, 6, 2) = sump*sinth
	Spinor(1, 6, 2) = DottedSpinor(1, 6, 2)

	DottedSpinor(1, 7, 1) = deltap*expIphi
	Spinor(1, 7, 1) = DCONJG(DottedSpinor(1, 7, 1))
	DottedSpinor(1, 7, 2) = deltap*sinth
	Spinor(1, 7, 2) = DottedSpinor(1, 7, 2)
	end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
* VecSet.F
* construct the CMS representation of a momentum and the
* three corresponding polarization vectors
* this file is part of FormCalc-9.0
* last modified 25 Mar 19 th

* i: the index of the momentum
* m, p: mass and three-momentum of the particle
* ex,ey,ez: the unit three-vector of the momentum

	subroutine VecSet_FC9(i, m, p, ex,ey,ez)
! output needed: vec0() = [spi(), EpsTR9(), EpsTL9(), EpsL9(), K9()]	
	implicit none
	integer i
	double precision m, p, ex,ey,ez

!#include "util.h"
	integer nvec0, nmom_max, minvec
! 12 = [K,Eps_i,Spi_j] with i = 1,3; j = 1, 8; 1+3+8 = 12
	parameter (nvec0 = 12)
	parameter (nmom_max = 20) ! max number of particles
	parameter (minvec = -3) ! FC9 convention

	double complex vec0(2,2,minvec:nvec0*nmom_max)
	common /vec0_FC9/ vec0

#define K9(x,y)     vec0(x,y, 1+nvec0*(i-1))
#define EpsTR9(x,y) vec0(x,y, 2+nvec0*(i-1))
#define EpsL9(x,y)  vec0(x,y, 3+nvec0*(i-1))
#define EpsTL9(x,y) vec0(x,y, 4+nvec0*(i-1))

#define Spi9(af,hel,om,d,x) vec0(x, om-5, af*2+(hel)+d+7+nvec0*(i-1))

* i.e., spi(x, 6, 5, i) = left-handed,  undotted, upper
*       spi(x, 7, 5, i) = left-handed,  undotted, lower
*       spi(x, 6, 6, i) = left-handed,  dotted,   upper
*       spi(x, 7, 6, i) = left-handed,  dotted,   lower
*       spi(x, 6, 7, i) = right-handed, undotted, upper
*       spi(x, 7, 7, i) = right-handed, undotted, lower
*       spi(x, 6, 8, i) = right-handed, dotted,   upper
*       spi(x, 7, 8, i) = right-handed, dotted,   lower

	double precision m2, p2, p0, deltap, sump
	double precision sinth, onePez, oneMez, r
	double complex expIphi, c

	double precision sqrt2
	parameter (sqrt2 = 1.4142135623730950488016887242096981D0)

	m2 = m**2
	if( m2 < 1D-14 ) then
	  p0 = p
	  deltap = 0
	else
	  p2 = p**2
	  p0 = sqrt(p2 + m2)
	  deltap = m2/(p0 + p)
	endif

	sinth = ex**2 + ey**2
	onePez = 1 + ez
	if( onePez < .5D0 ) onePez = sinth/(1 - ez)
	oneMez = 1 - ez
	if( oneMez < .5D0 ) oneMez = sinth/(1 + ez)

	if( sinth < 1D-14 ) then
* phi is irrelevant when theta = 0
	  expIphi = 1
	else
	  sinth = sqrt(sinth)
	  expIphi = DCMPLX(ex, ey)/sinth
	endif

	K9(1,1) = p0*onePez - deltap*ez
	K9(2,2) = p0*oneMez + deltap*ez
	K9(2,1) = p*DCMPLX(ex, ey)
	K9(1,2) = DCONJG(K9(2,1))

	if( m .ne. 0 ) then
	  EpsL9(1,1) = (p*onePez + deltap*ez)/m
	  EpsL9(2,2) = (p*oneMez - deltap*ez)/m
	  EpsL9(2,1) = p0/m*DCMPLX(ex, ey)
	  EpsL9(1,2) = DCONJG(EpsL9(2,1))
	endif

	EpsTR9(2,2) = 1/sqrt2*sinth
	EpsTL9(2,2) = EpsTR9(2,2)
	EpsTR9(1,1) = -EpsTR9(2,2)
	EpsTL9(1,1) = EpsTR9(1,1)
	EpsTR9(2,1) = -1/sqrt2*oneMez*expIphi
	EpsTL9(1,2) = DCONJG(EpsTR9(2,1))
	EpsTL9(2,1) = 1/sqrt2*onePez*expIphi
	EpsTR9(1,2) = DCONJG(EpsTL9(2,1))

* this is E^(I phi) cos(th/2) = 1/sqrt2*sqrt(1 + ez)*expIphi:
	expIphi = 1/sqrt2*sqrt(onePez)*expIphi
* this is sin(th/2):
	sinth = 1/sqrt2*sqrt(oneMez)

* 'big' components go as
	sump = sqrt(p0 + p)

	r = sump*sinth
	Spi9(-1,+1,7,1, 1) = r
	Spi9(-1,+1,7,2, 1) = r
	Spi9(+1,-1,7,1, 1) = r
	Spi9(+1,-1,7,2, 1) = r
	Spi9(+1,+1,6,1, 2) = r
	Spi9(+1,+1,6,2, 2) = r
	Spi9(-1,-1,6,1, 2) = -r
	Spi9(-1,-1,6,2, 2) = -r

	c = -sump*expIphi
	Spi9(-1,+1,7,1, 2) = c
	Spi9(+1,-1,7,1, 2) = c
	Spi9(-1,-1,6,2, 1) = c
	Spi9(-1,-1,6,1, 1) = DCONJG(c)
	Spi9(-1,+1,7,2, 2) = DCONJG(c)
	Spi9(+1,-1,7,2, 2) = DCONJG(c)
	Spi9(+1,+1,6,2, 1) = -c
	Spi9(+1,+1,6,1, 1) = DCONJG(-c)

* 'small' components go as
	deltap = sqrt(deltap)

	r = deltap*sinth
	Spi9(-1,+1,6,1, 1) = r
	Spi9(-1,+1,6,2, 1) = r
	Spi9(+1,-1,6,1, 1) = r
	Spi9(+1,-1,6,2, 1) = r
	Spi9(+1,+1,7,1, 2) = r
	Spi9(+1,+1,7,2, 2) = r
	Spi9(-1,-1,7,1, 2) = -r
	Spi9(-1,-1,7,2, 2) = -r

	c = -deltap*expIphi
	Spi9(-1,+1,6,1, 2) = c
	Spi9(+1,-1,6,1, 2) = c
	Spi9(-1,-1,7,2, 1) = c
	Spi9(-1,-1,7,1, 1) = DCONJG(c)
	Spi9(-1,+1,6,2, 2) = DCONJG(c)
	Spi9(+1,-1,6,2, 2) = DCONJG(c)
	Spi9(+1,+1,7,2, 1) = -c
	Spi9(+1,+1,7,1, 1) = DCONJG(-c)
	end

