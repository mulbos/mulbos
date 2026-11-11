* WARNING: process dependence (MQx2 and quark color factor appear), 
* becareful with fac_color (should be improved); 
* This works for EW and trivial QCD (one colour flow, with 
* only quark-gluon coupling). There is no colour correlation in general. 

* Ref0. Dittmaier hep-ph/9904440v2
* Ref1. Catani and Seymour hep-ph/9605323
* Ref2. Dittmaier and Huber, arXiv:0911.2329v2.
* Author: Le Duc Ninh (2011)

#define PDF_COUNTER

***********************************************************************
	double precision function FuncDPsing_coll(A0sq, fac, 
     &	 facx1, facx2, facx3, ifail)
        use global_vars
	implicit none

#include "model.h"

	double precision Cij
	external Cij

	integer i, j, ifail, nlegs_qqVV_max, nlegs
        parameter (nlegs_qqVV_max = 5) ! for qq > VVV maximally

	double precision charge(nlegs_qqVV_max), r
	
	double precision Jacb, fac, A0sq(10)
	double precision facx1, facx2, facx3
	common /Cijcom/ Jacb

      charge(:) = 0d0
      if (lNLO_EW) then
        charge(1) = charge_qed(1)
        charge(2) = charge_qed(2)
        charge(3) = charge_qed(3) + charge_qed(4)
        charge(4) = charge_qed(5) + charge_qed(6)
       if (lTriBoson) then
        charge(5) = charge_qed(7) + charge_qed(8)
       endif
      else
       charge(1) = charge_qcd(1)
       charge(2) = charge_qcd(2)
!       charge(3) = 0D0
!       charge(4) = 0D0
      endif

      if (lDiBoson) then
       nlegs = 4 ! qq > VV case
      else
       nlegs = 5 ! qq > VVV case
      endif
	
	FuncDPsing_coll = 0

	do i = 1, nlegs
	  do j = 1, nlegs

      if (lNLO_COLL_RAD) then
	    r = charge(i)*charge(j)
      endif

      if (lNLO_COLL_IND) then
C add a wrong minus sign to correct below
	    r = -charge(i)**2
      endif

	    if( r .ne. 0 ) then
	    Jacb = sign(1, 2 - i)*sign(1, 2 - j)*r  

      if (lNLO_QCD) then
* convention: add a factor CF = 4/3 compared to QED 
            Jacb = -2*Alfas/(3*pi)*Jacb*fac
      else if (lNLO_EW) then
            Jacb = -Alfa/(2*pi)*Jacb*fac
      endif

	    FuncDPsing_coll = FuncDPsing_coll + 
     &	    Jacb*Cij(i, j, A0sq, facx1, facx2, facx3, ifail)

c            print*,'Func =', i, j, FuncDPsing_coll ! Ninh

	    endif
	  enddo
	enddo

c        stop ! Ninh
	end


***********************************************************************
	double precision function 
     &  Cij(i, j, A0sq, facx1, facx2, facx3, ifail)
        use global_vars
        use global_mom_xpdf
	implicit none
	double precision facx1, facx2, facx3	
	integer i, j, ifail, iplus
        integer cfail

#include "model.h"

	double precision mi, mj, Pij2
	double precision mi2, mj2, MQx, MQx2

	double precision x, y, tau, z1, pp, kk
      double precision log_x, log_1mx, log_1mxOx, rat_1px2O1mx
      double precision pol_x2, pol_x1, x2

	double precision Jacb, A0sq(10)
	common /Cijcom/ Jacb
        common /plusdist/ iplus

      double precision qext(proc_nlegs,4)
      double precision delta, Kqqq, CF
      double precision scale_f2, ecm_sq

        Cij = 0 
        cfail = 1 ! initialize cut switch
	iplus = ifail

	if( i .eq. j ) return

      if (lNLO_QCD) then
       CF = 4D0/3
      else
       CF = 1D0
      endif
      
      scale_f2 = scale**2 ! factorization scale
      ecm_sq = ecm**2 ! square of ecm

* convention: i = emitter, j = spectator ***
C Important: we sum over q_bar q, q_bar g and q g here. 
C And we fix the Born amplitude is q_bar q. 
C So, for the gluon induced process, if the emitter is 1 
C then the NLO process is g q (hence pdfs = g*q = facx3); 
C if the emitter is 2 then the NLO process is q_bar g (facx2).

*** common kinematics ***
        if ( i.gt.2 ) then
	mi = mass_VV(i-2)
        else
        mi = MU
        endif
        if ( j.gt.2 ) then
	mj = mass_VV(j-2)
        else
        mj = MU
        endif

***[for MasReg
      if (lMassReg) then
* WARNING: this is for light quarks (u,d,c,s,b) in the initial state, 
* MQx is regarded as a regulator here.

	MQx = 1D-3
	MQx2 = MQx*MQx
      endif
***]for MasReg

	mi2 = mi*mi
	mj2 = mj*mj
	
	x = y_coll
	tau = tau_pdf**2

c        print*,'FuncDP x=',x 
	if( x .lt. 1D-10 .or. x .gt. (1 - 1D-10) ) return
      log_x = log(x)
      x2 = x*x
      pol_x1 = 1-x
      log_1mx = log(pol_x1)
      log_1mxOx = log(pol_x1/x)
      rat_1px2O1mx = (1 + x2)/pol_x1
      pol_x2 = pol_x1**2 + x2

***[initial - initial
	if (i .le. 2 .and. j .le. 2) then 

      if (lDimReg) then
***[for DimReg

***[for photon radiation
      if (lNLO_COLL_RAD) then

      if (leg(1).eq.leg(2)) then
* Not sure???
	print*,'Error in FuncDPcoll_had.F: 
     >  special case not done, code KernAA'
      endif

***[plus distribution
        pp = -rat_1px2O1mx
     >   *log(tau*ecm_sq/scale_f2)
        kk = -4*log_1mx/pol_x1 + 2*log_x/pol_x1
***]plus distribution

***[nonplus function
        if (ifail .eq. 1) then
         pp = pp - rat_1px2O1mx*log_x
         kk = kk + 2*(1 + x)*log_1mx - 1 + x 
     &   - log_x*(1+x)
	endif
***]nonplus function

***[delta distribution
        if (ifail .eq. 0) then
         z1 = pi**2/3D0 + 5 - pi**2
         kk = kk - z1
	endif
***]delta distribution

      if (lPDF_SCHEME_DIS) then
***[plus distribution
* for the DIS factorization scheme
	z1 = rat_1px2O1mx*(log_1mxOx - 3/4D0) + (9 + 5*x)/4
	kk = kk + z1
***]plus distribution
      endif

	Cij = pp + kk
* include PDFs
	Cij = -Cij*facx1*A0sq(1)

      endif
***]for photon radiation

***[for photon induced
* only nonplus functions appear here
      if (lNLO_COLL_IND) then
        if (ifail .eq. 1) then
	pp = pol_x2*(
     >  - log(x*tau*ecm_sq/scale_f2))
	kk = pol_x2*(
     >  - log_1mxOx - log_1mx) - 2*x*pol_x1 

      if (lPDF_SCHEME_DIS) then
* for the DIS factorization scheme
	kk = kk + pol_x2*log_1mxOx - 8*x2 + 8*x - 1
      endif

	z1 = pp + kk

      if (lNLO_QCD) then
* color factor for gluon induce (see arXiv:0806.1433)
        z1 = z1/2D0
* convert to the color factor CF=4/3 of radiation process
        z1 = z1*3/4D0
      endif

      if (lNLO_EW) then
* color factor for photon induce
        z1 = z1*3
      endif	

* include PDFs
	if ( i .gt. 1 ) then
	 z1 = -z1*facx2*A0sq(1)
	else
	 z1 = -z1*facx3*A0sq(1)
	endif

	Cij = Cij + z1

***[ the gam-gam contribution
      if (lWW.and.lNLO_EW) then
* Ref1. Eq. (10.25), (C.6)
	pp = (2/x + x - 2)*(
     >  - log(x*tau*ecm_sq/scale_f2))
* Ref1. Eq. (10.24): a = q, a' = photon 
* Ref1. Eq. (C.19):
	kk = -(2/x + x - 2)*log_1mx
* Ref1. Eq. (C.15):
	kk = kk - (2/x + x - 2)*log_1mxOx - x

	z1 = pp + kk
* undo the color factor
	z1 = z1/(fac_color)
* include PDFs
	if ( i .gt. 1 ) then
	 z1 = -z1*facx3*A0sq(2)
	else
	 z1 = -z1*facx2*A0sq(2)
	endif

	Cij = Cij + z1
      endif
***] the gam-gam contribution

	endif

***[ the gam-gam contribution
      if (lWW.and.lNLO_EW) then
      if (lPDF_SCHEME_DIS) then
***[plus distribution
* for the DIS factorization scheme
* Ref1. Eq. (C.25):
	z1 = rat_1px2O1mx*(log_1mxOx - 3/4D0) + (9 + 5*x)/4
* undo the color factor
	z1 = z1/(fac_color)
* include PDFs
	if ( i .gt. 1 ) then
	 z1 = z1*facx3*A0sq(2)
	else
	 z1 = z1*facx2*A0sq(2)
	endif
	Cij = Cij + z1
***]plus distribution
      endif
      endif
***] the gam-gam contribution

      endif
***]for photon induced

***]for DimReg
      else
***[for MasReg	
     
***[for photon radiation
      if (lNLO_COLL_RAD) then
        Cij = rat_1px2O1mx
     >   *log(tau*ecm_sq/MQx2) 
     >   - 2*x*(1D0+x)/(1D0-x2)

* include PDFs
	Cij = Cij*facx1*A0sq(1)

#ifdef PDF_COUNTER
	z1 = rat_1px2O1mx*log(scale_f2/MQx2)
     &   - rat_1px2O1mx*(2*log_1mx+1) 
      if (lPDF_SCHEME_DIS) then
* for the DIS factorization scheme
	z1 = z1 + rat_1px2O1mx*(log_1mxOx - 3/4D0) + (9 + 5*x)/4
      endif

* include PDFs
	z1 = z1*facx1*A0sq(1)

	Cij = Cij - z1
c	print*,"Cij in in=",i,j,Cij
#endif
      endif
***]for photon radiation

***[for photon induced
      if (lNLO_COLL_IND) then
        if (ifail .eq. 1) then
	z1 = pol_x2
     >   *log(tau*ecm_sq*pol_x1**2/MQx2) + 2*pol_x1*x

      if (lNLO_QCD) then
* color factor for gluon induce (see arXiv:0806.1433)
        z1 = z1/2D0
* convert to the color factor CF=4/3 of radiation process
        z1 = z1*3/4D0
      endif

      if (lNLO_EW) then
* color factor for photon induce
        z1 = z1*3
      endif	

* include PDFs
	if ( i .gt. 1 ) then
	 z1 = z1*facx2*A0sq(1)
	else
	 z1 = z1*facx3*A0sq(1)
	endif

	Cij = Cij + z1	 	      

#ifdef PDF_COUNTER
	z1 = log(scale_f2/MQx2)*pol_x2 

      if (lNLO_QCD) then
* color factor for gluon induce (see arXiv:0806.1433)
        z1 = z1/2D0
* convert to the color factor CF=4/3 of radiation process
        z1 = z1*3/4D0
      endif

      if (lPDF_SCHEME_DIS) then
* for the DIS factorization scheme
	z1 = z1 + pol_x2*log_1mxOx - 8*x2 + 8*x - 1
      endif

      if (lNLO_EW) then
* color factor for photon induce
        z1 = z1*3
      endif

* include PDFs
	if ( i .gt. 1 ) then
	 z1 = z1*facx2*A0sq(1)
	else
	 z1 = z1*facx3*A0sq(1)
	endif
	Cij = Cij - z1
#endif

***[ the gam-gam contribution
      if (lWW.and.lNLO_EW) then
	z1 = (2 + x2 - 2*x)/x
     >  *log(tau*ecm_sq*pol_x1**2/MQx2/x2) - 2*pol_x1/x

* undo the color factor
	z1 = z1/(fac_color)
* include PDFs
* watch the order: i = 1 means the first photon 
* is emitted from the first parton A, i = 2 means the second 
* photon is emitted from the second parton B. 
* The two parent processes are: A+photon -> X, photon+B -> X. 
	if ( i .gt. 1 ) then
	 z1 = z1*facx3*A0sq(2)
	else
	 z1 = z1*facx2*A0sq(2)
	endif
	Cij = Cij + z1

#ifdef PDF_COUNTER
* Ref2. Eq. (3.20)
	z1 = (log(scale_f2/MQx2) - 2*log_x -1)*(2 + x2 - 2*x)/x
* undo the color factor
	z1 = z1/(fac_color)
* include PDFs
	if ( i .gt. 1 ) then
	 z1 = z1*facx3*A0sq(2)
	else
	 z1 = z1*facx2*A0sq(2)
	endif
	Cij = Cij - z1
#endif

      endif
***] the gam-gam contribution

	endif

***[ the gam-gam contribution
      if (lWW.and.lNLO_EW) then
#ifdef PDF_COUNTER
      if (lPDF_SCHEME_DIS) then
* for the DIS factorization scheme
* This is a plus distribution
* Ref2. Eq. (3.25)
	z1 = - ( 
     &  rat_1px2O1mx*(log_1mxOx - 3/4D0) + (9 + 5*x)/4
     &   )
* undo the color factor
	z1 = z1/(fac_color)
* include PDFs
	if ( i .gt. 1 ) then
	 z1 = z1*facx3*A0sq(2)
	else
	 z1 = z1*facx2*A0sq(2)
	endif
	Cij = Cij - z1
      endif
#endif
      endif
***] the gam-gam contribution

      endif
***]for photon induced

***]for MasReg	
      endif

c        print*,'Cij before cut =',i,j,Cij   

	call GetMomenta_BORN(qext)
c	print*,"qext=",qext

***[impose cuts
C physical cuts
	call DoCut_parCM_coll(cfail, qext(1,1), i, ifail)

	if ( cfail .eq. 0 ) then
         Cij = 0
c         print*,'Cij after cut =',i,j,Cij   
         return
        endif
***]impose cuts

        if (lHist_Fill) then
C the momenta here has to be the OS-projected kinematics,
c to match with the subtraction term
        call FillHIS_coll(Cij, Jacb, qext, i, ifail)
        endif

	return
	endif 
***]initial - initial

      if (lNLO_COLL_RAD) then
***[final (i) - initial (j)
	if (i .ge. 3 .and. j .le. 2) then 

***[for MasReg
      if (lMassReg) then
! Eq. (A.2) of Dittmaier 1999
c	Pij2 = -2*DBLE(Pair(k(i), k(j)))+mi2
      Pij2 = oext(i,1)*oext(j,1)-oext(i,2)*oext(j,2)-oext(i,3)*oext(j,3)-oext(i,4)*oext(j,4)
      Pij2 = -2*Pij2+mi2
	z1 = mi2*x/(mi2-Pij2*pol_x1)
	
        Cij = 2*log((2-x-z1)/pol_x1)
     >	+0.5D0*(z1-1)*(3+z1-4*mi2*x/(Pij2-mi2)/pol_x1)
        Cij = Cij/pol_x1 
* include PDFs
	Cij = Cij*facx1*A0sq(1)
      endif
***]for MasReg
c        print*,'Cij =',i,j,Cij   

***[for DimReg
* Ninh checking
      if (lDimReg) then
* Ref2. Eq. (6.67) and Eq. (6.53)

	if ( ifail .eq. 1 ) then
c	Pij2 = 2*DBLE(Pair(k(i), k(j)))/x
        Pij2 = oext(i,1)*oext(j,1)-oext(i,2)*oext(j,2)-oext(i,3)*oext(j,3)-oext(i,4)*oext(j,4)
        Pij2 = 2*Pij2/x
	else
c	Pij2 = 2*DBLE(Pair(k(i), k(j)))
        Pij2 = oext(i,1)*oext(j,1)-oext(i,2)*oext(j,2)-oext(i,3)*oext(j,3)-oext(i,4)*oext(j,4)
        Pij2 = 2*Pij2
	endif

* Ref2. Eq. (6.55)
        kk = -Kqqq(x, Pij2, mi) 

* include PDFs
	Cij = -kk*facx1*A0sq(1)

      endif
***]for DimReg

	call GetMomenta_BORN(qext)

***[impose cuts
C physical cuts
	call DoCut_parCM_coll(cfail, qext(1,1), j, ifail)

	if ( cfail .eq. 0 ) then
         Cij = 0
         return
        endif
***]impose cuts

        if (lHist_Fill) then
C the momenta here has to be the OS-projected kinematics,
c to match with the subtraction term
        call FillHIS_coll(Cij, Jacb, qext, j, ifail) 
        endif

        return
        endif 
***]final (i) - initial (j)

***[initial (i) - final (j)
	if (i .le. 2 .and. j .ge. 3) then 

***[for MasReg
      if (lMassReg) then
c	Pij2 = -2*DBLE(Pair(k(i), k(j)))+mj2
        Pij2 = oext(i,1)*oext(j,1)-oext(i,2)*oext(j,2)-oext(i,3)*oext(j,3)-oext(i,4)*oext(j,4)
        Pij2 = -2*Pij2+mj2
	z1 = mj2*x/(mj2-Pij2*pol_x1)		
	
        Cij = (1+x2)/pol_x1*(log((mj2-Pij2)/(MQx2*x))+log(1-z1)-1) 
     >	-2/pol_x1*log(2-x-z1) + (1+x)*log_1mx+pol_x1
* include PDFs
	Cij = Cij*facx1*A0sq(1)

#ifdef PDF_COUNTER
	z1 = rat_1px2O1mx*log(scale_f2/MQx2)
     &   - rat_1px2O1mx*(2*log_1mx+1) 
      if (lPDF_SCHEME_DIS) then
* for the DIS factorization scheme
	z1 = z1 + rat_1px2O1mx*(log_1mxOx - 3/4D0) + (9 + 5*x)/4
      endif

* include PDFs
	z1 = z1*facx1*A0sq(1)

	Cij = Cij - z1
#endif
      endif
***]for MasReg
c        print*,'Cij =',i,j,Cij, z1, Pij2   

***[for DimReg
* Ninh checking
      if (lDimReg) then

* Ref2. Eq. (6.67) and Eq. (6.53)

	if ( ifail .eq. 1 ) then
c	Pij2 = 2*DBLE(Pair(k(i), k(j)))/x
        Pij2 = oext(i,1)*oext(j,1)-oext(i,2)*oext(j,2)-oext(i,3)*oext(j,3)-oext(i,4)*oext(j,4)
        Pij2 = 2*Pij2/x
	y = x
	else
c	Pij2 = 2*DBLE(Pair(k(i), k(j)))
        Pij2 = oext(i,1)*oext(j,1)-oext(i,2)*oext(j,2)-oext(i,3)*oext(j,3)-oext(i,4)*oext(j,4)
        Pij2 = 2*Pij2
	y = 1
	endif

***[plus*nonplus distribution
        pp = -rat_1px2O1mx
     >   *log(Pij2*y/scale_f2)
***]plus*nonplus distribution

* Ref2. Eq. (6.55)
       kk = 
     & 1/CF*(ifail*((-CF)*(1+x)*log(pol_x1*Pij2/(pol_x1*Pij2+mj**2))) + 
     & CF*3/2*delta(1-x)*(log((Pij2-2*mj*sqrt(Pij2+mj**2)+2*mj**2)/Pij2) + 
     & 2*mj/(sqrt(Pij2+mj**2)+mj))) 
       kk = -kk

***[PDF counter (K_bar and K_FS: Ref1. Eq. (10.24))
***[plus distribution
       z1 = -2/pol_x1*log_1mxOx
***]plus distribution

***[nonplus function
       if ( ifail .eq. 1 ) then
       z1 = z1 + (1+x)*log_1mxOx - 1 + x
       endif 
***]nonplus function

***[delta distribution
        if (ifail .eq. 0) then
         z1 = z1 - 5 + pi**2
	endif
***]delta distribution

	kk = kk + z1

      if (lPDF_SCHEME_DIS) then
***[plus distribution
* for the DIS factorization scheme
	z1 = rat_1px2O1mx*(log_1mxOx - 3/4D0) + (9 + 5*x)/4
	kk = kk + z1
***]plus distribution
      endif
***]PDF counter (K_bar and K_FS: Ref1. Eq. (10.24))

	Cij = pp + kk
* include PDFs
	Cij = -Cij*facx1*A0sq(1)
      endif
***]for DimReg

	call GetMomenta_BORN(qext)

***[impose cuts
C physical cuts
	call DoCut_parCM_coll(cfail, qext(1,1), i, ifail)

	if ( cfail .eq. 0 ) then
         Cij = 0
         return
        endif
***]impose cuts

        if (lHist_Fill) then
C the momenta here has to be the OS-projected kinematics,
c to match with the subtraction term
        call FillHIS_coll(Cij, Jacb, qext, i, ifail) 
        endif

        return 
        endif 		
***]initial (i) - final (j)

***[final - final
	if (i .ge. 3 .and. j .ge. 3) then

        Cij = 0 

        return
	endif
***]final - final
      endif

	end
***********************************************************************
* for spectator = quark
* Spliting: qq
      double precision function Kqqq(x, s, m)
      implicit none
      double precision x, s, m, m2, temp, plusd, delta
      double precision mu, JagQ
      integer iplus 
      common /plusdist/ iplus
      Kqqq = 0
      m2 = m**2
      mu = m/sqrt(s)
* Ref2. Eq. (6.58)
      Kqqq = 2*(plusd(log(1-x)/(1-x)) - log(2-x)/(1-x)*iplus) 
* for terms of the type: plusd(f[x])*g[x]
      if ( iplus .eq. 1 ) then
      temp = 2*plusd(1D0/(1-x))*log((2-x)*s/((2-x)*s+m2))
      else
      temp = 2*plusd(1D0/(1-x))*log((2-1)*s/((2-1)*s+m2))
      endif
      Kqqq = Kqqq + temp - 3D0/2*delta(1-x) 
     & + (m2/s*log(m2/(s+m2)) + 1D0/2*m2/(s+m2))*delta(1-x) 
     & + JagQ(x, mu)
      return
      end
***********************************************************************
      double precision function JagQ(x, mu)
      implicit none
      double precision x, mu, temp, plusd, delta
      integer iplus 
      common /plusdist/ iplus
      JagQ = 0
* Ref2. Eq. (5.58)
* the funtion J^a_{gQ}(x, mu)

* for terms of the type: plusd(f[x])*g[x]
      if ( iplus .eq. 1 ) then
      temp = plusd(2/(1-x))*log(2+mu**2-x)
      else
      temp = plusd(2/(1-x))*log(2+mu**2-1)
      endif
      temp = temp
     & + plusd((1-x)/(2*(1-x+mu**2)**2) - 2/(1-x)*(1+log(1-x+mu**2))) 
      JagQ = temp
      return
      end
***********************************************************************
      double precision function plusd(xx)
      implicit none
      double precision xx
      plusd=xx
      return
      end

      double precision function delta(yy)
      implicit none
      double precision yy
      integer iplus
      common /plusdist/ iplus
      if ( iplus .eq. 0 ) then
* Adding a minus here because we minus the delta(1-x) 
* contribution in xsection.F (minus minus = plus). 
      delta = -1d0
      else
      delta = 0d0
      endif
      return
      end
***********************************************************************
