* Author: Le Duc Ninh
* To fill in distributions for all cases
* last modifiaction Feb 2010

* For Tree, Virt, Dipole(Slicing)-Finite cases

	subroutine FillHIS(Gij, Jacb, fext)
	use global_vars
	use global_mom_xpdf
C Note: fext are in partonic center-of-mass system frame
	implicit none
	integer i, j

c#include "xsection.h"	

        double precision qext(proc_nlegs,4), kext, rgam, rbetgam
	double precision FUNC, Gij, Jacb, x1, x2
	double precision rext(proc_nlegs,4), fext(proc_nlegs,4)

* factor 1/2 since we will fill twice
        FUNC = Gij*Jacb/2

	do j = 1,4 
	 do i = 1,proc_nlegs
	  rext(i,j) = fext(i,j)
	  qext(i,j) = fext(i,j)
	 enddo
	enddo

***{boost parameters
* mq = mP = 0 here
* boost for the subprocess: P1/PARTON1 + P2/PARTON2 -> XXX
* Rules: OZ(x1), -OZ(x2)
	x2 = tau_pdf**2
	x2 = x2 + (1-x2)*x_pdf

        x1 = tau_pdf**2/x2

        rgam = 0.5D0*(x1 + x2)/sqrt(x1*x2)        

        rbetgam=0.5D0*(x2 - x1)/sqrt(x1*x2)
	do i = 1, proc_nlegs
	call boost0(rgam, rbetgam, 
     &  rext(i,1), rext(i,2), rext(i,3), rext(i,4))
	enddo
***}boost parameters
	call Fill(rext, FUNC)

***{boost parameters
* mq = mP = 0 here
* boost for the subprocess: P1/PARTON2 + P2/PARTON1 -> XXX
* Rules: OZ(x1), -OZ(x2)
* just exchange x1 <-> x2:

        rbetgam = 0.5D0*(x1 - x2)/sqrt(x1*x2)

	do i = 3, proc_nlegs
* rotation symmetry
	qext(i,2) = qext(i,2)
	qext(i,3) = -qext(i,3)
	qext(i,4) = -qext(i,4)
	enddo

	do i = 1, proc_nlegs
	call boost0(rgam, rbetgam, 
     &  qext(i,1), qext(i,2), qext(i,3), qext(i,4))
	enddo
***}boost parameters
	call Fill(qext, FUNC)

	end

* For NLO collinear (convolution) part.
* One has to boost since each initial particle can 
* radiate photon (gluon).

	subroutine FillHIS_coll(Gij, Jacb, fext, ie, ifail)
	use global_vars
	use global_mom_xpdf
C Note: fext are in partonic center-of-mass system frame
	implicit none
	integer i, j, ie, ifail

c#include "xsection.h"	

        double precision qext(proc_nlegs,4), kext
	double precision FUNC, Gij, Jacb
	double precision rext(proc_nlegs,4), fext(proc_nlegs,4)
        double precision rgam,rbetgam,x1,x2

* with factor 1/2 since we fill twice
        FUNC = Gij*Jacb/2D0

***{boost parameters 1
* ME = 0 here
* to be in the CMS of tau_pdf*ecm

        if (ie .gt. 2) then
	 print*,'FillHIS.F error: ie > 2'
	 print*,'ie =', ie
	 stop
	endif 

* Either
	if ( ifail .eq. 1 ) then
        x1 = y_coll
	else
        x1 = 1D0
	endif
* Or
* This sets x1 = 1 (Boost off) automatically 
* for the [x1=1]-term in the plus distribution. 
c        x1 = (sqrtS/(tau_pdf*ecm))**2

        rgam = 0.5D0*(x1+1)/sqrt(x1)        

        if (ie .eq. 1) then 
        rbetgam = 0.5D0*(1-x1)/sqrt(x1)
        else     
        rbetgam = -0.5D0*(1-x1)/sqrt(x1)
	endif
***}boost parameters 1

	do j = 1,4 
	 do i = 1,proc_nlegs
	  qext(i,j) = fext(i,j)
	 enddo
	enddo

        do i=1,proc_nlegs
        call boost0(rgam, rbetgam, 
     &  qext(i,1),qext(i,2),qext(i,3),qext(i,4))
	enddo	

	do j = 1,4 
	 do i = 1,proc_nlegs
	  rext(i,j) = qext(i,j)
	 enddo
	enddo

***{boost parameters 2a
* mq = mP = 0 here
* boost for the subprocess: P1/PARTON1 + P2/PARTON2 -> XXX
* Rules: OZ(x1), -OZ(x2)
	x2 = tau_pdf**2
	x2 = x2 + (1-x2)*x_pdf

        x1 = tau_pdf**2/x2

        rgam = 0.5D0*(x1 + x2)/sqrt(x1*x2)        

        rbetgam=0.5D0*(x2 - x1)/sqrt(x1*x2)
	do i = 1, proc_nlegs
	call boost0(rgam, rbetgam, 
     &  rext(i,1), rext(i,2), rext(i,3), rext(i,4))
	enddo
***}boost parameters 2a
	call Fill(rext, FUNC)

***{boost parameters 2b
* mq = mP = 0 here
* boost for the subprocess: P1/PARTON2 + P2/PARTON1 -> XXX
* Rules: OZ(x1), -OZ(x2)
* just exchange x1 <-> x2

        rbetgam=0.5D0*(x1 - x2)/sqrt(x1*x2)

	do i = 3, proc_nlegs
* rotation symmetry
	qext(i,2) = qext(i,2)
	qext(i,3) = -qext(i,3)
	qext(i,4) = -qext(i,4)
	enddo

	do i = 1, proc_nlegs
	call boost0(rgam, rbetgam, 
     &  qext(i,1), qext(i,2), qext(i,3), qext(i,4))
	enddo
***}boost parameters 2b
	call Fill(qext, FUNC)

        end
******************************************************
	subroutine Fill(fext, FUNC)
C Note: fext are in LAB frame with order q1 q2 > l3 l4 ... jet
        use, intrinsic :: ieee_arithmetic
	use global_vars
	implicit none
        double precision fext(proc_nlegs,4)
        double precision FUNC

***[check for NaN and Infinity
        if (ISNAN(FUNC) .or. .not. ieee_is_finite(FUNC)) then
         return
        endif
***]check for NaN and Infinity

!        print*,"hello from Fill()"

	if (lDiBoson) then
	 call Fill_DiBoson(fext, FUNC)
	else
	 call Fill_TriBoson(fext, FUNC)
	endif
        
	end subroutine Fill
******************************************************
	subroutine Fill_TriBoson(fext, FUNC)
C Note: fext are in LAB frame with order q1 q2 > l3 l4 ... jet
	use global_vars
	use global_mom_xpdf
	implicit none
      double precision fext(proc_nlegs,4)
      double precision qext(proc_nlegs,4)
      double precision FUNC
      integer i, j

        do i = 1, 4
         do j = 1, proc_nlegs   
         qext(j,i) = fext(j,i)
         enddo
        enddo

!      print*,"hello from Fill_TriBoson()"

***[Fill Histograms
      if (lWWW) then
       call Fill_WWW(qext(1,1), FUNC)
      endif
***]Fill Histograms

        end subroutine Fill_TriBoson
******************************************************
	subroutine Fill_DiBoson(fext, FUNC)
C Note: fext are in LAB frame with order q1 q2 > l3 l4 ... jet
        use, intrinsic :: ieee_arithmetic
	use global_vars
	use global_mom_xpdf
	use global_cuts
	implicit none
	integer i, j

      double precision fext(proc_nlegs,4)
      double precision qext(proc_nlegs,4)
      double precision FUNC

      double precision phi_lep(Nlepton), phi_W, phi_Z
      double precision phot(7) ! 7 for E, px, py, pz, pt, pseudo-rap, phi
      double precision clep_bared(7,Nlep_charged), cLep_dressed(7,Nlep_charged)
      
!      print*,"hello from FillHIS, Nlep_charged =",Nlep_charged

        do i = 1, 4
         do j = 1, proc_nlegs   
         qext(j,i) = fext(j,i)
         enddo
        enddo

***[photon-chargedlepton recombination
c Final-state photons with rapidity |y| > 5 are
c considered as lost in the beam pipe and excluded from recombination. 
      IF (lNLO_FIN_RAD) THEN
      IF (lNLO_EW) THEN

      do i = 1, 4
       phot(i) = qext(proc_nlegs,i)
      enddo

      call lego( phot(1), phot(5), phot(6), phot(7) )
      If ( abs(phot(6)) .le. py_phot_recomb ) then

      do j = 1,4
       if (lWZ) then
        clep_bared(j,1) = qext(3,j) ! electron
        clep_bared(j,2) = qext(5,j) ! muon +
        clep_bared(j,3) = qext(6,j) ! muon -
       else if (lWW) then
        clep_bared(j,1) = qext(3,j) ! electron +
        clep_bared(j,2) = qext(5,j) ! muon -
       else if (lZZ) then
        clep_bared(j,1) = qext(3,j) ! electron +
        clep_bared(j,2) = qext(4,j) ! electron -
        clep_bared(j,3) = qext(5,j) ! muon +
        clep_bared(j,4) = qext(6,j) ! muon -
       endif 
      enddo

      do i = 1, Nlep_charged
      call lego(clep_bared(1,i),clep_bared(5,i),clep_bared(6,i),clep_bared(7,i))
      enddo

      call PhotLepComb(Nlep_charged,DRLA,phot(1),clep_bared(1,1),clep_dressed(1,1))

      do j = 1, 4
       if (lWZ) then
        qext(3,j) = clep_dressed(j,1)
        qext(5,j) = clep_dressed(j,2)
        qext(6,j) = clep_dressed(j,3)
       else if (lWW) then
        qext(3,j) = clep_dressed(j,1)
        qext(5,j) = clep_dressed(j,2)
       else if (lZZ) then
        qext(3,j) = clep_dressed(j,1)
        qext(4,j) = clep_dressed(j,2)
        qext(5,j) = clep_dressed(j,3)
        qext(6,j) = clep_dressed(j,4)
       endif
      enddo

      Endif
      ENDIF
      ENDIF
***]photon-chargedlepton recombination

***[Fill Histograms
      if (lWZ) then
       call Fill_WZ(qext, FUNC)
      else if (lWW) then
       call Fill_WW(qext, FUNC)
      else if (lZZ) then
       call Fill_ZZ(qext, FUNC)
      endif
***]Fill Histograms

	end subroutine Fill_DiBoson
***********************************************************************
        subroutine GetMomNeutrino(pe,pn,MW2,pn_new)
        implicit none
        double precision pe(4),pn(4),pn_new(4),MW2,pz
        double precision pe0,pex,pey,pez,pnx,pny
        double precision a,b,Delta

        double complex sqrt_Del, pz1,pz2
        integer i
C No need to use this because ATLAS unfold results to the true neutrino momentum.
cc#define USE_RECONSTRUCTED_NEU_MOM

#ifdef USE_RECONSTRUCTED_NEU_MOM
        pe0 = pe(1)
        pex = pe(2)
        pey = pe(3)
        pez = pe(4)
        pnx = pn(2)
        pny = pn(3)

        a = 8*(pe0**2 - pez**2)
        b = pez*(-4*MW2 + 4*pe0**2 - 4*pex**2 - 4*pey**2 - 4*pez**2 - 8*pex*pnx - 8*pey*pny)
        Delta = pez**2*(-4*MW2 + 4*pe0**2 - 4*pex**2 - 4*pey**2 - 4*pez**2 - 8*pex*pnx - 8*pey*pny)**2 - 
     -      4*(4*pe0**2 - 4*pez**2)*(-MW2**2 + 2*MW2*pe0**2 - pe0**4 - 2*MW2*pex**2 + 2*pe0**2*pex**2 - pex**4 - 2*MW2*pey**2 + 
     -         2*pe0**2*pey**2 - 2*pex**2*pey**2 - pey**4 - 2*MW2*pez**2 + 2*pe0**2*pez**2 - 2*pex**2*pez**2 - 2*pey**2*pez**2 - 
     -         pez**4 - 4*MW2*pex*pnx + 4*pe0**2*pex*pnx - 4*pex**3*pnx - 4*pex*pey**2*pnx - 4*pex*pez**2*pnx + 4*pe0**2*pnx**2 - 
     -         4*pex**2*pnx**2 - 4*MW2*pey*pny + 4*pe0**2*pey*pny - 4*pex**2*pey*pny - 4*pey**3*pny - 4*pey*pez**2*pny - 
     -         8*pex*pey*pnx*pny + 4*pe0**2*pny**2 - 4*pey**2*pny**2)
        sqrt_Del = sqrt(DCMPLX(Delta))

        pz1 = (-b - sqrt_Del)/a
        pz2 = (-b + sqrt_Del)/a
c        print*,"check pz1   =",pz1
c        print*,"check pz2   =",pz2
c        print*,"check MW    =",sqrt(MW2)
c        stop

        if (abs(DBLE(pz1)).le.abs(DBLE(pz2))) then
         pz = DBLE(pz1)
        else
         pz = DBLE(pz2)
        endif

        pn_new(2) = pn(2)
        pn_new(3) = pn(3)
        pn_new(4) = pz
        pn_new(1) = sqrt(pn_new(2)**2 + pn_new(3)**2 + pn_new(4)**2)
#else
        do i = 1, 4
         pn_new(i) = pn(i)
        enddo
#endif
        end
***********************************************************************
	subroutine theta_star(qlep, qW, c_theta_star, phi_star)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Refs: https://arxiv.org/abs/1103.5445 (read Fig 10)
c and also https://arxiv.org/abs/1204.6427
c this routine works for both W and Z
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        use global_vars
	use global_mom_xpdf
	implicit none
	double precision qW(4)
	integer i, j, ipart

       double precision qlep(4,2) ! only 2 decay leptons here
       double precision c_theta_star(2), phi_star(2), qW_inv(4), 
     & qlep_Wrest(4,2), dot_space
      double precision proton_lab(4,2), proton_Wrest(4,2),
     & lep_P_Wrest(2), W_P_Wrest
      double precision proton_P_Wrest(2), c_theta_proton(2),
     & qlep_Wrest_xyz(2:4,2), vectorY(2:4)
       double precision temp(2), tmp(2), vectorX(2:4)

#include "model.h"	

***[cos_theta_star calculation
      qW_inv(1) = qW(1)
      do i = 2, 4
      qW_inv(i) = -qW(i)
      enddo

      W_P_Wrest = sqrt(qW(2)**2 + qW(3)**2 + qW(4)**2)

C boost from the Lab (proton-proton CMS) frame to the W-rest frame
      do i = 1, 2
      call BOOSTN(qlep(1,i),qW_inv(1),qlep_Wrest(1,i))
      lep_P_Wrest(i) = sqrt(qlep_Wrest(2,i)**2 + qlep_Wrest(3,i)**2 + qlep_Wrest(4,i)**2)
      c_theta_star(i) = dot_space(qlep_Wrest(2,i),qW(2))/(lep_P_Wrest(i)*W_P_Wrest)
      enddo
***]cos_theta_star calculation

***[phi_star calculation
      do i = 1, 2
      qlep_Wrest_xyz(4,i) = lep_P_Wrest(i)*c_theta_star(i)
      enddo

      proton_lab(1,1) = ecm/2
      proton_lab(2,1) = 0
      proton_lab(3,1) = 0
      proton_lab(4,1) = proton_lab(1,1)

      proton_lab(1,2) = ecm/2
      proton_lab(2,2) = 0
      proton_lab(3,2) = 0
      proton_lab(4,2) = -proton_lab(1,2)
C Note: the y component is still 0 after this boost
      do i = 1, 2
       call BOOSTN(proton_lab(1,i),qW_inv(1),proton_Wrest(1,i))
      enddo

C warning: print this with Ncall <= 200 to avoid large-size output
c      do i = 1,2
c      print*,"proton =",i,proton_Wrest(1,i),proton_Wrest(2,i),proton_Wrest(3,i),proton_Wrest(4,i)
c      enddo
c      stop

***[find the proton pointing to the positive X direction
      do i = 1,2
      proton_P_Wrest(i) = sqrt(proton_Wrest(2,i)**2 + proton_Wrest(3,i)**2 + proton_Wrest(4,i)**2)
      temp(i) = W_P_Wrest*proton_P_Wrest(i)
      c_theta_proton(i) = dot_space(proton_Wrest(2,i),qW(2))/temp(i)
      enddo

       if ( c_theta_proton(1) .gt. c_theta_proton(2) ) then
        ipart = 1
       else
        ipart = 2
       endif
***]find the proton pointing to the positive X direction

***[find the y-axis vector
C Y = Z x X (right-handed system), Z = qW, X = proton(ipart)
       vectorY(2) = qW(3)*proton_Wrest(4,ipart) - qW(4)*proton_Wrest(3,ipart)
       vectorY(3) = qW(4)*proton_Wrest(2,ipart) - qW(2)*proton_Wrest(4,ipart)
       vectorY(4) = qW(2)*proton_Wrest(3,ipart) - qW(3)*proton_Wrest(2,ipart)
***]find the y-axis vector

***[find the X-axis vector
C X = Y x Z (right-handed system)
       vectorX(2) = vectorY(3)*qW(4) - vectorY(4)*qW(3)
       vectorX(3) = vectorY(4)*qW(2) - vectorY(2)*qW(4)
       vectorX(4) = vectorY(2)*qW(3) - vectorY(3)*qW(2)
***]find the X-axis vector

      do i = 1, 2
C component X
      temp(i) = dot_space(qlep_Wrest(2,i),vectorX(2))
c in the xyz frame, vectorX_x = |vectorX|
      qlep_Wrest_xyz(2,i) = 
     & temp(i)/sqrt(vectorX(2)**2 + vectorX(3)**2 + vectorX(4)**2)

C component Y
      temp(i) = dot_space(qlep_Wrest(2,i),vectorY(2))
c in the xyz frame, vectorY_y = |vectorY|
      qlep_Wrest_xyz(3,i) = 
     & temp(i)/sqrt(vectorY(2)**2 + vectorY(3)**2 + vectorY(4)**2)

C warning: print this with Ncall <= 200 to avoid large-size output
c      print*,"Y comp proj =",i,qlep_Wrest_xyz(3,i)
c      print*,"Y comp magn =",sqrt(lep_P_Wrest(i)**2 - qlep_Wrest_xyz(4,i)**2 - qlep_Wrest_xyz(2,i)**2)
c      stop

       phi_star(i) = atan2( qlep_Wrest_xyz(3,i), qlep_Wrest_xyz(2,i) ) ! from -pi to +pi
       if ( phi_star(i) .lt. 0D0 ) phi_star(i) = phi_star(i) + 2*pi ! from 0 to 2*pi
       enddo

C check for NaN
#ifdef CHECK_NAN
        do j = 1, 2
        if (ISNAN( phi_star(j) )) then

	call PrintPSP()
	call PrintPSP_OS()
	print*,'phi_star(j) =',j,phi_star(j), qlep_Wrest_xyz(3,j), qlep_Wrest_xyz(2,j),
     &  proton_P_Wrest(ipart), ipart, c_theta_proton(1), c_theta_proton(2)
        print*,"qW(4) =",qW(1),qW(2),qW(3),qW(4)
        do i = 1, 2
        print*,"proton_Wrest(4,i) =",i,proton_Wrest(1,i),proton_Wrest(2,i),proton_Wrest(3,i),proton_Wrest(4,i)
        enddo
        stop

c         result(1) = 0
c         return
        endif
        enddo
#endif
***]phi_star calculation
      END
***********************************************************************
	subroutine theta_vvv(iV,plep,qV_Lab,c_theta)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Refs: https://arxiv.org/abs/1103.5445 (read Fig 10)
c and also https://arxiv.org/abs/1204.6427
C Diff from theta_star(): the Z' axis is p_V in the VVV center of mass system
C Note: theta_star() is the angle between the lepton momentum (as observed 
C in the V-rest frame relative to the VVV frame) and the V momentum (as observed 
C in the VVV frame). Note that the V-rest frame is not unique. 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        use global_vars
	use global_mom_xpdf
	implicit none
	integer iV, i
	double precision c_theta(2), dot_space
	double precision plep(4,Nlepton) ! all leptons here (incl neutrinos)
	double precision qlep(4,2) ! the leptons corresp. to the selected boson
	double precision qlep_VVV(4,2), qlep_Wrest(4,2)
	double precision qV_Lab(4,3), qX(4), qX_inv(4)
	double precision qW_Lab(4) ! mom of the selected resonance
	double precision qW(4) ! mom in the CMS_VVV frame
	double precision qW_inv_VVV(4), W_P_Wrest, lep_P_Wrest(2)

#include "model.h"

***[boost from Lab frame to the VVV-c.m.s frame
! get the qW in the VVV frame
      do i = 1, 4
       qX(i) = qV_Lab(i,1) + qV_Lab(i,2) + qV_Lab(i,3) ! in the Lab frame
      enddo
      qX_inv(1) = qX(1)
      do i = 2, 4
       qX_inv(i) = -qX(i)
      enddo
      qW_Lab(:) = qV_Lab(:,iV) 
      call BOOSTN(qW_Lab(1),qX_inv(1),qW(1))
***]boost from Lab frame to the VVV-c.m.s frame

***[cos_theta_star calculation
! get the lep mom in the W-rest frame
      qW_inv_VVV(1) = qW(1)
      do i = 2, 4
      qW_inv_VVV(i) = -qW(i)
      enddo

      W_P_Wrest = sqrt(qW(2)**2 + qW(3)**2 + qW(4)**2)

!      print*,"Lab pW:",qW_Lab(:)
!      print*,"VVV  pW:",qW(:)

      do i = 1,2
       qlep(:,i) = plep(:,i+(iV-1)*2)
      enddo

C boost from the Lab (proton-proton CMS) frame to the W-rest frame
      do i = 1, 2
C boost: Lab > VVV frame > V-rest frame
      call BOOSTN(qlep(1,i),qX_inv(1),qlep_VVV(1,i))
      call BOOSTN(qlep_VVV(1,i),qW_inv_VVV(1),qlep_Wrest(1,i))
!      print*,"Lab VVV Rest qlep:",i,qlep_Wrest(:,i)

CC direct boost: Lab > V-rest frame (different results)
cc      call BOOSTN(qlep(1,i),qW_inv(1),qlep_Wrest(1,i))

      lep_P_Wrest(i) = sqrt(qlep_Wrest(2,i)**2 + qlep_Wrest(3,i)**2 + qlep_Wrest(4,i)**2)
      c_theta(i) = dot_space(qlep_Wrest(2,i),qW(2))/(lep_P_Wrest(i)*W_P_Wrest)
      enddo

!      stop
***]cos_theta_star calculation

       end subroutine theta_vvv

***********************************************************************
	subroutine theta_wz(qlep, qW_Lab, qZ_Lab, c_theta_star, phi_star)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Refs: https://arxiv.org/abs/1103.5445 (read Fig 10)
c and also https://arxiv.org/abs/1204.6427
c this routine works for both W and Z
C Diff from theta_star(): the Z' axis is p_V in the WZ center of mass system
C Note: theta_star() is the angle between the lepton momentum (as observed 
C in the V-rest frame relative to the WZ frame) and the V momentum (as observed 
C in the WZ frame). Note that the V-rest frame is not unique. 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        use global_vars
	use global_mom_xpdf
	implicit none
	double precision qW_Lab(4), qZ_Lab(4)
	integer i, j, ipart

       double precision qlep(4,2) ! only 2 decay leptons here
       double precision c_theta_star(2), phi_star(2), qW_inv(4), 
     & qlep_Wrest(4,2), dot_space
      double precision proton_lab(4,2), proton_Wrest(4,2),
     & lep_P_Wrest(2), W_P_Wrest
      double precision proton_P_Wrest(2), c_theta_proton(2),
     & qlep_Wrest_xyz(2:4,2), vectorY(2:4)
       double precision temp(2), tmp(2), vectorX(2:4)
       double precision qWZ(4), qWZ_inv(4)
	double precision qW(4) ! in the WZ c.m.s frame

       double precision qlep_WZ(4,2), qW_WZ_inv(4)

#include "model.h"

***[boost from Lab frame to the WZ-c.m.s frame
      do i = 1, 4
       qWZ(i) = qW_Lab(i) + qZ_Lab(i) ! in the Lab frame
      enddo
      qWZ_inv(1) = qWZ(1)
      do i = 2, 4
       qWZ_inv(i) = -qWZ(i)
      enddo
      call BOOSTN(qW_Lab(1),qWZ_inv(1),qW(1))
***]boost from Lab frame to the WZ-c.m.s frame

***[cos_theta_star calculation
      qW_inv(1) = qW_Lab(1)
      do i = 2, 4
      qW_inv(i) = -qW_Lab(i)
      enddo

      qW_WZ_inv(1) = qW(1)
      do i = 2, 4
      qW_WZ_inv(i) = -qW(i)
      enddo

      W_P_Wrest = sqrt(qW(2)**2 + qW(3)**2 + qW(4)**2)

c      print*,"Lab pW:",qW_Lab(:)
c      print*,"WZ  pW:",qW(:)

C boost from the Lab (proton-proton CMS) frame to the W-rest frame
      do i = 1, 2
C boost: Lab > WZ frame > V-rest frame
      call BOOSTN(qlep(1,i),qWZ_inv(1),qlep_WZ(1,i))
      call BOOSTN(qlep_WZ(1,i),qW_WZ_inv(1),qlep_Wrest(1,i))
c      print*,"Lab WZ Rest qlep:",i,qlep_Wrest(:,i)
CC direct boost: Lab > V-rest frame
cc      call BOOSTN(qlep(1,i),qW_inv(1),qlep_Wrest(1,i))
cc      print*,"Lab di Rest qlep:",i,qlep_Wrest(:,i)
      lep_P_Wrest(i) = sqrt(qlep_Wrest(2,i)**2 + qlep_Wrest(3,i)**2 + qlep_Wrest(4,i)**2)
      c_theta_star(i) = dot_space(qlep_Wrest(2,i),qW(2))/(lep_P_Wrest(i)*W_P_Wrest)
      enddo
c      stop
***]cos_theta_star calculation

***[phi_star calculation
      do i = 1, 2
      qlep_Wrest_xyz(4,i) = lep_P_Wrest(i)*c_theta_star(i)
      enddo

      proton_lab(1,1) = ecm/2
      proton_lab(2,1) = 0
      proton_lab(3,1) = 0
      proton_lab(4,1) = proton_lab(1,1)

      proton_lab(1,2) = ecm/2
      proton_lab(2,2) = 0
      proton_lab(3,2) = 0
      proton_lab(4,2) = -proton_lab(1,2)
C Note: the y component is still 0 after this boost
      do i = 1, 2
       call BOOSTN(proton_lab(1,i),qW_inv(1),proton_Wrest(1,i))
      enddo

C warning: print this with Ncall <= 200 to avoid large-size output
c      do i = 1,2
c      print*,"proton =",i,proton_Wrest(1,i),proton_Wrest(2,i),proton_Wrest(3,i),proton_Wrest(4,i)
c      enddo
c      stop

***[find the proton pointing to the positive X direction
      do i = 1,2
      proton_P_Wrest(i) = sqrt(proton_Wrest(2,i)**2 + proton_Wrest(3,i)**2 + proton_Wrest(4,i)**2)
      temp(i) = W_P_Wrest*proton_P_Wrest(i)
      c_theta_proton(i) = dot_space(proton_Wrest(2,i),qW(2))/temp(i)
      enddo

       if ( c_theta_proton(1) .gt. c_theta_proton(2) ) then
        ipart = 1
       else
        ipart = 2
       endif
***]find the proton pointing to the positive X direction

***[find the y-axis vector
C Y = Z x X (right-handed system), Z = qW, X = proton(ipart)
       vectorY(2) = qW(3)*proton_Wrest(4,ipart) - qW(4)*proton_Wrest(3,ipart)
       vectorY(3) = qW(4)*proton_Wrest(2,ipart) - qW(2)*proton_Wrest(4,ipart)
       vectorY(4) = qW(2)*proton_Wrest(3,ipart) - qW(3)*proton_Wrest(2,ipart)
***]find the y-axis vector

***[find the X-axis vector
C X = Y x Z (right-handed system)
       vectorX(2) = vectorY(3)*qW(4) - vectorY(4)*qW(3)
       vectorX(3) = vectorY(4)*qW(2) - vectorY(2)*qW(4)
       vectorX(4) = vectorY(2)*qW(3) - vectorY(3)*qW(2)
***]find the X-axis vector

      do i = 1, 2
C component X
      temp(i) = dot_space(qlep_Wrest(2,i),vectorX(2))
c in the xyz frame, vectorX_x = |vectorX|
      qlep_Wrest_xyz(2,i) = 
     & temp(i)/sqrt(vectorX(2)**2 + vectorX(3)**2 + vectorX(4)**2)

C component Y
      temp(i) = dot_space(qlep_Wrest(2,i),vectorY(2))
c in the xyz frame, vectorY_y = |vectorY|
      qlep_Wrest_xyz(3,i) = 
     & temp(i)/sqrt(vectorY(2)**2 + vectorY(3)**2 + vectorY(4)**2)

C warning: print this with Ncall <= 200 to avoid large-size output
c      print*,"Y comp proj =",i,qlep_Wrest_xyz(3,i)
c      print*,"Y comp magn =",sqrt(lep_P_Wrest(i)**2 - qlep_Wrest_xyz(4,i)**2 - qlep_Wrest_xyz(2,i)**2)
c      stop

       phi_star(i) = atan2( qlep_Wrest_xyz(3,i), qlep_Wrest_xyz(2,i) ) ! from -pi to +pi
       if ( phi_star(i) .lt. 0D0 ) phi_star(i) = phi_star(i) + 2*pi ! from 0 to 2*pi
       enddo

C check for NaN
#ifdef CHECK_NAN
        do j = 1, 2
        if (ISNAN( phi_star(j) )) then

	call PrintPSP()
	call PrintPSP_OS()
	print*,'phi_star(j) =',j,phi_star(j), qlep_Wrest_xyz(3,j), qlep_Wrest_xyz(2,j),
     &  proton_P_Wrest(ipart), ipart, c_theta_proton(1), c_theta_proton(2)
        print*,"qW(4) =",qW(1),qW(2),qW(3),qW(4)
        do i = 1, 2
        print*,"proton_Wrest(4,i) =",i,proton_Wrest(1,i),proton_Wrest(2,i),proton_Wrest(3,i),proton_Wrest(4,i)
        enddo
        stop

c         result(1) = 0
c         return
        endif
        enddo
#endif
***]phi_star calculation
      END
***********************************************************************
	subroutine theta_w_wz(qW_Lab, qZ_Lab, c_theta_w)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C calculate cos_theta_w in the WZ frame
C The positive Z axis is the direction of p1 (first proton) in the LAB frame
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	implicit none
	double precision qW_Lab(4), qZ_Lab(4)
	integer i
        double precision c_theta_w, qqW
        double precision qWZ(4), qWZ_inv(4)
	double precision qW(4) ! in the WZ c.m.s frame

c#include "xsection.h"

***[boost from Lab frame to the WZ-c.m.s frame
      do i = 1, 4
       qWZ(i) = qW_Lab(i) + qZ_Lab(i) ! in the Lab frame
      enddo
      qWZ_inv(1) = qWZ(1)
      do i = 2, 4
       qWZ_inv(i) = -qWZ(i)
      enddo
      call BOOSTN(qW_Lab(1),qWZ_inv(1),qW(1))
***]boost from Lab frame to the WZ-c.m.s frame

      qqW = sqrt(qW(2)**2 + qW(3)**2 + qW(4)**2)
      c_theta_w = qW(4)/qqW      
      END      
***********************************************************************
#ifdef USE_COLLINS_SOPER_RESULTS
	subroutine theta_cs(qlep, qW, c_theta_cs, phi_cs)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Ninh's implementation from the analytical results in the paper and the book
c Ref1: Collins, Soper PRD16, 2219 (1977) (read Introduction and Section II)
c Ref2: Book by Yorikiyo Nagashima, Foundations of the Standard Model (2013) p. 77
c this routine works for both W and Z 
c Notation: W = W, Z
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        use global_vars
	implicit none
	double precision qW(4), dot, proton_lab(4,2)
	integer i, j

       double precision qlep(4,2) ! only 2 decay leptons here

      double precision c_theta_cs(2), phi_cs(2), x, y
      double precision Q2, QT2, qlep_plus(2), qlep_minus(2)
      double precision Delta(4), vectorR(2:4), hatRT(2:3), hatQT(2:3), t_phi_cs(2)

#include "xsection.h"	

      proton_lab(1,1) = ecm/2
      proton_lab(2,1) = 0
      proton_lab(3,1) = 0
      proton_lab(4,1) = proton_lab(1,1)

      proton_lab(1,2) = ecm/2
      proton_lab(2,2) = 0
      proton_lab(3,2) = 0
      proton_lab(4,2) = -proton_lab(1,2)

***[cos_theta, phi calculation: Collins-Soper's method
c Ref1. Collin Eq. (2.1)
      Q2 = dot(qW(1),qW(1))
      QT2 = qW(2)**2 + qW(3)**2
c Ref2: Book by Yorikiyo Nagashima, Foundations of the Standard Model (2013) p. 77
      do i = 1, 2
      qlep_plus(i) = (qlep(1,i) + qlep(4,i))/sqrt2
      qlep_minus(i) = (qlep(1,i) - qlep(4,i))/sqrt2
      enddo

      c_theta_cs(2) = 2*(qlep_plus(2)*qlep_minus(1) - qlep_plus(1)*qlep_minus(2))
     & /sqrt(Q2*(Q2 + QT2))
      c_theta_cs(1) = -c_theta_cs(2)

c phi: 
      do i = 1, 4
      Delta(i) = qlep(i,2) - qlep(i,1)
      enddo
C R = P_A x Q, p_A = proton(1), Q = qW
       vectorR(2) = proton_lab(3,1)*qW(4) - proton_lab(4,1)*qW(3)
       vectorR(3) = proton_lab(4,1)*qW(2) - proton_lab(2,1)*qW(4)
       vectorR(4) = proton_lab(2,1)*qW(3) - proton_lab(3,1)*qW(2)

       do i = 2, 3
       hatRT(i) = vectorR(i)/sqrt(vectorR(2)**2 + vectorR(3)**2) ! unit transverse vector
       hatQT(i) = qW(i)/sqrt(qW(2)**2 + qW(3)**2) ! unit transverse vector
       enddo

C see (F.8) of Book by Yorikiyo Nagashima, Foundations of the Standard Model (2013)
       y = (Delta(2)*hatRT(2) + Delta(3)*hatRT(3))*sqrt((Q2 + QT2)/Q2)
       x = Delta(2)*hatQT(2) + Delta(3)*hatQT(3)
       phi_cs(2) = atan2(y, x) ! from -pi to +pi
       if ( phi_cs(2) .lt. 0D0 ) phi_cs(2) = phi_cs(2) + 2*pi ! from 0 to 2*pi
       if (phi_cs(2) .lt. pi) then
        phi_cs(1) = phi_cs(2) + pi
       else
        phi_cs(1) = phi_cs(2) - pi 
       endif
***]cos_theta, phi calculation: Collins-Soper's method

      END
#endif
#ifndef USE_COLLINS_SOPER_RESULTS
***********************************************************************
	subroutine theta_cs(qlep, qW, c_theta, phi)
c#define CHECK_WITH_COLLINS_SOPER ! watch the print out
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Ninh's implementation from the definition
c Ref1: Collins, Soper PRD16, 2219 (1977) (read Introduction and Section II)
c Ref2: Book by Yorikiyo Nagashima, Foundations of the Standard Model (2013) p. 77
c this routine works for both W and Z 
c Notation: W = W, Z
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        use global_vars
	implicit none
	double precision qW(4)
	integer i, j

       double precision qlep(4,2) ! only 2 decay leptons here
       double precision c_theta(2), phi(2), qW_inv(4), 
     & qlep_Wrest(4,2), lep_P_Wrest(2), dot_space, dot
      double precision proton_lab(4,2), proton_Wrest(4,2), oz_Wrest(2:4), oz_P_Wrest
      double precision vectorY(2:4), proton_P_Wrest(2)
      double precision c_theta_proton(2), qlep_Wrest_xyz(2:4,2), temp(2), tmp
      double precision directionX(2:4), vectorX(2:4)

      double precision c_theta_cs(2), phi_cs(2), x, y
      double precision Q2, QT2, qlep_plus(2), qlep_minus(2)
      double precision Delta(4), vectorR(2:4), hatRT(2:3), hatQT(2:3), t_phi_cs(2)

#include "model.h"	

      qW_inv(1) = qW(1)
      do i = 2, 4
      qW_inv(i) = -qW(i)
      enddo

C boost from the Lab (proton-proton CMS) frame to the W-rest frame
      do i = 1, 2
      call BOOSTN(qlep(1,i),qW_inv(1),qlep_Wrest(1,i))
      enddo

      proton_lab(1,1) = ecm/2
      proton_lab(2,1) = 0
      proton_lab(3,1) = 0
      proton_lab(4,1) = proton_lab(1,1)

      proton_lab(1,2) = ecm/2
      proton_lab(2,2) = 0
      proton_lab(3,2) = 0
      proton_lab(4,2) = -proton_lab(1,2)
C Note: the y component is still 0 after this boost
      do i = 1, 2
       call BOOSTN(proton_lab(1,i),qW_inv(1),proton_Wrest(1,i))
      enddo

***[calculate the OZ' axis in Wrest frame as seen in the Lab coordinate system
      do i = 1,2
      proton_P_Wrest(i) = sqrt(proton_Wrest(2,i)**2 
     & + proton_Wrest(3,i)**2 + proton_Wrest(4,i)**2)
      enddo

      do i = 2, 4
      oz_Wrest(i) = proton_Wrest(i,1)/proton_P_Wrest(1) 
     & - proton_Wrest(i,2)/proton_P_Wrest(2)
      enddo
      oz_P_Wrest = sqrt(oz_Wrest(2)**2 + oz_Wrest(3)**2 + oz_Wrest(4)**2)
      do i = 2, 4
      oz_Wrest(i) = oz_Wrest(i)/oz_P_Wrest ! normalized
      enddo

      tmp = dot_space(oz_Wrest(2), qW(2))
      if ( tmp .lt. 0 ) then
      do i = 2, 4
       oz_Wrest(i) = -oz_Wrest(i) ! choose the Z' direction the same as qW in lab frame
      enddo
      endif
***]calculate the OZ' axis in Wrest frame as seen in the Lab coordinate system

***[cos_theta calculation: Ninh's method
      do i = 1, 2
C boost from the Lab (proton-proton CMS) frame to the W-rest frame
      call BOOSTN(qlep(1,i),qW_inv(1),qlep_Wrest(1,i))
      lep_P_Wrest(i) = sqrt(qlep_Wrest(2,i)**2 
     & + qlep_Wrest(3,i)**2 + qlep_Wrest(4,i)**2)
      c_theta(i) = dot_space(qlep_Wrest(2,i),oz_Wrest(2))/lep_P_Wrest(i)
      enddo
***]cos_theta calculation: Ninh's method

***[phi calculation: Ninh's method
       do i = 2, 4
       directionX(i) = -proton_Wrest(i,1) - proton_Wrest(i,2) ! follow 1407.2940
       enddo

***[find the y-axis vector: OY'
C Y' = Z' x X' (right-handed system), Z' = oz_Wrest
       vectorY(2) = oz_Wrest(3)*directionX(4) - oz_Wrest(4)*directionX(3)
       vectorY(3) = oz_Wrest(4)*directionX(2) - oz_Wrest(2)*directionX(4)
       vectorY(4) = oz_Wrest(2)*directionX(3) - oz_Wrest(3)*directionX(2)
***]find the y-axis vector: OY'

***[find the X-axis vector: OX'
C X' = Y' x Z' (right-handed system)
       vectorX(2) = vectorY(3)*oz_Wrest(4) - vectorY(4)*oz_Wrest(3)
       vectorX(3) = vectorY(4)*oz_Wrest(2) - vectorY(2)*oz_Wrest(4)
       vectorX(4) = vectorY(2)*oz_Wrest(3) - vectorY(3)*oz_Wrest(2)
***]find the X-axis vector: OX'

      do i = 1, 2
      qlep_Wrest_xyz(4,i) = lep_P_Wrest(i)*c_theta(i)

C component X
      temp(i) = dot_space(qlep_Wrest(2,i),vectorX(2))
c in the xyz frame, vectorX_x = |vectorX|
      qlep_Wrest_xyz(2,i) = 
     & temp(i)/sqrt(vectorX(2)**2 + vectorX(3)**2 + vectorX(4)**2)
C component Y
      temp(i) = dot_space(qlep_Wrest(2,i),vectorY(2))
c in the xyz frame, vectorY_y = |vectorY|
      qlep_Wrest_xyz(3,i) = 
     & temp(i)/sqrt(vectorY(2)**2 + vectorY(3)**2 + vectorY(4)**2)

C warning: print this with Ncall <= 200 to avoid large-size output
c      print*,"Y comp proj =",qlep_Wrest_xyz(3,i)
c      print*,"Y comp magn =",sqrt(lep_P_Wrest(i)**2 - qlep_Wrest_xyz(4,i)**2 - qlep_Wrest_xyz(2,i)**2)
c      stop

       phi(i) = atan2( qlep_Wrest_xyz(3,i), qlep_Wrest_xyz(2,i) ) ! from -pi to +pi
       if ( phi(i) .lt. 0D0 ) phi(i) = phi(i) + 2*pi ! from 0 to 2*pi
       enddo
***]phi calculation: Ninh's method

#ifdef CHECK_WITH_COLLINS_SOPER
***[cos_theta, phi calculation: Collins-Soper's method
c Ref1. Collin Eq. (2.1)
      Q2 = dot(qW(1),qW(1))
      QT2 = qW(2)**2 + qW(3)**2
c Ref2: Book by Yorikiyo Nagashima, Foundations of the Standard Model (2013) p. 77
      do i = 1, 2
      qlep_plus(i) = (qlep(1,i) + qlep(4,i))/sqrt2
      qlep_minus(i) = (qlep(1,i) - qlep(4,i))/sqrt2
      enddo

      c_theta_cs(2) = 2*(qlep_plus(2)*qlep_minus(1) - qlep_plus(1)*qlep_minus(2))
     & /sqrt(Q2*(Q2 + QT2))
      c_theta_cs(1) = -c_theta_cs(2)

c phi: 
      do i = 1, 4
      Delta(i) = qlep(i,2) - qlep(i,1)
      enddo
C R = P_A x Q, p_A = proton(1), Q = qW
       vectorR(2) = proton_lab(3,1)*qW(4) - proton_lab(4,1)*qW(3)
       vectorR(3) = proton_lab(4,1)*qW(2) - proton_lab(2,1)*qW(4)
       vectorR(4) = proton_lab(2,1)*qW(3) - proton_lab(3,1)*qW(2)

       do i = 2, 3
       hatRT(i) = vectorR(i)/sqrt(vectorR(2)**2 + vectorR(3)**2) ! unit transverse vector
       hatQT(i) = qW(i)/sqrt(qW(2)**2 + qW(3)**2) ! unit transverse vector
       enddo

C see (F.8) of Book by Yorikiyo Nagashima, Foundations of the Standard Model (2013)
       y = (Delta(2)*hatRT(2) + Delta(3)*hatRT(3))*sqrt((Q2 + QT2)/Q2)
       x = Delta(2)*hatQT(2) + Delta(3)*hatQT(3)
       phi_cs(2) = atan2(y, x) ! from -pi to +pi
       if ( phi_cs(2) .lt. 0D0 ) phi_cs(2) = phi_cs(2) + 2*pi ! from 0 to 2*pi
       if (phi_cs(2) .lt. pi) then
        phi_cs(1) = phi_cs(2) + pi
       else
        phi_cs(1) = phi_cs(2) - pi 
       endif
***]cos_theta, phi calculation: Collins-Soper's method

***[checking
      do i = 1, 2
       print*,"c_theta_ninh, phi =", i, c_theta(i), phi(i)
       print*,"c_theta_sope, phi =", i, c_theta_cs(i), phi_cs(i)
      enddo
c      stop
***]checking
#endif

      END
#endif
