! Ref0: Dittmaier 1507.04676

	double precision function FuncDPsing_decayW(fac)
	use global_vars
	implicit none

#include "model.h"

	double precision GdW
	external GdW

	double precision Jacb, fac
!	double precision qV(4), MVint2(2)
	
C the off-shell gauge-boson invariant masses
C Note: pext(6,4) are the true qq -> 4l + gamma momenta
!        do i = 1, 4
!         qV(i) = pext(3,i) + pext(4,i) + pext(7,i)
!        enddo
!	MVint2(1) = qV(1)**2 - qV(2)**2 - qV(3)**2 - qV(4)**2

!	MVint2(2) = pext(5,1)*pext(6,1) - pext(5,2)*pext(6,2) 
!     &  - pext(5,3)*pext(6,3) - pext(5,4)*pext(6,4)
!        MVint2(2) = MVint2(2)*2

	Jacb = -1D0*fac

      if (lNLO_QCD) then
       !Jacb = -16*pi*Jacb/3D0
       print*,"error: hadronic decay of W not implemented"
       stop
      else
       Jacb = -4*pi*Jacb
      endif

       FuncDPsing_decayW = Jacb*GdW(Jacb)
	end

************************************************************************
	double precision function GdW(Jacb)
        use global_vars
        use global_mom_xpdf
	implicit none
#include "model.h"

	integer i, a, k, n, i1

        double precision Jacb
	double precision ma2, bb
	double precision temp, temp1
	double precision Pik, Sik, Pin, Pnk
	double precision yia, zia, gia, Ria

        double precision tp(2,4) ! W > e6 + nu_e7 + phot, phot = oext(5)
!        double precision MVint2(2) ! the off-shell invariant mass of the gauge bosons (N+1 kin)
!        double precision MVint2_x(2) ! the off-shell invariant mass of the gauge bosons (reduced kin)
      integer i_off, n_off, k_off
      double precision tp_off(2,4), pW_off(4)
      double precision ma2_off, bb_off
      double precision PiK_off,Sik_off,Pin_off,Pnk_off,yia_off,zia_off
      double precision Ria_off, gia_off

! oext: q(1)q(2) > W(3)Z(4)gamma(5) > e(6)nu_e(7)mu+(8)mu-(9)
! pext: q(1)q(2) > e(3)nu_e(4)mu+(5)mu-(6)gamma(7)

*** convention: i = final, a = initial ***
        k = 5 ! on-shell photon
        k_off = 7 ! off-shell photon
      if (lDECAYV1_NLO_EW) then
        a = 3 ! W
        i = 6 ! electron
        n = 7 ! neutrino
! off-shell correspondence
        i_off = 3
        n_off = 4
      else
        a = 4 ! W
        i = 8 ! muon
        n = 9 ! neutrino
! off-shell correspondence
        i_off = 5
        n_off = 6
      endif

        ma2 = MW2

        Pik = oext(i,1)*oext(k,1)-oext(i,2)*oext(k,2)
     &  -oext(i,3)*oext(k,3)-oext(i,4)*oext(k,4)
        Sik = 2*Pik

        Pin = oext(i,1)*oext(n,1)-oext(i,2)*oext(n,2)
     &  -oext(i,3)*oext(n,3)-oext(i,4)*oext(n,4)

        Pnk = oext(n,1)*oext(k,1)-oext(n,2)*oext(k,2)
     &  -oext(n,3)*oext(k,3)-oext(n,4)*oext(k,4)

C Eq. (3.15) of Ref0
        yia = Sik/ma2
        zia = Pin/(Pin + Pnk)

        bb = ma2/(ma2 - Sik)

	do i1 = 1,4
        tp(2,i1)=bb*oext(n,i1) ! neutrino
	tp(1,i1)=oext(a,i1) - tp(2,i1) ! electron
	enddo
	
! must use the on-shell PiK here, 
! Note: PiK = yia*ma2/2 = yia_off*ma2/2, 
! yia is the integration variable when integrating the dipole term,
! ma2 = MW2 = 80^2 to match the on-shell W > e nu gamma decay
        Ria = 1 - yia	
	gia = (2d0/(1 - zia*Ria)-1-zia)/(Pik*Ria)
        gia = gia-4d0/(ma2*Ria*(1-zia*Ria)**2)

!!![off-shell correspondence
	do i1 = 1,4
        pW_off(i1) = pext(i_off,i1) + pext(n_off,i1) + pext(k_off,i1)
        enddo
        ma2_off = pW_off(1)**2-pW_off(2)**2-pW_off(3)**2-pW_off(4)**2

      Pik_off=pext(i_off,1)*pext(k_off,1)-pext(i_off,2)*pext(k_off,2)
     &  -pext(i_off,3)*pext(k_off,3)-pext(i_off,4)*pext(k_off,4)
      Sik_off = 2*Pik_off

      Pin_off=pext(i_off,1)*pext(n_off,1)-pext(i_off,2)*pext(n_off,2)
     &  -pext(i_off,3)*pext(n_off,3)-pext(i_off,4)*pext(n_off,4)

      Pnk_off=pext(n_off,1)*pext(k_off,1)-pext(n_off,2)*pext(k_off,2)
     &  -pext(n_off,3)*pext(k_off,3)-pext(n_off,4)*pext(k_off,4)

C Eq. (3.15) of Ref0
        yia_off = Sik_off/ma2_off
        zia_off = Pin_off/(Pin_off + Pnk_off)
        bb_off = ma2_off/(ma2_off - Sik_off)

! should be equal if OS mapping Denner ZZ 2021 is 
! used for W > e nu_e gamma decay 
c        print*,"yia_off =",yia, yia_off
c        print*,"zia_off =",zia, zia_off
c        print*,"bb_off =",bb, bb_off

	do i1 = 1,4
        tp_off(2,i1)=bb_off*pext(n_off,i1) ! neutrino
	tp_off(1,i1)=pW_off(i1) - tp_off(2,i1) ! electron
	enddo
!!!]off-shell correspondence

      call DoTildeAmpHis_decayW_OffS(tp,tp_off,1,gia,Jacb)
c      print*,"gia =", gia

        GdW = gia
	end
