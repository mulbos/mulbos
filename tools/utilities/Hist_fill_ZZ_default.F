	subroutine Fill_ZZ(qext, FUNC)
! Note: qext are in LAB frame
	use global_vars
	implicit none
	integer i, j, Ncout, ipart

      double precision pi
      parameter (pi = 3.1415926535897932384626433832795029D0)
      double precision MWsq ! W mass squared

      double precision qext(proc_nlegs,4)
      double precision FUNC, xdis, pT2(2)
      double precision temp, tmp(2)

      double precision qlepW(4,2), qlepZ(4,2)
      double precision qW(4), pt_W, qZ(4), pt_Z, pt_miss
      double precision m_Z, y_Z, eta_Z
      double precision qWZ(4), pt_WZ
      double precision pt_lep(Nlepton), py_lep(Nlepton), y_lep(Nlepton)
      double precision DRLL_W, DyLL_W, DphiLL_W
      double precision DRLL_Z, DyLL_Z, DphiLL_Z
      double precision DRLL_35, DyLL_35, DphiLL_35
      double precision DRLL_36, DyLL_36, DphiLL_36

      double precision c_theta_star_W(2), phi_star_W(2) ! for two decay leptons
      double precision c_theta_star_Z(2), phi_star_Z(2) ! for two decay leptons
      double precision c_theta_cs_W(2), phi_cs_W(2), 
     &c_theta_cs_Z(2), phi_cs_Z(2)

      double precision phi_lep(Nlepton), phi_W, phi_Z
      double precision Dy_Ze
      double precision c_theta_wz_W(2), phi_wz_W(2) ! for two decay leptons
      double precision c_theta_wz_Z(2), phi_wz_Z(2) ! for two decay leptons
      double precision qW_atlas(4), pneu_atlas(4), qlepW_atlas(4,2) !reconstructing the neutrino and W momenta as in ATLAS
      double precision pZ_abs_space,c_theta_Z,pW_abs_space,c_theta_W
      double precision y_W,eta_W
      double precision pp_lep(Nlepton) ! magnitude of the space components (Nlepton = 4)
      double precision dot_space

! qext(proc_nlegs,4): p_1 p_2 > L_3 L_4 L_5 L_6 (order of particles)

        do i = 1, 4
         do j = 1, 2
         qlepW(i,j) = qext(j+2,i)
         qlepZ(i,j) = qext(j+4,i)
         enddo
        enddo

      do i = 1, 4
       qW(i) = qlepW(i,1) + qlepW(i,2)
       qZ(i) = qlepZ(i,1) + qlepZ(i,2)
       qWZ(i) = qW(i) + qZ(i)
      enddo
       pt_W = sqrt(qW(2)**2 + qW(3)**2)
       pt_Z = sqrt(qZ(2)**2 + qZ(3)**2)
       pt_WZ = sqrt(qWZ(2)**2 + qWZ(3)**2)

!      print*,"hello from FillHIS 2"

***[get the reconstructed neutrino and W momenta as in ATLAS-CONF-2018-034
      MWsq = mass_VV(1)**2
      call GetMomNeutrino(qlepW(1,1),qlepW(1,2),MWsq,pneu_atlas)
C for W_atlas
      do i = 1, 4
       qW_atlas(i) = qlepW(i,1) + pneu_atlas(i)
       qlepW_atlas(i,1) = qlepW(i,1)
       qlepW_atlas(i,2) = pneu_atlas(i)
      enddo
***]get the reconstructed neutrino and W momenta as in ATLAS-CONF-2018-034

!      print*,"hello from FillHIS 3"

***[check momenta
c       print*,"check momenta in FillHIS.F, FUNC=",FUNC
c       do i = 1, proc_nlegs
c        print*,'[E, px, py, pz] of particle',i
c        print*,qext(i,1),qext(i,2),qext(i,3),qext(i,4)
c        print*,'mass2=',qext(i,1)**2-qext(i,2)**2-qext(i,3)**2-qext(i,4)**2
c       enddo
c       stop
***]check momenta

      do i = 1, Nlepton
      ipart = i + 2
      xdis=qext(ipart,2)**2+qext(ipart,3)**2+qext(ipart,4)**2
      pp_lep(i)=sqrt(xdis)
!      print*,"pp_lep =",i,pp_lep(i)
      enddo
!      print*,"FillHIS individual leptons"
***[individual leptons
       pt_miss=sqrt(qext(4,2)**2+qext(4,3)**2)

       do i = 1, Nlepton
       ipart = i + 2
       pt_lep(i) = sqrt(qext(ipart,2)**2 + qext(ipart,3)**2)
       enddo

       do i = 1, Nlepton
       ipart = i + 2
        xdis=qext(ipart,2)**2+qext(ipart,3)**2+qext(ipart,4)**2
	xdis=sqrt(xdis)
	py_lep(i)=log((xdis+qext(ipart,4))/(xdis-qext(ipart,4)))/2D0
        y_lep(i)=py_lep(i)
        phi_lep(i) = atan2(qext(ipart,3), qext(ipart,2)) ! from -pi to +pi
       enddo
***]individual leptons      
!      print*,"FillHIS gauge bosons"
***[gauge bosons
        phi_W = atan2(qW(3), qW(2)) ! from -pi to +pi
        phi_Z = atan2(qZ(3), qZ(2)) ! from -pi to +pi

        pZ_abs_space=sqrt(qZ(2)**2+qZ(3)**2+qZ(4)**2)
	c_theta_Z=qZ(4)/pZ_abs_space
	y_Z=log((qZ(1)+qZ(4))/(qZ(1)-qZ(4)))/2D0
	eta_Z=(pZ_abs_space-qZ(4))/(pZ_abs_space+qZ(4))
        if (eta_Z .le. 1D-10) then
          eta_Z = 1D2 ! infinity
        else
	  eta_Z=-log(eta_Z)/2D0 
        endif

        pW_abs_space=sqrt(qW(2)**2+qW(3)**2+qW(4)**2)
        c_theta_W=qW(4)/pW_abs_space
        y_W=log((qW(1) + qW(4))/(qW(1) - qW(4)))/2D0
C        eta_W=log((pW_abs_space+qW(4))/(pW_abs_space-qW(4)))/2D0 ! can be singular here (e.g. pW_x, pW_y very small)
***]gauge bosons
!      print*,"FillHIS histograms and scatter plot"
***{ Fill histograms and scatter plot
      Ncout = 1
! azimuthal separation between e and mu_plus
! CALL XHINIT(Ncout,0D0,pi,Nbin,'Distribution Dphi(3,5)')
       DphiLL_35 = abs(phi_lep(1) - phi_lep(3))
       if ( DphiLL_35 .gt. pi ) DphiLL_35 = 2*pi - DphiLL_35
       CALL XHFILL(Ncout,DphiLL_35,FUNC)
       Ncout = Ncout + 1
       
! azimuthal separation between e and mu_minus
! CALL XHINIT(Ncout,0D0,pi,Nbin,'Distribution Dphi(3,6)')
       DphiLL_36 = abs(phi_lep(1) - phi_lep(4))
       if ( DphiLL_36 .gt. pi ) DphiLL_36 = 2*pi - DphiLL_36
       CALL XHFILL(Ncout,DphiLL_36,FUNC)
       Ncout = Ncout + 1
       
! azimuthal separation between mu_plus and mu_minus
! CALL XHINIT(Ncout,0D0,pi,Nbin,'Distribution Dphi(5,6)')       
       DphiLL_Z = abs(phi_lep(3) - phi_lep(4))
       if ( DphiLL_Z .gt. pi ) DphiLL_Z = 2*pi - DphiLL_Z
       CALL XHFILL(Ncout,DphiLL_Z,FUNC)
       Ncout = Ncout + 1
       
! azimuthal separation between e and missing transverse momentum (neutrino)
! CALL XHINIT(Ncout,0D0,pi,Nbin,'Distribution Dphi(3,4)')
       DphiLL_W = abs(phi_lep(1) - phi_lep(2))
       if ( DphiLL_W .gt. pi ) DphiLL_W = 2*pi - DphiLL_W
       CALL XHFILL(Ncout,DphiLL_W,FUNC)
       Ncout = Ncout + 1
       
!       print*,"FillHIS upper mid"
       
! cosine of the angle between e and mu_plus:
! CALL XHINIT(Ncout,-1D0,1D0,Nbin,'Distribution cos(theta_35)')
      xdis = dot_space(qlepW(2,1),qlepZ(2,1))/(pp_lep(1)*pp_lep(3))
      CALL XHFILL(Ncout,xdis,FUNC)
      Ncout = Ncout + 1
      
!      print*,"FillHIS upper mid 0",pp_lep(1),pp_lep(4)

! cosine of the angle between e and mu_minus:
! CALL XHINIT(Ncout,-1D0,1D0,Nbin,'Distribution cos(theta_36)')       
      xdis = dot_space(qlepW(2,1),qlepZ(2,2))/(pp_lep(1)*pp_lep(4))
      CALL XHFILL(Ncout,xdis,FUNC)
      Ncout = Ncout + 1 
      
!      print*,"FillHIS upper mid 05"
      
! cosine of the theta_wz angle of e (W frame) and W boson (WZ rest frame), truth level
! CALL XHINIT(Ncout,-1D0,1D0,Nbin,'Distribution cos(theta_wz_W_e)')  
       call theta_wz(qlepW_atlas, qW_atlas, qZ, c_theta_wz_W(1), phi_wz_W(1)) ! for electron
       CALL XHFILL(Ncout,c_theta_wz_W(1),FUNC)
       Ncout = Ncout + 1    

!       print*,"FillHIS upper mid 1"
       
! cosine of the theta_wz angle of mu_minus (Z frame) and Z boson (WZ rest frame), truth level
! CALL XHINIT(Ncout,-1D0,1D0,Nbin,'Distribution cos(theta_wz_Z_mu6)')       
       call theta_wz(qlepZ, qZ, qW_atlas, c_theta_wz_Z(1), phi_wz_Z(1))
       CALL XHFILL(Ncout,c_theta_wz_Z(2),FUNC) ! for mu6 = mu-
       Ncout = Ncout + 1
       
!       print*,"FillHIS upper mid 2"
       
! invariant mass of the four-lepton system, truth level: [0, 2000 GeV]
! CALL XHINIT(Ncout,0D0,2D3,Nbin,'Distribution m(3,4,5,6)')       
       xdis=qWZ(1)**2-qWZ(2)**2-qWZ(3)**2-qWZ(4)**2
       xdis=sqrt(xdis)
       CALL XHFILL(Ncout,xdis,FUNC)	
       Ncout = Ncout + 1
       
! rapidity of the electron: [−2.5, 2.5]
! CALL XHINIT(Ncout,-2.5D0,2.5D0,Nbin,'Distribution y(e)')       
       CALL XHFILL(Ncout,y_lep(1),FUNC)
       Ncout = Ncout + 1       

!      print*,"FillHIS histograms middle"

! rapidity of the neutrino: [−2.5, 2.5]
! CALL XHINIT(Ncout,-2.5D0,2.5D0,Nbin,'Distribution y(neu)')
       CALL XHFILL(Ncout,y_lep(2),FUNC)
       Ncout = Ncout + 1       
       
! rapidity of the mu_plus: [−2.5, 2.5]
! CALL XHINIT(Ncout,-2.5D0,2.5D0,Nbin,'Distribution y(mu+)')
       CALL XHFILL(Ncout,y_lep(3),FUNC)
       Ncout = Ncout + 1       
       
! rapidity of the mu_minus: [−2.5, 2.5]
! CALL XHINIT(Ncout,-2.5D0,2.5D0,Nbin,'Distribution y(mu-)')       
       CALL XHFILL(Ncout,y_lep(4),FUNC)
       Ncout = Ncout + 1   
       
! rapidity of the Z boson: [−2.5, 2.5]
! CALL XHINIT(Ncout,-2.5D0,2.5D0,Nbin,'Distribution y(Z)')           
       CALL XHFILL(Ncout,y_Z,FUNC)
       Ncout = Ncout + 1
       
! rapidity of the W boson: [−2.5, 2.5], truth level
! CALL XHINIT(Ncout,-2.5D0,2.5D0,Nbin,'Distribution y(W)')       
       xdis=log((qW(1)+qW(4))/(qW(1)-qW(4)))/2D0
       CALL XHFILL(Ncout,xdis,FUNC)
       Ncout = Ncout + 1
       
! rapidity of WZ system, truth level: [−2.5, 2.5]
! CALL XHINIT(Ncout,-2.5D0,2.5D0,Nbin,'Distribution y(WZ)')       
       xdis=log((qWZ(1)+qWZ(4))/(qWZ(1)-qWZ(4)))/2D0
       CALL XHFILL(Ncout,xdis,FUNC)
       Ncout = Ncout + 1
       
! rapidity separation (in absolute value) between e and mu_plus: [0, 5]
! CALL XHINIT(Ncout,0D0,5D0,Nbin,'Distribution |Dy(3,5)|')
       DyLL_35 = abs(y_lep(3) - y_lep(1))
       CALL XHFILL(Ncout,DyLL_35,FUNC)
       Ncout = Ncout + 1       

! rapidity separation (in absolute value) between e and mu_minus: [0, 5]
! CALL XHINIT(Ncout,0D0,5D0,Nbin,'Distribution |Dy(3,6)|')
       DyLL_36 = abs(y_lep(4) - y_lep(1))
       CALL XHFILL(Ncout,DyLL_36,FUNC)
       Ncout = Ncout + 1
       
!       print*,"FillHIS mid end"
       
! scattering angle (in absolute value) of the single boson in the WZ CM frame computed w.r.t.
! the positive Z-axis direction of the LAB, truth level: [0, 1]
! CALL XHINIT(Ncout,0D0,1D0,Nbin,'Distribution |cos(W_WZ)|')
       CALL theta_w_wz(qW_atlas, qZ, xdis)
       CALL XHFILL(Ncout,abs(xdis),FUNC) ! absolute value
       Ncout = Ncout + 1  
       
! transverse mass of the W boson: [0, 200 GeV]
! CALL XHINIT(Ncout,0D0,2D2,Nbin,'Distribution mT_W')            
      xdis = cos(phi_lep(1) - phi_lep(2))
      xdis = sqrt(2*pt_lep(1)*pt_lep(2)*(1 - xdis))
      CALL XHFILL(Ncout,xdis,FUNC)
      Ncout = Ncout + 1
      
! transverse mass of the WZ system: [0, 1000 GeV]
! CALL XHINIT(Ncout,0D0,1D3,Nbin,'Distribution mT_WZ') 
! qext(4,mu): neutrino momentum
      temp = 0
      tmp(1) = 0
      tmp(2) = 0
      do i = 1, Nlepton
       temp = temp + pt_lep(i)
       ipart = i + 2
       tmp(1) = tmp(1) + qext(ipart,2) ! p_x
       tmp(2) = tmp(2) + qext(ipart,3) ! p_y
      enddo
      xdis = sqrt(temp**2 - tmp(1)**2 - tmp(2)**2)
      CALL XHFILL(Ncout,xdis,FUNC)
      Ncout = Ncout + 1 
      
! transverse momentum of electron: [0, 500 GeV]
! CALL XHINIT(Ncout,0D0,500D0,Nbin,'Distribution pT(e)')  
       CALL XHFILL(Ncout,pt_lep(1),FUNC)
       Ncout = Ncout + 1
       
! transverse momentum of muon plus: [0, 500 GeV]
! CALL XHINIT(Ncout,0D0,500D0,Nbin,'Distribution pT(mu+)')           
       CALL XHFILL(Ncout,pt_lep(3),FUNC)
       Ncout = Ncout + 1
       
! transverse momentum of muon minus: [0, 500 GeV]
! CALL XHINIT(Ncout,0D0,500D0,Nbin,'Distribution pT(mu-)')       
       CALL XHFILL(Ncout,pt_lep(4),FUNC)
       Ncout = Ncout + 1      

! transverse momentum of neutrino: [0, 500 GeV]
! CALL XHINIT(Ncout,0D0,500D0,Nbin,'Distribution pT(neu)')
       CALL XHFILL(Ncout,pt_lep(2),FUNC)
       Ncout = Ncout + 1      

! transverse momentum of the W boson: [0, 700 GeV]
! CALL XHINIT(Ncout,0D0,700D0,Nbin,'Distribution pT(W)')
       CALL XHFILL(Ncout,pt_W,FUNC)
       Ncout = Ncout + 1
       
! transverse momentum of the Z boson: [0, 700 GeV]
! CALL XHINIT(Ncout,0D0,700D0,Nbin,'Distribution pT(Z)')       
       CALL XHFILL(Ncout,pt_Z,FUNC)
       Ncout = Ncout + 1
       
! rapidity separation (in absolute value) between e and the Z boson (muons system): [0, 5] 
! CALL XHINIT(Ncout,0D0,5D0,Nbin,'Distribution Dy(e,Z)')       
       Dy_Ze = abs(y_Z - y_lep(1))
       CALL XHFILL(Ncout,Dy_Ze,FUNC)
       Ncout = Ncout + 1
       
!       print*,"FillHIS near end"
       
! rapidity separation (in absolute value) between W boson (truth level) and Z: [0, 5]
! CALL XHINIT(Ncout,0D0,5D0,Nbin,'Distribution Dy(W,Z)')       
       xdis = abs(y_Z - y_W)
       CALL XHFILL(Ncout,xdis,FUNC)
       Ncout = Ncout + 1
       
! transverse momentum of the WZ system: [0, 200 GeV]
! CALL XHINIT(Ncout,0D0,200D0,Nbin,'Distribution pT_WZ')       
       CALL XHFILL(Ncout,pt_WZ,FUNC)
       Ncout = Ncout + 1

! R_21 ratio = pT_V2/pT_V1 (where V1 is the boson with hardest pT): [0, 1]
! CALL XHINIT(Ncout,0D0,1D0,Nbin,'Distribution pT_V2/pT_V1') 
       xdis = pt_Z/pt_W
       if (xdis .gt. 1D0) then
        xdis = 1D0/xdis
       endif   
       CALL XHFILL(Ncout,xdis,FUNC)
       Ncout = Ncout + 1
       
! (pT_ep+pT_em+pT_mup+pT_mum)/4
! CALL XHINIT(Ncout,0D0,500D0,Nbin,'Distribution average of pTs(lep) ')
       xdis = (pt_lep(1) + pt_lep(2) + pt_lep(3) + pt_lep(4))/4D0
       CALL XHFILL(Ncout,xdis,FUNC)
       !Ncout = Ncout + 1
***} Fill histograms and scatter plot
!      print*,"FillHIS end"
      end subroutine
