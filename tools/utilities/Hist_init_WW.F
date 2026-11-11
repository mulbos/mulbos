      subroutine hist_init_WW()
	use global_dists
	implicit none
      double precision pi
      parameter (pi = 3.1415926535897932384626433832795029D0)
      integer Ncout

C ATLAS-2022 wish list
! particle ordering: p(1) p(2) > e(3) nu_e(4) mu+(5) mu-(6) 
      Ncout = 1
      ! azimuthal separation between e_3 and mu_5
      CALL XHINIT(Ncout,0D0,pi,Nbin,'Distribution Dphi(3,5)')
      Ncout = Ncout + 1
      
      ! azimuthal separation between e and missing transverse momentum
      CALL XHINIT(Ncout,0D0,pi,Nbin,'Distribution Dphi(3,miss)')
      Ncout = Ncout + 1      
      
      ! azimuthal separation between mu and missing transverse momentum
      CALL XHINIT(Ncout,0D0,pi,Nbin,'Distribution Dphi(5,miss)')
      Ncout = Ncout + 1

      ! cosine of the angle between e_3 and mu_5:
      CALL XHINIT(Ncout,-1D0,1D0,Nbin,'Distribution cos(theta_35)')
      Ncout = Ncout + 1
      
      ! cosine of the theta_wz angle of e_3 (W frame) and W boson (WZ rest frame), truth level
      CALL XHINIT(Ncout,-1D0,1D0,Nbin,'Distribution cos(theta_wz_W_e3)')
      Ncout = Ncout + 1
      
      ! cosine of the theta_wz angle of mu_5 (Z frame) and Z boson (WZ rest frame)
      CALL XHINIT(Ncout,-1D0,1D0,Nbin,'Distribution cos(theta_wz_Z_mu5)')
      Ncout = Ncout + 1      

!!![phi_wz distribution
! for phi_wz distribution using WZ-helicity frame
      CALL XHINIT(Ncout,0D0,2*pi,Nbin,'Distribution phi_wz_W_e3(0:2pi)')
      Ncout = Ncout + 1
      CALL XHINIT(Ncout,0D0,2*pi,Nbin,'Distribution phi_wz_Z_mu5(0:2pi)')
      Ncout = Ncout + 1
      CALL XHINIT(Ncout,-pi,pi,Nbin,'Distribution phi_wz_W_e3(-pi:pi)')
      Ncout = Ncout + 1
      CALL XHINIT(Ncout,-pi,pi,Nbin,'Distribution phi_wz_Z_mu5(-pi:pi)')
      Ncout = Ncout + 1
!!!]phi_wz distribution
      
      ! invariant mass of the four-lepton system, truth level: [0, 2000 GeV]
      CALL XHINIT(Ncout,0D0,2D3,Nbin,'Distribution m(3,4,5,6)')
      Ncout = Ncout + 1
      
      ! invariant mass of W+
      CALL XHINIT(Ncout,55D0,105D0,Nbin,'Distribution m(3,4)')
      Ncout = Ncout + 1
      
      ! invariant mass of e-mu system
      CALL XHINIT(Ncout,0D0,1D3,Nbin,'Distribution m(3,5)')
      Ncout = Ncout + 1
      
      ! rapidity of the electron_3: [−2.5, 2.5]
      CALL XHINIT(Ncout,-2.5D0,2.5D0,Nbin,'Distribution y(e3)')
      Ncout = Ncout + 1
      
      ! rapidity of the neutrino_4: [−2.5, 2.5]
      CALL XHINIT(Ncout,-2.5D0,2.5D0,Nbin,'Distribution y(neu4)')
      Ncout = Ncout + 1      
      
      ! rapidity of the mu_5: [−2.5, 2.5]
      CALL XHINIT(Ncout,-2.5D0,2.5D0,Nbin,'Distribution y(mu5)')
      Ncout = Ncout + 1
      
      ! rapidity of the neutrino_6: [−2.5, 2.5]
      CALL XHINIT(Ncout,-2.5D0,2.5D0,Nbin,'Distribution y(neu6)')
      Ncout = Ncout + 1
      
      ! rapidity of the W2 boson: [−2.5, 2.5]
      CALL XHINIT(Ncout,-2.5D0,2.5D0,Nbin,'Distribution y(W2)')
      Ncout = Ncout + 1
      
      ! rapidity of the W1 boson: [−2.5, 2.5], truth level
      CALL XHINIT(Ncout,-2.5D0,2.5D0,Nbin,'Distribution y(W1)')
      Ncout = Ncout + 1
      
      ! rapidity of WW system, truth level: [−2.5, 2.5]
      CALL XHINIT(Ncout,-2.5D0,2.5D0,Nbin,'Distribution y(WW)')
      Ncout = Ncout + 1
      
      ! rapidity separation (in absolute value) between e3 and mu5: [0, 5]
      CALL XHINIT(Ncout,0D0,5D0,Nbin,'Distribution |Dy(3,5)|')
      Ncout = Ncout + 1
      
      ! scattering angle (in absolute value) of the single boson in the WW CM frame computed w.r.t.
!the positive Z-axis direction of the LAB, truth level: [0, 1]
      CALL XHINIT(Ncout,0D0,1D0,Nbin,'Distribution |cos(W1_WW)|')
      Ncout = Ncout + 1
      
      ! transverse momentum of electron3: [0, 500 GeV]
      CALL XHINIT(Ncout,0D0,500D0,Nbin,'Distribution pT(3)')
      Ncout = Ncout + 1                  
      
      ! transverse momentum of muon5: [0, 500 GeV]
      CALL XHINIT(Ncout,0D0,500D0,Nbin,'Distribution pT(5)')
      Ncout = Ncout + 1
      
      ! transverse momentum of missing energy (neutrinos): [0, 500 GeV]
      CALL XHINIT(Ncout,0D0,500D0,Nbin,'Distribution pT(miss)')
      Ncout = Ncout + 1
      
      ! transverse momentum of the W1 boson: [0, 500 GeV]
      CALL XHINIT(Ncout,0D0,500D0,Nbin,'Distribution pT(W1)')
      Ncout = Ncout + 1      
      
      ! transverse momentum of the W2 boson: [0, 500 GeV]
      CALL XHINIT(Ncout,0D0,500D0,Nbin,'Distribution pT(W2)')
      Ncout = Ncout + 1
      
      ! rapidity separation (in absolute value) between e3 and the W2 boson: [0, 5] 
      CALL XHINIT(Ncout,0D0,5D0,Nbin,'Distribution Dy(3,W2)')
      Ncout = Ncout + 1
      
      ! rapidity separation (in absolute value) between W1 and W2: [0, 5]
      CALL XHINIT(Ncout,0D0,5D0,Nbin,'Distribution Dy(W1,W2)')
      Ncout = Ncout + 1
      
      ! transverse momentum of the WW system: [0, 200 GeV]
      CALL XHINIT(Ncout,0D0,200D0,Nbin,'Distribution pT_WW')
      Ncout = Ncout + 1
      
      ! R_21 ratio = pT_V2/pT_V1 (where V1 is the boson with hardest pT): [0, 1]
      CALL XHINIT(Ncout,0D0,1D0,Nbin,'Distribution pT_V2/pT_V1')
      Ncout = Ncout + 1
      
      ! (pT_e+pT_mu)/2
      CALL XHINIT(Ncout,0D0,500D0,Nbin,'Distribution average of pTs(lep) ')
      Ncout = Ncout + 1
      
!!![helicity coordinate system
      ! cosine of the theta_star angle of e_3 (W frame) and W boson (Lab frame), truth level
      CALL XHINIT(Ncout,-1D0,1D0,Nbin,'Distribution cos(theta_star_W_e3)')
      Ncout = Ncout + 1
      
      ! cosine of the theta_star angle of mu_5 (Z frame) and Z boson (Lab frame)
      CALL XHINIT(Ncout,-1D0,1D0,Nbin,'Distribution cos(theta_star_Z_mu5)')
      Ncout = Ncout + 1      
      
!!![phi_star distribution
      CALL XHINIT(Ncout,0D0,2*pi,Nbin,'Distribution phi_star_W_e3(0:2pi)')
      Ncout = Ncout + 1
      CALL XHINIT(Ncout,0D0,2*pi,Nbin,'Distribution phi_star_Z_mu5(0:2pi)')
      Ncout = Ncout + 1
      CALL XHINIT(Ncout,-pi,pi,Nbin,'Distribution phi_star_W_e3(-pi:pi)')
      Ncout = Ncout + 1
      CALL XHINIT(Ncout,-pi,pi,Nbin,'Distribution phi_star_Z_mu5(-pi:pi)')
      ! Ncout = Ncout + 1
!!!]phi_star distribution
!!!]helicity coordinate system

      End subroutine
