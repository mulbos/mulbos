      subroutine hist_init_ZZ()
	use global_dists
	implicit none
      double precision pi
      parameter (pi = 3.1415926535897932384626433832795029D0)
      integer Ncout

C ATLAS-2022 wish list
! particle ordering: p(1) p(2) > e(3) nu_e(4) mu+(5) mu-(6) 
      Ncout = 1
      ! azimuthal separation between e and mu_plus
      CALL XHINIT(Ncout,0D0,pi,Nbin,'Distribution Dphi(3,5)')
      Ncout = Ncout + 1
      
      ! azimuthal separation between e and mu_minus
      CALL XHINIT(Ncout,0D0,pi,Nbin,'Distribution Dphi(3,6)')
      Ncout = Ncout + 1
      
      ! azimuthal separation between mu_plus and mu_minus
      CALL XHINIT(Ncout,0D0,pi,Nbin,'Distribution Dphi(5,6)')
      Ncout = Ncout + 1      
      
      ! azimuthal separation between e and missing transverse momentum
      CALL XHINIT(Ncout,0D0,pi,Nbin,'Distribution Dphi(3,4)')
      Ncout = Ncout + 1      

      ! cosine of the angle between e and mu_plus:
      CALL XHINIT(Ncout,-1D0,1D0,Nbin,'Distribution cos(theta_35)')
      Ncout = Ncout + 1
      
      ! cosine of the angle between e and mu_minus:
      CALL XHINIT(Ncout,-1D0,1D0,Nbin,'Distribution cos(theta_36)')
      Ncout = Ncout + 1
      
      ! cosine of the theta_wz angle of e (W frame) and W boson (WZ rest frame), truth level
      CALL XHINIT(Ncout,-1D0,1D0,Nbin,'Distribution cos(theta_wz_W_e)')
      Ncout = Ncout + 1
      
      ! cosine of the theta_wz angle of mu_minus (Z frame) and Z boson (WZ rest frame), truth level
      CALL XHINIT(Ncout,-1D0,1D0,Nbin,'Distribution cos(theta_wz_Z_mu6)')
      Ncout = Ncout + 1      
      
      ! invariant mass of the four-lepton system, truth level: [0, 2000 GeV]
      CALL XHINIT(Ncout,0D0,2D3,Nbin,'Distribution m(3,4,5,6)')
      Ncout = Ncout + 1
      
      ! rapidity of the electron: [−2.5, 2.5]
      CALL XHINIT(Ncout,-2.5D0,2.5D0,Nbin,'Distribution y(e)')
      Ncout = Ncout + 1
      
      ! rapidity of the neutrino: [−2.5, 2.5]
      CALL XHINIT(Ncout,-2.5D0,2.5D0,Nbin,'Distribution y(neu)')
      Ncout = Ncout + 1      
      
      ! rapidity of the mu_plus: [−2.5, 2.5]
      CALL XHINIT(Ncout,-2.5D0,2.5D0,Nbin,'Distribution y(mu+)')
      Ncout = Ncout + 1
      
      ! rapidity of the mu_minus: [−2.5, 2.5]
      CALL XHINIT(Ncout,-2.5D0,2.5D0,Nbin,'Distribution y(mu-)')
      Ncout = Ncout + 1
      
      ! rapidity of the Z boson: [−2.5, 2.5]
      CALL XHINIT(Ncout,-2.5D0,2.5D0,Nbin,'Distribution y(Z)')
      Ncout = Ncout + 1
      
      ! rapidity of the W boson: [−2.5, 2.5], truth level
      CALL XHINIT(Ncout,-2.5D0,2.5D0,Nbin,'Distribution y(W)')
      Ncout = Ncout + 1
      
      ! rapidity of WZ system, truth level: [−2.5, 2.5]
      CALL XHINIT(Ncout,-2.5D0,2.5D0,Nbin,'Distribution y(WZ)')
      Ncout = Ncout + 1
      
      ! rapidity separation (in absolute value) between e and mu_plus: [0, 5]
      CALL XHINIT(Ncout,0D0,5D0,Nbin,'Distribution |Dy(3,5)|')
      Ncout = Ncout + 1
      
      ! rapidity separation (in absolute value) between e and mu_minus: [0, 5]
      CALL XHINIT(Ncout,0D0,5D0,Nbin,'Distribution |Dy(3,6)|')
      Ncout = Ncout + 1
      
      ! scattering angle (in absolute value) of the single boson in the WZ CM frame computed w.r.t.
!the positive Z-axis direction of the LAB, truth level: [0, 1]
      CALL XHINIT(Ncout,0D0,1D0,Nbin,'Distribution |cos(W_WZ)|')
      Ncout = Ncout + 1
      
      ! transverse mass of the W boson: [0, 200 GeV]
      CALL XHINIT(Ncout,0D0,2D2,Nbin,'Distribution mT_W')
      Ncout = Ncout + 1
      
      ! transverse mass of the WZ system: [0, 1000 GeV]
      CALL XHINIT(Ncout,0D0,1D3,Nbin,'Distribution mT_WZ')
      Ncout = Ncout + 1            
      
      ! transverse momentum of electron: [0, 500 GeV]
      CALL XHINIT(Ncout,0D0,500D0,Nbin,'Distribution pT(e)')
      Ncout = Ncout + 1                  
      
      ! transverse momentum of muon plus: [0, 500 GeV]
      CALL XHINIT(Ncout,0D0,500D0,Nbin,'Distribution pT(mu+)')
      Ncout = Ncout + 1
      
      ! transverse momentum of muon minus: [0, 500 GeV]
      CALL XHINIT(Ncout,0D0,500D0,Nbin,'Distribution pT(mu-)')
      Ncout = Ncout + 1
      
      ! transverse momentum of neutrino: [0, 500 GeV]
      CALL XHINIT(Ncout,0D0,500D0,Nbin,'Distribution pT(neu)')
      Ncout = Ncout + 1
      
      ! transverse momentum of the W boson: [0, 700 GeV]
      CALL XHINIT(Ncout,0D0,700D0,Nbin,'Distribution pT(W)')
      Ncout = Ncout + 1      
      
      ! transverse momentum of the Z boson: [0, 700 GeV]
      CALL XHINIT(Ncout,0D0,700D0,Nbin,'Distribution pT(Z)')
      Ncout = Ncout + 1
      
      ! rapidity separation (in absolute value) between e and the Z boson (muons system): [0, 5] 
      CALL XHINIT(Ncout,0D0,5D0,Nbin,'Distribution Dy(e,Z)')
      Ncout = Ncout + 1
      
      ! rapidity separation (in absolute value) between W boson (truth level) and Z: [0, 5]
      CALL XHINIT(Ncout,0D0,5D0,Nbin,'Distribution Dy(W,Z)')
      Ncout = Ncout + 1
      
      ! transverse momentum of the WZ system: [0, 200 GeV]
      CALL XHINIT(Ncout,0D0,200D0,Nbin,'Distribution pT_WZ')
      Ncout = Ncout + 1
      
      ! R_21 ratio = pT_V2/pT_V1 (where V1 is the boson with hardest pT): [0, 1]
      CALL XHINIT(Ncout,0D0,1D0,Nbin,'Distribution pT_V2/pT_V1')
      Ncout = Ncout + 1
      
      ! (pT_e+pT_mup+pT_mum)/3
      CALL XHINIT(Ncout,0D0,500D0,Nbin,'Distribution average of pTs(lep) ')
      !Ncout = Ncout + 1
      End subroutine
