      subroutine hist_init_ZZ()
	use global_dists
	implicit none
      double precision pi
      parameter (pi = 3.1415926535897932384626433832795029D0)
      integer Ncout

C ATLAS-2022 wish list
! particle ordering: p(1) p(2) > e+(3) e-(4) mu+(5) mu-(6) 
      Ncout = 1

! invariant mass
      
      ! invariant mass of the four-lepton system, truth level: [0, 2000 GeV]
      CALL XHINIT(Ncout,0D0,2D3,50,'Distribution m(3,4,5,6)')
      Ncout = Ncout + 1
      
      CALL XHINIT(Ncout,0D0,1D3,50,'Distribution m(3,4,5)')
      Ncout = Ncout + 1
      
      CALL XHINIT(Ncout,0D0,8D2,40,'Distribution m(3,5)')
      Ncout = Ncout + 1
      
      CALL XHINIT(Ncout,81D0,101D0,40,'Distribution m(3,4)')
      Ncout = Ncout + 1

! transverse momentum

      ! transverse momentum of the 4L system: [0, 300 GeV]
      CALL XHINIT(Ncout,0D0,300D0,30,'Distribution pT(3,4,5,6)')
      Ncout = Ncout + 1

      ! transverse momentum of the (3,5) system: [0, 700 GeV]
      CALL XHINIT(Ncout,0D0,600D0,40,'Distribution pT(3,5)')
      Ncout = Ncout + 1

      ! transverse momentum of the Z boson: [0, 700 GeV]
      CALL XHINIT(Ncout,0D0,600D0,40,'Distribution pT(3,4)')
      Ncout = Ncout + 1      
      
      ! transverse momentum of electron: [0, 500 GeV]
      CALL XHINIT(Ncout,0D0,400D0,40,'Distribution pT(3)')
      Ncout = Ncout + 1                  
      
      ! transverse momentum of muon minus: [0, 500 GeV]
      CALL XHINIT(Ncout,0D0,400D0,40,'Distribution pT(6)')
      Ncout = Ncout + 1
      
! rapidity      
      
      ! rapidity of the electron: [−2.5, 2.5]
      CALL XHINIT(Ncout,-2.5D0,2.5D0,20,'Distribution y(3)')
      Ncout = Ncout + 1
      
      ! rapidity of the mu_minus: [−2.5, 2.5]
      CALL XHINIT(Ncout,-2.5D0,2.5D0,20,'Distribution y(6)')
      Ncout = Ncout + 1
      
      ! rapidity of the Z boson: [−2.5, 2.5]
      CALL XHINIT(Ncout,-3D0,3D0,24,'Distribution y(3,4)')
      Ncout = Ncout + 1
      
! rapidity separation      
      
      ! rapidity separation (in absolute value) between e+ and e-: [0, 5]
      CALL XHINIT(Ncout,0D0,5D0,20,'Distribution |Dy(3,4)|')
      Ncout = Ncout + 1
      
      ! rapidity separation (in absolute value) between e+ and mu+: [0, 5]
      CALL XHINIT(Ncout,0D0,5D0,20,'Distribution |Dy(3,5)|')
      Ncout = Ncout + 1
      
      ! rapidity separation (in absolute value) between e+ and the (5,6) boson (muons system): [0, 5] 
      CALL XHINIT(Ncout,0D0,5D0,20,'Distribution Dy(e,Zmu)')
      Ncout = Ncout + 1
      
! azimuthal separation      
      ! azimuthal separation between e+ and e-
      CALL XHINIT(Ncout,0D0,180D0,30,'Distribution Dphi(3,4)')
      Ncout = Ncout + 1
      
      ! azimuthal separation between e+ and mu+
      CALL XHINIT(Ncout,0D0,180D0,30,'Distribution Dphi(3,5)')
      Ncout = Ncout + 1

! R distance
      ! R separation between e+ and e-
      CALL XHINIT(Ncout,0D0,6D0,30,'Distribution DR(3,4)')
      Ncout = Ncout + 1
      
      ! R separation between e+ and mu+
      CALL XHINIT(Ncout,0D0,6D0,30,'Distribution DR(3,5)')
      Ncout = Ncout + 1

! cosine of the polar decay angle

      ! cosine of the theta_wz angle of e (W frame) and W boson (WZ rest frame), truth level
      CALL XHINIT(Ncout,-1D0,1D0,40,'Distribution cos(theta_wz_W_e)')
      Ncout = Ncout + 1
      
      ! cosine of the theta_wz angle of mu_minus (Z frame) and Z boson (WZ rest frame), truth level
      CALL XHINIT(Ncout,-1D0,1D0,40,'Distribution cos(theta_wz_Z_mu6)')
      Ncout = Ncout + 1      

! cosine of the scattering angle
! scattering angle (in absolute value) of the single boson in the WZ CM frame computed w.r.t.
!the positive Z-axis direction of the LAB, truth level: [0, 1]
      CALL XHINIT(Ncout,0D0,1D0,40,'Distribution |cos(W_WZ)|')
      Ncout = Ncout + 1
      
! cossine of angular separation

      ! cosine of the angle between e+ and e-:
      CALL XHINIT(Ncout,-1D0,1D0,40,'Distribution cos(theta_34)')
      Ncout = Ncout + 1
      
      ! cosine of the angle between e+ and mu+:
      CALL XHINIT(Ncout,-1D0,1D0,40,'Distribution cos(theta_35)')
      !Ncout = Ncout + 1

      End subroutine
      
      
      
