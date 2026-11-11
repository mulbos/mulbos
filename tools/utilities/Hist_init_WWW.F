      subroutine hist_init_WWW()
	use global_dists
	implicit none
      double precision pi
      parameter (pi = 3.1415926535897932384626433832795029D0)
      integer Ncout

!      print*,"hello from hist_init_WWW"

! particle ordering: p(1) p(2) > e(3) nu_e(4) mu+(5) mu-(6) 
      Ncout = 1

      ! transverse momentum of electron: [0, 500 GeV]
      CALL XHINIT(Ncout,0D0,500D0,Nbin,'Distribution pT(3)')
      Ncout = Ncout + 1                  

      ! transverse momentum of muon: [0, 500 GeV]
      CALL XHINIT(Ncout,0D0,500D0,Nbin,'Distribution pT(5)')
      Ncout = Ncout + 1                  

      ! transverse momentum of tauon: [0, 500 GeV]
      CALL XHINIT(Ncout,0D0,500D0,Nbin,'Distribution pT(7)')
      Ncout = Ncout + 1                  

      ! cosine of the theta_www angle of e (W frame) and W boson (WWW rest frame), truth level
      CALL XHINIT(Ncout,-1D0,1D0,Nbin,'Distribution cos(theta_www_W1_e)')
      Ncout = Ncout + 1

      ! cosine of the theta_www angle of mu (W frame) and W boson (WWW rest frame), truth level
      CALL XHINIT(Ncout,-1D0,1D0,Nbin,'Distribution cos(theta_www_W2_mu)')
      Ncout = Ncout + 1

      ! cosine of the theta_www angle of tau (W frame) and W boson (WWW rest frame), truth level
      CALL XHINIT(Ncout,-1D0,1D0,Nbin,'Distribution cos(theta_www_W3_tau)')
      Ncout = Ncout + 1


      End subroutine hist_init_WWW
