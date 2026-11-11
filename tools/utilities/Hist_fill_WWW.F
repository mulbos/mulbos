	subroutine Fill_WWW(qext, FUNC)
! Note: qext are in LAB frame with order q1 q2 > l3 l4 ... jet
	use global_vars
	implicit none
	integer i, j, Ncout, ipart

      double precision pi
      parameter (pi = 3.1415926535897932384626433832795029D0)

      double precision qext(proc_nlegs,4)
      double precision rext(4,proc_nlegs),qV(4,3)
      double precision FUNC, xdis
      double precision temp
      double precision dot_space, dot, dot_p2, inv_mass, inv_mass_p
      double precision pt(proc_nlegs), py(proc_nlegs), phi(proc_nlegs)
      double precision q_clep(8,Nlep_charged) ! charged lep mom and particle id number
      double precision c_theta_VVV_V1(2),c_theta_VVV_V2(2),c_theta_VVV_V3(2)

! qext(proc_nlegs,4): p_1 p_2 > e_3 neu_4 mu_5 neu_6 (order of particles)

!      print*,"hello from Fill_WWW"

        do i = 1, 4
         do j = 1, proc_nlegs   
         rext(i,j) = qext(j,i)
         enddo
        enddo

***[all final-sate particles
      do i = 3, proc_nlegs
       call lego( rext(1,i), pt(i), py(i), phi(i) )
      enddo
***]all final-sate particles

      call Get_charged_leptons(proc_nlegs,rext(1,1),pt(1),py(1),phi(1),q_clep(1,1))

***[the resonances
      qV(:,1) = rext(:,3)+rext(:,4)
      qV(:,2) = rext(:,5)+rext(:,6)
      qV(:,3) = rext(:,7)+rext(:,8)
***]the resonances

***{ Fill histograms and scatter plot
      Ncout = 1
! transverse momentum of electron3: [0, 500 GeV]
! CALL XHINIT(Ncout,0D0,500D0,Nbin,'Distribution pT(e)')  
      do i = 1, Nlep_charged
       CALL XHFILL(Ncout,q_clep(5,i),FUNC)
       Ncout = Ncout + 1
      enddo
       
      ! cosine of the theta_www angle of e (W frame) and W boson (WWW rest frame), truth level
      ! CALL XHINIT(Ncout,-1D0,1D0,Nbin,'Distribution cos(theta_www_W1_e)')
       call theta_vvv(1,rext(1,3),qV(1,1),c_theta_VVV_V1(1))
       CALL XHFILL(Ncout,c_theta_VVV_V1(1),FUNC)
       Ncout = Ncout + 1    

      ! cosine of the theta_www angle of mu (W frame) and W boson (WWW rest frame), truth level
       call theta_vvv(2,rext(1,3),qV(1,1),c_theta_VVV_V2(1))
       CALL XHFILL(Ncout,c_theta_VVV_V2(1),FUNC)
       Ncout = Ncout + 1    

      ! cosine of the theta_www angle of tau (W frame) and W boson (WWW rest frame), truth level
       call theta_vvv(3,rext(1,3),qV(1,1),c_theta_VVV_V3(1))
       CALL XHFILL(Ncout,c_theta_VVV_V3(1),FUNC)
       Ncout = Ncout + 1    


***} Fill histograms and scatter plot
      end subroutine Fill_WWW
