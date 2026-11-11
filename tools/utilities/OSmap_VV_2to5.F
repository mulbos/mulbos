      subroutine OSmap_VV_2to5(fac)
      use global_vars
      use global_mom_xpdf

      implicit none
      integer i, j
      double precision fac
      double precision qX(0:3), qX_inv(0:3), p_Xrest(0:3,proc_nlegs)
      double precision q_Xrest(0:3,proc_nlegs_os)
      double precision k_Xrest(0:3,proc_nlegs_os)      
      double precision Vext(0:3,proc_nlegs)
      double precision sext(0:3,proc_nlegs_os)
      double precision rm2(2)
      
      rm2(1) = mass_VV(1)**2        ! mass squared input
      rm2(2) = mass_VV(2)**2        ! mass squared input
      
      call translate_P2V(Vext(0,1),proc_nlegs) ! get Vext()

***[DPA approximation with on-shell projection
***[ boost to the X (= V1 + V2) rest frame
      do i = 0, 3
       qX(i) = Vext(i,1) + Vext(i,2) - Vext(i,3)
      enddo
      qX_inv(0) = qX(0)
      do i = 1, 3
      qX_inv(i) = -qX(i)
      enddo
      do i = 1, proc_nlegs
      call BOOSTN(Vext(0,i),qX_inv(0),p_Xrest(0,i))
      enddo
***] boost to the X (= V1 + V2) rest frame

! initial-state particles and the final jet
       do i = 1, 3
       do j = 0, 3
       q_Xrest(j,i) = p_Xrest(j,i)
       enddo
       enddo
! Lepton momenta: p_Xrest(0,4:7)
! the ordering: q1 q2 > g3 V4 V5 > l6 l7 l8 l9:
      call OSProject(osmap_case,p_Xrest(0,4),rm2(1),q_Xrest(0,4))
      
!!![re-ordering of particles
      if (lNLO_FIN_RAD) then
!!![ for the radiation process
!C the ordering: q1 q2 > V3 V4 g5 > g5 l6 l7 l8 l9, 
!C note, this ordering has to match the calculation of 
!C N+1 amplitude by FormCalc. 
      ! q_Xrest(0:3,proc_nlegs_os)  
      do i = 1, 2
       k_Xrest(:,i) = q_Xrest(:,i)
      enddo
      
      do i = 3, 4
       k_Xrest(:,i) = q_Xrest(:,i+1)
      enddo
      
      k_Xrest(:,5) = q_Xrest(:,3) ! final jet

      do i = 6, 9
       k_Xrest(:,i) = q_Xrest(:,i)
      enddo
!!!] for the radiation process
      else
       k_Xrest = q_Xrest
      endif
!!!]re-ordering of particles

!!![store os momenta in the VV CMS
! oext_vv(nlegs_os_max,4)! on-shell momenta in the VV CMS
      do i = 1, proc_nlegs_os
      do j = 0, 3
       oext_vv(i,j+1) = k_Xrest(j,i)
      enddo
      enddo
!!!]store os momenta in the VV CMS

!!![store os momenta in the partonic CMS
      do i = 1, proc_nlegs_os
! boost back from VV CMS to partonic CMS
       call BOOSTN(k_Xrest(0,i),qX(0),sext(0,i))
      enddo
! oext(nlegs_os_max,4)! on-shell momenta in the partonic CMS
      do i = 1, proc_nlegs_os
      do j = 0, 3
       oext(i,j+1) = sext(j,i)
      enddo
      enddo
!!!]store os momenta in the partonic CMS

!!![store os momenta in the Lab frame
      do i = 1, proc_nlegs_os
! boost from partonic CMS to the LAB frame
      call boost0(rgam_boost, -rbetgam_boost,
     &  sext(0,i), sext(1,i), sext(2,i), sext(3,i))
      enddo
! oext_pp(nlegs_os_max,4)! on-shell momenta in the pronton-proton CMS (Lab frame)
      do i = 1, proc_nlegs_os
      do j = 0, 3
       oext_pp(i,j+1) = sext(j,i)
      enddo
      enddo
!!!]store os momenta in the Lab frame
        return
***]DPA approximation with on-shell projection
	end subroutine
