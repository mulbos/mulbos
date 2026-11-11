      subroutine OSmap_VV_2to6(fac)
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
      
      ! pext(nlegs_max,4) ! exact momenta in the partonic CMS
      do i = 1, proc_nlegs
       Vext(:,i) = pext(i,:)
      enddo

***[DPA approximation with on-shell projection
***[ boost to the X (= V1 + V2) rest frame
      do i = 0, 3
       qX(i) = Vext(i,3) + Vext(i,4) + Vext(i,5) + Vext(i,6)
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
       do i = 1, 2
       q_Xrest(:,i) = p_Xrest(:,i)
       enddo
       q_Xrest(:,9) = p_Xrest(:,7)
       q_Xrest(:,10)= p_Xrest(:,8)
       
! Lepton momenta: p_Xrest(0,4:7)
! the ordering: l3 - l6 > W3W4l5l6l7l8
      call OSProject(osmap_case,p_Xrest(0,3),rm2(1),q_Xrest(0,3))
      
!!![re-ordering of particles
!C the ordering: q1 q2 > V3 V4 b5 b6 > V3 V4 l5 l6 l7 l8 b9 b10  
      ! q_Xrest(0:3,proc_nlegs_os)
      
       k_Xrest = q_Xrest

!!!]re-ordering of particles

!!![store os momenta in the VV CMS
! oext_vv(nlegs_os_max,4)! on-shell momenta in the VV CMS
      do i = 1, proc_nlegs_os
       oext_vv(i,:) = k_Xrest(:,i)
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
