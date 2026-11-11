	subroutine Do_VecSet_DPA_FC9()
	use global_vars
	use global_mom_xpdf
! Input: 
! Output:
	implicit none
	integer i
        double precision mext(proc_nlegs_os), tmp2
        double precision p_os(nlegs_os_max,4)! on-shell momenta

        mext(:) = 0D0
        if (lNLO_FIN_IND) then
          mext(4) = mass_VV(1)
          mext(5) = mass_VV(2)
        else
          mext(3) = mass_VV(1)
          mext(4) = mass_VV(2)
        endif
        
        if (frame_polar.eq.1) then ! VV cms frame
         p_os = oext_vv
        else if (frame_polar.eq.2) then ! lab frame
         p_os = oext_pp
        else if (frame_polar.eq.3) then ! partonic frame
         p_os = oext
        endif

        do i = 1, proc_nlegs_os
        tmp2 = sqrt(p_os(i,2)**2 + p_os(i,3)**2 + p_os(i,4)**2)
        call VecSet_FC9(i, mext(i), tmp2, p_os(i,2)/tmp2, p_os(i,3)/tmp2, p_os(i,4)/tmp2)
        enddo
        
        return
	end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine Do_VecSet_FULL_FC9()
	use global_vars
	use global_mom_xpdf
! Input: 
! Output:
	implicit none
	integer i
        double precision tmp2

        do i = 1, proc_nlegs
        tmp2 = sqrt(pext(i,2)**2 + pext(i,3)**2 + pext(i,4)**2)
        call VecSet_FC9(i, 0D0, tmp2, pext(i,2)/tmp2, pext(i,3)/tmp2, pext(i,4)/tmp2)
        enddo
        return

	end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine Do_VecSet_DPA_TW_FC9()
	use global_vars
	use global_mom_xpdf
! Input: 
! Output:
	implicit none
	integer i, j
        double precision mext(proc_nlegs_os), tmp2
        double precision p_os(nlegs_os_max,4)! on-shell momenta
        double precision pT(4), pW_rem(4), pW_fromT(4)

        mext(:) = 0D0
        mext(4) = mass_VV(1)
        mext(5) = mass_VV(2)
        
        if (frame_polar.eq.1) then ! VV cms frame
         p_os = oext_vv
        else if (frame_polar.eq.2) then ! lab frame
         p_os = oext_pp
        else if (frame_polar.eq.3) then ! partonic frame
         p_os = oext
        endif

!        do i = 1, proc_nlegs_os
        do i = 1, 2
        tmp2 = sqrt(p_os(i,2)**2 + p_os(i,3)**2 + p_os(i,4)**2)
        call VecSet_FC9(i, mext(i), tmp2, p_os(i,2)/tmp2, p_os(i,3)/tmp2, p_os(i,4)/tmp2)
        enddo
        ! top quark
      if (leg(1).eq.-5) then ! b_bar gluon process
       pT(:) = p_os(3,:) + p_os(5,:) ! Tbar
       pW_rem(:) = p_os(4,:) ! Wp
       pW_fromT(:) = p_os(5,:) ! Wm
      else ! b gluon process
       pT(:) = p_os(3,:) + p_os(4,:) ! Top
       pW_rem(:) = p_os(5,:) ! Wm
       pW_fromT(:) = p_os(4,:) ! Wp
      endif
      tmp2 = sqrt(pT(2)**2 + pT(3)**2 + pT(4)**2)
      call VecSet_FC9(3, mass_Top, tmp2, pT(2)/tmp2, pT(3)/tmp2, pT(4)/tmp2) ! top
      tmp2 = sqrt(pW_rem(2)**2 + pW_rem(3)**2 + pW_rem(4)**2)
      call VecSet_FC9(4, mext(4), tmp2, pW_rem(2)/tmp2, pW_rem(3)/tmp2, pW_rem(4)/tmp2) ! W_rem
      tmp2 = sqrt(pW_fromT(2)**2 + pW_fromT(3)**2 + pW_fromT(4)**2)
      call VecSet_FC9(5, mext(5), tmp2, pW_fromT(2)/tmp2, pW_fromT(3)/tmp2, pW_fromT(4)/tmp2) ! W from Top
      tmp2 = sqrt(p_os(3,2)**2 + p_os(3,3)**2 + p_os(3,4)**2)
      call VecSet_FC9(6, 0d0, tmp2, p_os(3,2)/tmp2, p_os(3,3)/tmp2, p_os(3,4)/tmp2) ! bottom
        do i = 6, 9
        j = i + 1
        tmp2 = sqrt(p_os(i,2)**2 + p_os(i,3)**2 + p_os(i,4)**2)
        call VecSet_FC9(j, 0d0, tmp2, p_os(i,2)/tmp2, p_os(i,3)/tmp2, p_os(i,4)/tmp2)
        enddo
        return
	end subroutine
