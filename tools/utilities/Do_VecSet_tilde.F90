	subroutine Do_VecSet_tilde()
	use global_vars
	use global_mom_xpdf
! Input: 
! Output:
	implicit none
	integer i, j
        double precision temp
        double precision p_os(nlegs_os_max,4) ! on-shell momenta
        double precision p, ex,ey,ez, psq
        
        ! oext_tilde(nlegs_os_max,4), p_os(nlegs_os_max,4)
        if (frame_polar.eq.1) then ! VV cms frame
         p_os = oext_tilde_vv
        else if (frame_polar.eq.2) then ! lab frame
         p_os = oext_tilde_pp
        else if (frame_polar.eq.3) then ! parent partonic frame
         p_os = oext_tilde
        endif
        
!***[for the q1 q2 > V3 V4 > l5 l6 l7 l8
! run over 2 initial particles:
        do i = 1, 2
        psq = p_os(i,2)**2+p_os(i,3)**2+p_os(i,4)**2
	p = sqrt(psq)
	if (p .eq. 0) then
	ex = 0
	ey = 0
	ez = 1
	else
	ex = p_os(i,2)/p
	ey = p_os(i,3)/p
	ez = p_os(i,4)/p
	endif
	
	call VecSet(i, 0d0, p, ex,ey,ez) ! massless initial quarks
	enddo

! run over resonances:
        do j = 1, nlegs_res 
        i = j + 2
        psq = p_os(i,2)**2+p_os(i,3)**2+p_os(i,4)**2
	p = sqrt(psq)
	if (p .eq. 0) then
	ex = 0
	ey = 0
	ez = 1
	else
	ex = p_os(i,2)/p
	ey = p_os(i,3)/p
	ez = p_os(i,4)/p
	endif
	
	! calculate the mass from the momentum:
	!temp = p_os(i,1)**2-psq
	!if (abs(temp) .lt. 1D-8*psq) then
	 !temp = 0
	!else
	 !temp = sqrt(temp)
	!endif   
	temp = mass_VV(j)
	call VecSet(i, temp, p, ex,ey,ez)
	enddo

! run over final state leptons (Nlepton):
        do j = 1, Nlepton
        i = 2 + nlegs_res + j
        psq = p_os(i,2)**2+p_os(i,3)**2+p_os(i,4)**2
	p = sqrt(psq)
	if (p .eq. 0) then
	ex = 0
	ey = 0
	ez = 1
	else
	ex = p_os(i,2)/p
	ey = p_os(i,3)/p
	ez = p_os(i,4)/p
	endif
	
	call VecSet(i, 0d0, p, ex,ey,ez) ! massless leptons
	enddo
!***]for the q1 q2 > V3 V4 > l5 l6 l7 l8
        
        return
	end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
