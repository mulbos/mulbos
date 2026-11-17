	subroutine scaleFR(lSub)
	use global_vars
	use global_mom_xpdf
! Input: lSub
! Output: factorization and renormalization scales calculated
	implicit none
	logical lSub
	integer i
	integer scale_ID ! for adding a new dynamical scale
        double precision tmp, scale0(2:n_scales_max) ! 0, 1 are reserved for fix scales
        double precision p(4,nlegs_max), q(4)
        double precision inv_mass_p, inv_mass_p2, inv_mass, inv_mass2

        p(:,:) = 0d0
        do i = 1, proc_nlegs
        if (lSub) then ! dipole subtraction
! pext_tilde(nlegs_max,4) ! exact momenta in the parent partonic CMS (for cut, hist), same particle order as LO
         p(:,i) = pext_tilde(i,:)
        else
! vext_jend(5,nlegs_max) ! exact momenta in the partonic CMS (for cut, hist)
         p(:,i) = vext_jend(1:4,i)
        endif
        enddo
        
!!![default dynamical scale choice with ID = 2
        scale_ID = 2 ! default dynamical scale 
        if (scale_muf_case.eq.scale_ID.or.scale_mur_case.eq.scale_ID) then
          q(:) = 0d0
          do i = 3, Nlepton + 2
           q(:) = q(:) + p(:,i)
          enddo
! scale0 = m_inv_leptons/2:
        tmp = sqrt(q(2)**2+q(3)**2+q(4)**2)
        scale0(scale_ID) = 0.5D0*sqrt((q(1)-tmp)*(q(1)+tmp))
        endif
!!!]default dynamical scale choice with ID = 2

!!![new dynamical scale choice with ID > 2
        ! scale_ID = XXX ! choose an integer value > 2 and has not been used 
        if (scale_muf_case.gt.2.or.scale_mur_case.gt.2) then
        print*,"Error: this dynamical scale choice has not been implemented"
        print*,"scale_muf_case =",scale_muf_case
        print*,"scale_mur_case =",scale_mur_case
        print*,"check your input scale choice in input.dat.in"
        print*,"or implement the new scale here"
        stop
        endif
!!!]new dynamical scale choice with ID > 2

        ! whether to update scales:
        if (lscaleF_dyn) then
         scale = scale0(scale_muf_case)*scale_facF
        endif
        if (lscaleR_dyn) then
         scaleren = scale0(scale_mur_case)*scale_facR
        endif

        return
	end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

