	program main
	!use ifport ! to open a directory for distributions (ifx only)
	use global_vars
	use global_cuts
	implicit none

	integer fail
	character*512 file_input, file_cuts
	character*200 path_home
        !logical lopen_dir

        path_home = "./"
c        file_input = path_home(1:lnblnk(path_home))//'input.dat'
        file_input = trim(path_home) // 'input.dat'
c        file_cuts = path_home(1:lnblnk(path_home))//'cuts.dat'
        file_cuts = trim(path_home) // 'cuts.dat'
        call Read_Input(file_input)
        print*,"Read_Input done"
        call Read_Cuts(file_cuts)
        print*,"Read_Cuts done"
c        dir_save = path_home(1:lnblnk(path_home))//dir_name
        dir_save = trim(path_home) // dir_name
c        dir_save = dir_save(1:lnblnk(dir_save))//'/'
        dir_save = trim(dir_save) // '/'
        call set_run_params()
        print*,"set_run_params done"
        call get_sm_para_ext() ! calculate the common block: sm_para_ext

	call ProcessIni(fail)
	if( fail .ne. 0 ) goto 999
	
        !lopen_dir=makedirqq(dir_save) ! intel ifx
        call makedir_mb(dir_save)
        
        call IntegratedCS

	goto 1

999	continue
	print*, "Error: cannot initiate the process"
	stop

1	continue
	end
*******************************************************************************************
      subroutine makedir_mb(dir)
      implicit none
      character(len=*), intent(in) :: dir
      integer :: ierr
c      ierr = system('mkdir -p ' // trim(dir))
      ierr = 0
      call execute_command_line('mkdir -p ' // trim(dir), exitstat=ierr)
      if (ierr == 0) then
        print *, 'Directory created successfully, dir =',trim(dir)
      else
        print *, 'Failed to create directory, dir =',trim(dir)
      end if
      end subroutine makedir_mb
*******************************************************************************************	
	subroutine linking_helper
	! This routine calls VecSet() and SquaredME() so that 
	! xsection.F and Integrand.F in utilities/ know them. 
	! This routine is never called, and hence has no numerical effects on the results.
	use global_vars
	implicit none
        double precision temp(2), q_LO(nlegs_max,4), MVint2(nlegs_res_max)
     	! calling vecset(): needed for utilities.a to find vecset()
	call VecSet(1, 0D0, 1D0, 0D0, 0D0, 1D0)
	! calling squaredme(): needed for utilities.a to find squaredme()
	q_LO = 0D0
	!q_NLO = 0D0
	MVint2 = 0D0
        call SquaredME(temp(1), q_LO(1,1), 0) ! LO or NLO-fin
        call P1SquaredME(temp(1), q_LO(1,1), 0) ! photon-photon induced contribution
        call P0SquaredME(temp(1), MVint2(1), 0) ! dipole subtraction for (n+1)
       end subroutine
