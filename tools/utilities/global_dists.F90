module global_dists
integer, parameter :: Nbin = 50
! this has to match with bshbok.f in bases
! 200 files for normal distributions, 
! 200 files for scatter plots
integer, parameter :: NHISMAX = 400
character(len=200) :: filem(NHISMAX)
integer, parameter :: N_dists_2D = 200 ! maximal value
integer, parameter :: N_dists_3D = 200 ! scatter plots (maximal value)
        
! histogram ranges:
double precision dis2Dmin(N_dists_2D),dis2Dmax(N_dists_2D)
end module
