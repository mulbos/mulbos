subroutine KinDistSetup
use global_vars
use global_dists
implicit none
integer i, j, offset_3D
character(len=1024) :: format_string

!!![name of distribution files
format_string = "(A4,I3)" ! 4 letters and 3-digit integer as below
do i=1, N_dists_2D
 j = N_dists_2D + i
 ! 2D plots:
 write (filem(i),format_string) "dist", 400 + i
 ! 3D plots (600 because the limit of N_dists_2D is 200):
 offset_3D = max(400 + N_dists_2D, 600)
 write (filem(j),format_string) "dist", offset_3D + i
!        print *, trim(filem(i))
enddo
!!!]name of distribution files

end subroutine
