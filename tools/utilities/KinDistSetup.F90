subroutine KinDistSetup
use global_vars
use global_dists
implicit none
integer i, j
character(len=1024) :: format_string

!!![name of distribution files
format_string = "(A4,I3)" ! 4 letters and 3-digit integer as below
do i=1, N_dists_2D
 j = N_dists_2D + i
 write (filem(i),format_string) "dist", 400+i
 write (filem(j),format_string) "dist", 500+i
!        print *, trim(filem(i))
enddo
!!!]name of distribution files

end subroutine
