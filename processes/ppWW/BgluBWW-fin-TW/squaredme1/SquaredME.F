	subroutine P1SquaredME(result, qext, flags)
        use global_vars
	implicit none
	double precision result(*)
	integer nlegs, flags
        parameter (nlegs = 6)
        double precision qext(nlegs,4)

	result(1) = 0
	result(2) = 0
	if (flags.eq.0) then
	 return
	endif

	end

