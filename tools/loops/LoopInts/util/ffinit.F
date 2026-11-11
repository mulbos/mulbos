#include "defs.h"

! user-interface routines:

	subroutine loopints_init

	implicit none
	integer init
	save init
	data init /0/

	external DefaultDP, DefaultQP
	external SetPrecDP, SetPrecQP

***[authors
*	check whether there is anything to do
	if ( init .ne. 0 ) return
	init = 1
	print *,'===================================================='
	print *,' LoopInts: a package to evaluate one-loop integrals'
	print *,' written by LE Duc Ninh'
	print *,' Note: some routines (cache; 3,4-point tensor reduction)'
	print *,' are taken from LoopTools [hep-ph/9807565]'
	print *,'===================================================='
****]authors

	call clearcacheDP
	call clearcacheQP

***[precision etc:
	call SetPrecDP
	call SetPrecQP
***]precision etc:

***[default parameters
	call DefaultDP
	call DefaultQP
***]default parameters
	end
************************************************************************
        subroutine loopints_info
        implicit none
#include "ff.h"

#ifndef NO_QUAD
        call ffexiQP
#else
        print*,'LoopInts information:'
#ifdef NO_CACHE
        print*,'Caching is off.'
#else
        print*,'Caching is on.'
#endif
        print*,'UV regulator: delta  =',delta
        print*,'IR regulator: lambda =',lambda
        print*,'DR regulator: mudim  =',mudim
        print*,'xrat switch : xrat   =',xrat
        print*,'3pt maximal rank =',rankmax3
        print*,'4pt maximal rank =',rankmax4
        print*,'5pt maximal rank =',rankmax5
#endif
        call ffexiDP

        end
************************************************************************
	subroutine clearcache
	implicit none
	call clearcacheDP
	call clearcacheQP
	end

	subroutine markcache
	implicit none
	call markcacheDP
	call markcacheQP
	end

	subroutine restorecache
	implicit none
	call restorecacheDP
	call restorecacheQP
	end

	subroutine setmudim(mudim_)
	implicit none
	RVAR mudim_
	RVAR16 x
	x = mudim_
	call setmudimDP(mudim_)
	call setmudimQP(x)
	end
	
	subroutine setmudim_uv(mudim_)
	implicit none
	RVAR mudim_
	RVAR16 x
	x = mudim_
	call setmudimDP(mudim_)
	call setmudimQP(x)
	end

	subroutine setmudim_ir(mudim_)
	implicit none
	RVAR mudim_
	RVAR16 x
	x = mudim_
	call setmudimDP(mudim_)
	call setmudimQP(x)
	end

	subroutine setdelta(delta_)
	implicit none
	RVAR delta_
	RVAR16 x
	x = delta_
	call setdeltaDP(delta_)
	call setdeltaQP(x)
	end

	subroutine setlambda(lambda_)
	implicit none
	RVAR lambda_
	RVAR16 x
	x = lambda_
	call setlambdaDP(lambda_)
	call setlambdaQP(x)
	end

	subroutine setstab0(stab0_)
	implicit none
	RVAR stab0_
	RVAR16 x
	if ( stab0_*(1 - stab0_) .gt. 0 ) then 
	x = stab0_
	call setstab0DP(stab0_)
	call setstab0QP(x*x)
	endif
	end

	subroutine setsta03(stab0_)
	implicit none
	RVAR stab0_
	RVAR16 x
	if ( stab0_*(1 - stab0_) .gt. 0 ) then 
	x = stab0_
	call setstab0DP(stab0_)
	call setstab0QP(x*x)
	endif
	end

	subroutine setsta04(stab0_)
	implicit none
	RVAR stab0_
	RVAR16 x
	if ( stab0_*(1 - stab0_) .gt. 0 ) then 
	x = stab0_
	call setstab0DP(stab0_)
	call setstab0QP(x*x)
	endif
	end

	subroutine setstab2(stab2_)
	implicit none
	RVAR stab2_
	RVAR16 x
	x = stab2_
	call setstab2DP(stab2_)
	call setstab2QP(x)
	end

	subroutine setstab3(stab3_)
	implicit none
	RVAR stab3_
	RVAR16 x
	x = stab3_
	call setstab3DP(stab3_)
	call setstab3QP(x)
	end

	subroutine setstab4(stab4_)
	implicit none
	RVAR stab4_
	RVAR16 x
	x = stab4_
	call setstab4DP(stab4_)
	call setstab4QP(x)
	end

	subroutine setstab5(stab5_)
	implicit none
	RVAR stab5_
	RVAR16 x
	x = stab5_
	call setstab5DP(stab5_)
	call setstab5QP(x)
	end

	subroutine setstab6(stab6_)
	implicit none
	RVAR stab6_
	RVAR16 x
	x = stab6_
	call setstab6DP(stab6_)
	call setstab6QP(x)
	end

	subroutine setdebugkey(debugkey_)
	implicit none
	integer debugkey_
	call setdebugkeyDP(debugkey_)
	call setdebugkeyQP(debugkey_)
	end

	subroutine setversionkey(versionkey_)
	implicit none
	integer versionkey_
	call setversionkeyDP(versionkey_)
	call setversionkeyQP(versionkey_)
	end

	subroutine setmaxdev(maxdev_)
	implicit none
	RVAR maxdev_
	RVAR16 x
	x = maxdev_
	call setmaxdevDP(maxdev_)
	call setmaxdevQP(x)
	end

        subroutine setrunmode(key_)
        implicit none
        integer key_
        call setrunmodeDP(key_)
        call setrunmodeQP(key_)
        end

	subroutine setrankmax(rank3, rank4, rank5)
	implicit none
	integer rank3, rank4, rank5
! save it to the DP common block
	call setrankmaxDP(rank3, rank4, rank5)
! save it to the QP common block
	call setrankmaxQP(rank3, rank4, rank5)
	end        
        
        
        
        
