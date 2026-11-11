* A0s.F
* the scalar one-point function

#include "defs.h"

	CVAR function CYPA0s(m)
	implicit none
	CVAR m

#include "lt.h"

	if( lambda .lt. 0 ) then
	  CYPA0s = 0
	  return
	endif

	if( abs(m) .le. 0 ) then
	  CYPA0s = 0
	else
	  CYPA0s = m*(1 - CYPLOG(m/mudim, -1D0) + delta)
	endif
	end

****************************************************************
	CVAR function CYPA00s(m)
	implicit none
	CVAR m

#include "lt.h"

	CVAR CYPA0s
	external CYPA0s

        if( lambda .lt. 0 ) then
          CYPA00s = 0
          return
        endif

	CYPA00s = .25D0*m*(CYPA0s(m) + .5D0*m)
	end
