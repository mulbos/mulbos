* ini.F

#include "defs.h"

	subroutine YPclearcache
	implicit none

#include "lt.h"

#ifndef NO_CACHE
	integer i

	do i = 1, ncaches
#ifndef QUADPREC
	  cacheptr(1,1,i) = 0
	  cacheptr(1,2,i) = 0
#else
	  cacheptrQP(1,i) = 0
	  cacheptrQP(2,i) = 0
#endif
	enddo

#else
        iBget = 0
        iCget = 0
        iDget = 0
        iEget = 0
        iFget = 0 
#endif
	end

************************************************************************

	subroutine YPmarkcache
	implicit none

#include "lt.h"

#ifndef NO_CACHE
	integer i

	do i = 1, ncaches
#ifndef QUADPREC
	  savedptr(1,i) = cacheptr(1,1,i)
	  savedptr(2,i) = cacheptr(1,2,i)
#else
	  savedptrQP(1,i) = cacheptrQP(1,i)
	  savedptrQP(2,i) = cacheptrQP(2,i)
#endif
	enddo

#else
        iBgetsv = iBget
        iCgetsv = iCget
        iDgetsv = iDget
        iEgetsv = iEget
        iFgetsv = iFget 	
#endif
	end

************************************************************************

	subroutine YPrestorecache
	implicit none

#include "lt.h"

#ifndef NO_CACHE
	integer i

	do i = 1, ncaches
#ifndef QUADPREC
	  cacheptr(1,1,i) = savedptr(1,i)
	  cacheptr(1,2,i) = savedptr(2,i)
#else
	  cacheptrQP(1,i) = savedptrQP(1,i)
	  cacheptrQP(2,i) = savedptrQP(2,i)
#endif
	enddo

#else
        iBget = iBgetsv
        iCget = iCgetsv
        iDget = iDgetsv
        iEget = iEgetsv
        iFget = iFgetsv
#endif
	end

************************************************************************

	subroutine YPsetmudim(mudim_)
	implicit none
	RVAR mudim_

#include "lt.h"

	mudim = mudim_
	call YPclearcache
	end

************************************************************************

	RVAR function YPgetmudim()
	implicit none

#include "lt.h"

	YPgetmudim = mudim
	end

************************************************************************

	subroutine YPsetdelta(delta_)
	implicit none
	RVAR delta_

#include "lt.h"

	delta = delta_
	call YPclearcache
	end

************************************************************************

	RVAR function YPgetdelta()
	implicit none

#include "lt.h"

	YPgetdelta = delta
	end

************************************************************************

	subroutine YPsetlambda(lambda_)
	implicit none
	RVAR lambda_

#include "lt.h"

	lambda = lambda_
	if ( lambda .lt. 0 ) xrat = 0
	call YPclearcache
	end

************************************************************************

	RVAR function YPgetlambda()
	implicit none

#include "lt.h"

	YPgetlambda = lambda
	end

************************************************************************

	subroutine YPsetstab0(stab0_)
	implicit none
	RVAR stab0_

#include "lt.h"
	if ( stab0_*(1 - stab0_) .gt. 0 ) then 
	stab0 = stab0_
	call YPclearcache
	endif
	end

************************************************************************

	subroutine YPsetstab2(stab2_)
	implicit none
	RVAR stab2_

#include "lt.h"

	stab2 = stab2_
	call YPclearcache
	end

************************************************************************

	subroutine YPsetstab3(stab3_)
	implicit none
	RVAR stab3_

#include "lt.h"

	stab3 = stab3_
	call YPclearcache
	end

************************************************************************

	subroutine YPsetstab4(stab4_)
	implicit none
	RVAR stab4_

#include "lt.h"

	stab4 = stab4_
	call YPclearcache
	end

************************************************************************

	subroutine YPsetstab5(stab5_)
	implicit none
	RVAR stab5_

#include "lt.h"

	stab5 = stab5_
	call YPclearcache
	end

************************************************************************

	subroutine YPsetstab6(stab6_)
	implicit none
	RVAR stab6_

#include "lt.h"

	stab6 = stab6_
	call YPclearcache
	end

************************************************************************

	subroutine YPsetdebugkey(debugkey_)
	implicit none
	integer debugkey_

#include "lt.h"

	debugkey = debugkey_
	end

************************************************************************

	subroutine YPsetversionkey(versionkey_)
	implicit none
	integer versionkey_

#include "lt.h"

	versionkey = versionkey_
	call YPclearcache
	end

************************************************************************

	subroutine YPsetmaxdev(maxdev_)
	implicit none
	RVAR maxdev_

#include "lt.h"

	maxdev = maxdev_
	end

************************************************************************

	subroutine YPsetrankmax(rank3, rank4, rank5)
	implicit none
	integer rank3, rank4, rank5

#include "lt.h"

	rankmax3 = rank3
	rankmax4 = rank4
	rankmax5 = rank5
	end

************************************************************************

        subroutine YPsetrunmode(key_)
        implicit none
        integer key_

#include "lt.h"

        runmode = key_
        end

************************************************************************

        subroutine YPffexi
        implicit none
#include "lt.h"

#ifdef QUADPREC
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

        print*,'precision para   : prec   =',prec
        print*,'small-mass apprx : stab0  =',stab0

#ifndef QUADPREC
cc      print*,'2-pt cutoff : stab2  =',stab2
        print*,'3-pt cutoff : stab3  =',stab3
        print*,'4-pt cutoff : stab4  =',stab4
        print*,'5-pt cutoff : stab5  =',stab5
        print*,'6-pt cutoff : stab6  =',stab6
#endif
        end

************************************************************************

	block data YPNameData
	implicit none
	integer i

	character*6 paraname(Pff,2:6)
	common /YPparanames/ paraname

	character*8 coeffname(Nff,2:6)
	common /YPcoeffnames/ coeffname

	data (paraname(i,2), i = 1, Pbb) /
     &    "p", "m1", "m2" /

	data (paraname(i,3), i = 1, Pcc) /
     &    "p1", "p2", "p1p2", "m1", "m2", "m3" /

	data (paraname(i,4), i = 1, Pdd) /
     &    "p1", "p2", "p3", "p4", "p1p2", "p2p3",
     &    "m1", "m2", "m3", "m4" /

	data (paraname(i,5), i = 1, Pee) /
     &    "p1", "p2", "p3", "p4", "p5",
     &    "p1p2", "p2p3", "p3p4", "p4p5", "p5p1",
     &    "m1", "m2", "m3", "m4", "m5" /

	data (paraname(i,6), i = 1, Pff) /
     &    "p1", "p2", "p3", "p4", "p5", "p6",
     &    "p1p2", "p2p3", "p3p4", "p4p5", "p5p6", "p6p1", 
     &    "p1p2p3", "p2p3p4", "p3p4p5", 
     &    "m1", "m2", "m3", "m4", "m5", "m6" /

	data (coeffname(i,2), i = 1, Nbb) /
     &  "bb0", "bb1", "bb00", "bb11", "bb001", "bb111",
     &  "dbb0", "dbb1", "dbb00", "dbb11","bb0000", "bb0011", "bb1111",
     &  "bb00001", "bb00111", "bb11111"/

	data (coeffname(i,3), i = 1, Ncc) /
     &    "cc0", "cc1", "cc2", "cc00", "cc11", "cc12", "cc22",
     &    "cc001", "cc002", "cc111", "cc112", "cc122", "cc222",
     &    "cc0000", "cc0011", "cc0012", "cc0022", "cc1111",
     &    "cc1112", "cc1122", "cc1222", "cc2222",
     &    "cc00001", "cc00002", "cc00111", "cc00112",
     &    "cc00122", "cc00222", "cc11111", "cc11112",
     &    "cc11122", "cc11222", "cc12222", "cc22222" /

	data (coeffname(i,4), i = 1, Ndd) /
     &    "dd0", "dd1", "dd2", "dd3", "dd00", "dd11", "dd12",
     &    "dd13", "dd22", "dd23", "dd33", "dd001", "dd002", "dd003",
     &    "dd111", "dd112", "dd113", "dd122", "dd123", "dd133",
     &    "dd222", "dd223", "dd233", "dd333", "dd0000", "dd0011",
     &    "dd0012", "dd0013", "dd0022", "dd0023", "dd0033", "dd1111",
     &    "dd1112", "dd1113", "dd1122", "dd1123", "dd1133", "dd1222",
     &    "dd1223", "dd1233", "dd1333", "dd2222", "dd2223", "dd2233",
     &    "dd2333", "dd3333", "dd00001", "dd00002", "dd00003",
     &    "dd00111", "dd00112", "dd00113", "dd00122", "dd00123",
     &    "dd00133", "dd00222", "dd00223", "dd00233", "dd00333",
     &    "dd11111", "dd11112", "dd11113", "dd11122", "dd11123",
     &    "dd11133", "dd11222", "dd11223", "dd11233", "dd11333",
     &    "dd12222", "dd12223", "dd12233", "dd12333", "dd13333",
     &    "dd22222", "dd22223", "dd22233", "dd22333", "dd23333",
     &    "dd33333" /

	data (coeffname(i,5), i = 1, Nee) /
     &    "ee0", "ee1", "ee2", "ee3", "ee4", "ee00", "ee11",
     &    "ee12", "ee13", "ee14", "ee22", "ee23", "ee24", "ee33",
     &    "ee34", "ee44", "ee001", "ee002", "ee003", "ee004",
     &    "ee111", "ee112", "ee113", "ee114", "ee122", "ee123",
     &    "ee124", "ee133", "ee134", "ee144", "ee222", "ee223",
     &    "ee224", "ee233", "ee234", "ee244", "ee333", "ee334",
     &    "ee344", "ee444", "ee0000", "ee0011", "ee0012", "ee0013",
     &    "ee0014", "ee0022", "ee0023", "ee0024", "ee0033", "ee0034",
     &    "ee0044", "ee1111", "ee1112", "ee1113", "ee1114", "ee1122",
     &    "ee1123", "ee1124", "ee1133", "ee1134", "ee1144", "ee1222",
     &    "ee1223", "ee1224", "ee1233", "ee1234", "ee1244", "ee1333",
     &    "ee1334", "ee1344", "ee1444", "ee2222", "ee2223", "ee2224",
     &    "ee2233", "ee2234", "ee2244", "ee2333", "ee2334", "ee2344",
     &    "ee2444", "ee3333", "ee3334", "ee3344", "ee3444", "ee4444", 
     &    "ee00001", "ee00002", "ee00003", "ee00004", "ee00111", 
     &    "ee00112", "ee00122", "ee00222", "ee00113", "ee00123", 
     &    "ee00223", "ee00133", "ee00233", "ee00333", "ee00114", 
     &    "ee00124", "ee00224", "ee00134", "ee00234", "ee00334", 
     &    "ee00144", "ee00244", "ee00344", "ee00444", "ee11111", 
     &    "ee11112", "ee11122", "ee11222", "ee12222", "ee22222", 
     &    "ee11113", "ee11123", "ee11223", "ee12223", "ee22223", 
     &    "ee11133", "ee11233", "ee12233", "ee22233", "ee11333", 
     &    "ee12333", "ee22333", "ee13333", "ee23333", "ee33333", 
     &    "ee11114", "ee11124", "ee11224", "ee12224", "ee22224", 
     &    "ee11134", "ee11234", "ee12234", "ee22234", "ee11334", 
     &    "ee12334", "ee22334", "ee13334", "ee23334", "ee33334", 
     &    "ee11144", "ee11244", "ee12244", "ee22244", "ee11344", 
     &    "ee12344", "ee22344", "ee13344", "ee23344", "ee33344", 
     &    "ee11444", "ee12444", "ee22444", "ee13444", "ee23444", 
     &    "ee33444", "ee14444", "ee24444", "ee34444", "ee44444" /

	data (coeffname(i,6), i = 1, Nff) /
     &    "ff0", "ff1", "ff2", "ff3", "ff4", "ff5", "ff00", "ff11",
     &    "ff12", "ff13", "ff14", "ff15", "ff22", "ff23",
     &    "ff24", "ff25", "ff33", "ff34", "ff35", "ff44",
     &    "ff45", "ff55", "ff001", "ff002", "ff003", "ff004",
     &    "ff005", "ff111", "ff112", "ff113", "ff114",
     &    "ff115", "ff122", "ff123", "ff124", "ff125",
     &    "ff133", "ff134", "ff135", "ff144", "ff145",
     &    "ff155", "ff222", "ff223", "ff224", "ff225",
     &    "ff233", "ff234", "ff235", "ff244", "ff245",
     &    "ff255", "ff333", "ff334", "ff335", "ff344",
     &    "ff345", "ff355", "ff444", "ff445", "ff455",
     &    "ff555", "ff0000", "ff0011", "ff0012", "ff0013",
     &    "ff0014", "ff0015", "ff0022", "ff0023", "ff0024",
     &    "ff0025", "ff0033", "ff0034", "ff0035", "ff0044",
     &    "ff0045", "ff0055", "ff1111", "ff1112", "ff1113",
     &    "ff1114", "ff1115", "ff1122", "ff1123", "ff1124",
     &    "ff1125", "ff1133", "ff1134", "ff1135", "ff1144",
     &    "ff1145", "ff1155", "ff1222", "ff1223", "ff1224",
     &    "ff1225", "ff1233", "ff1234", "ff1235", "ff1244",
     &    "ff1245", "ff1255", "ff1333", "ff1334", "ff1335",
     &    "ff1344", "ff1345", "ff1355", "ff1444", "ff1445",
     &    "ff1455", "ff1555", "ff2222", "ff2223", "ff2224",
     &    "ff2225", "ff2233", "ff2234", "ff2235", "ff2244",
     &    "ff2245", "ff2255", "ff2333", "ff2334", "ff2335",
     &    "ff2344", "ff2345", "ff2355", "ff2444", "ff2445",
     &    "ff2455", "ff2555", "ff3333", "ff3334", "ff3335",
     &    "ff3344", "ff3345", "ff3355", "ff3444", "ff3445",
     &    "ff3455", "ff3555", "ff4444", "ff4445", "ff4455",
     &    "ff4555", "ff5555", 
     &    "ff00001", "ff00002", "ff00003", "ff00004", "ff00005", 
     &    "ff00111", "ff00112", "ff00122", "ff00222", "ff00113", 
     &    "ff00123", "ff00223", "ff00133", "ff00233", "ff00333", 
     &    "ff00114", "ff00124", "ff00224", "ff00134", "ff00234", 
     &    "ff00334", "ff00144", "ff00244", "ff00344", "ff00444", 
     &    "ff00115", "ff00125", "ff00225", "ff00135", "ff00235", 
     &    "ff00335", "ff00145", "ff00245", "ff00345", "ff00445", 
     &    "ff00155", "ff00255", "ff00355", "ff00455", "ff00555", 
     &    "ff11111", "ff11112", "ff11122", "ff11222", "ff12222", 
     &    "ff22222", "ff11113", "ff11123", "ff11223", "ff12223", 
     &    "ff22223", "ff11133", "ff11233", "ff12233", "ff22233", 
     &    "ff11333", "ff12333", "ff22333", "ff13333", "ff23333", 
     &    "ff33333", "ff11114", "ff11124", "ff11224", "ff12224", 
     &    "ff22224", "ff11134", "ff11234", "ff12234", "ff22234", 
     &    "ff11334", "ff12334", "ff22334", "ff13334", "ff23334", 
     &    "ff33334", "ff11144", "ff11244", "ff12244", "ff22244", 
     &    "ff11344", "ff12344", "ff22344", "ff13344", "ff23344", 
     &    "ff33344", "ff11444", "ff12444", "ff22444", "ff13444", 
     &    "ff23444", "ff33444", "ff14444", "ff24444", "ff34444", 
     &    "ff44444", "ff11115", "ff11125", "ff11225", "ff12225", 
     &    "ff22225", "ff11135", "ff11235", "ff12235", "ff22235", 
     &    "ff11335", "ff12335", "ff22335", "ff13335", "ff23335", 
     &    "ff33335", "ff11145", "ff11245", "ff12245", "ff22245", 
     &    "ff11345", "ff12345", "ff22345", "ff13345", "ff23345", 
     &    "ff33345", "ff11445", "ff12445", "ff22445", "ff13445", 
     &    "ff23445", "ff33445", "ff14445", "ff24445", "ff34445", 
     &    "ff44445", "ff11155", "ff11255", "ff12255", "ff22255", 
     &    "ff11355", "ff12355", "ff22355", "ff13355", "ff23355", 
     &    "ff33355", "ff11455", "ff12455", "ff22455", "ff13455", 
     &    "ff23455", "ff33455", "ff14455", "ff24455", "ff34455", 
     &    "ff44455", "ff11555", "ff12555", "ff22555", "ff13555", 
     &    "ff23555", "ff33555", "ff14555", "ff24555", "ff34555", 
     &    "ff44555", "ff15555", "ff25555", "ff35555", "ff45555", 
     &    "ff55555" /
	end
