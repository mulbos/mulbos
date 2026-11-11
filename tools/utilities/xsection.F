* xsection.F
* routines to compute the cross-section
* parts of this file are taken from FormCalc-6.0 of Thomas Hahn
* last modified 21 Apr 2010

************************************************************************
	subroutine ProcessIni(fail)
	use global_vars
	use global_xsection
	implicit none
	integer fail
	character*20 pol
	integer i, c, bits, df

***[todo: cleaned up
#define SCALAR 0
#define FERMION 1
#define PHOTON 2
#define GLUON PHOTON
#define VECTOR 3
	integer type_gen(4)
	data type_gen /0, 1, 2, 3/
***]todo: cleaned up

	integer type(proc_nlegs)
	
	if( proc_nlegs .eq. 4 ) then
	  pol = 'UUUU'
	else if ( proc_nlegs .eq. 5 ) then
	  pol = 'UUUUU'
	else if ( proc_nlegs .eq. 6 ) then
	  pol = 'UUUUUU'
	else if ( proc_nlegs .eq. 7 ) then
	  pol = 'UUUUUUU'
	else if ( proc_nlegs .eq. 8 ) then
	  pol = 'UUUUUUUU'
	else if ( proc_nlegs .eq. 9 ) then
	  pol = 'UUUUUUUUU'
	else if ( proc_nlegs .eq. 10 ) then
	  pol = 'UUUUUUUUUU'
	else
	  print*,"invalid value of proc_nlegs:",proc_nlegs
	  stop
	endif

	df = 2**proc_nlegs
!	print*,"kinematic factor 1, df =",df
	do i = 1, proc_nlegs
	 type(i) = type_gen(2) 
        enddo

	do i = 1, proc_nlegs
	  if( type(i) .eq. SCALAR ) then
	    c = ichar('S')
	    bits = 1
	  else
	    c = ichar(pol(i:i))
	    if( c .eq. ichar('+') ) then
	      bits = 4
	    else if( c .eq. ichar('-') ) then
	      bits = 1
	    else
	      c = iand(c, 223)
	      if( c .eq. ichar('U') ) then
	        bits = 7
	      else if( c .eq. ichar('T') ) then
	        bits = 5
	      else if( c .eq. ichar('L') ) then
	        bits = 2
	      else
	        bits = 0
	      endif

	      if( type(i) .ne. VECTOR ) then
	        bits = iand(bits, -3)
	        c = ior(c, 1)
	      endif

c	      if( bits .eq. 0 ) then
c	        Error("Invalid polarization for leg "/ /Digit(i))
c	        return
c	      endif

	      if( i .le. 2 ) df = df*(ibits(bits, 0, 1) +
     &          ibits(bits, 1, 1) + ibits(bits, 2, 1))
	    endif
	  endif

	  pol(i:i) = char(c)
	  helicities = helicities*8 + bits

      if (lNLO_FIN_IND) then
          if (i .ne. 3) helicities0 = helicities0*8 + bits
          if (i .ne. 4) helicities1 = helicities1*8 + bits
      else
          if (i .lt. proc_nlegs) helicities0 = helicities0*8 + bits 	  
      endif
	enddo

C particles
c        Nlepton = 4
c        Nlep_charged = 3

	avgfac = 2**proc_nlegs*fac_color*fac_identical/df
!	print*,"kinematic factor 2, df =",df

      if (lVIRT.or.lGLUGLU) then
	call loopints_init
      endif

	call ModelDefaults

! on-shell mapping method: osmap_case
! possible values for NLO QCD: 1 (default, Denner ZZ 2021), 2 (Denner RACOONWW 2000)
! possible values for NLO EW: 1 (default, Denner ZZ 2021)

      call InitPdfs
      call KinDistSetup

	fail = 0
	end
************************************************************************
** IntegratedCS computes the integrated cross-section.

	subroutine IntegratedCS
	use global_vars
	use global_xsection
	use global_virt
	implicit none

#include "model.h"
	external DifferentialCS

	double precision result(NCOMPO), error(NCOMPO)
	integer fail, ndim

        afcutcount = 0
        badcount = 0
        nancount = 0
        gramcount4 = 0 ! number of bad 4-point gram determinants

      if (lVIRT.or.lGLUGLU) then
       call set_virt_params()
      endif

	call ModelConstIni(fail)
	if( fail .ne. 0 ) goto 999

* Ninh: flags (scale, parameters): old stuff, 
* flags1 (on scale, parameters).
* Warning: these flags are OK for NCOMPO=1 or 2.
	flags = 1 + 2*(NCOMPO-1)
	flags1 = flags

***[number of integration variables
        ndim = (proc_nlegs-2)*3-3 ! LO, VIRT, FIN case
	if (lNLO_COLL_RAD.or.lNLO_COLL_IND) then
         ndim = ndim+1 ! coll case
        endif
***]number of integration variables
	
c	call SetEnergy(fail, ecm)
c	if( fail .ne. 0 ) goto 999

	call ModelDigest
#ifdef YES_Coupl
* first call to setup the paramaters
	print*,'Use MADGRAPH'
#ifdef YES_AUTODIPOLE
	print*,'Use AUTODIPOLE'
#endif
#ifdef YES_MadDipole
	print*,'Use MadDipole'
#endif
        call ModelCoupl(fail)
c        call setpara('param_card.dat',.true.)
c        call printout()
#endif

c#define SHOW print 100,
c100	format("|# ", A, "=", F10.4, SP, F10.4, " I")

      if (lVIRT) then
       if (lNLO_QCD) then
        call CalcRenConst_QCD
       else if (lNLO_EW) then
        call CalcRenConst_EW
       endif
      endif
      
      if (lVIRT.or.lGLUGLU) then
       call clearcache
      endif

! phase space integration
      call BASES2(ndim, result, error)

***[print out
c102	format("|+   ", NCOMPO G24.15)
c	print 102, result
c	print 102, error

      if (lVIRT.or.lGLUGLU) then
	call loopints_info
      endif

	call flush(6)
***]print out

999	continue
	end
***********************************************************************
	subroutine GetMomenta(pp, nn)
	use global_vars
	use global_mom_xpdf
	implicit none
	integer i, j, nn
	double precision pp(0:3,nn), temp(0:3)

        do i=1, nn

         do j = 0,3
          pp(j,i) = pext(i,j+1)
         enddo

***[fix momenta
#ifdef CHECK_POINT_MIRROR
* space-time reflection:
	pp(0,i) = pp(0,i)
	pp(1,i) = -pp(1,i)
	pp(2,i) = -pp(2,i)
	pp(3,i) = -pp(3,i)
#endif
***]fix momenta
	enddo

***[fix momenta
#ifdef CHECK_POINT_CP
* CP symmetry:
* for U_bar(p1) + U(p2) -> W+(p3) + W-(p4) + Z(p5)
* CP conjugated process: U(-p1) + U_bar(-p2) -> W-(-p3) + W+(-p4) + Z(-p5).
* In the CMSqq we have p1 = -p2
	temp(0) = pp(0,3)
	pp(0,3) = pp(0,4)
	pp(0,4) = temp(0)
	pp(0,5) = pp(0,5)
	do i = 1, 3
	temp(i) = pp(i,3)
	pp(i,3) = -pp(i,4)
	pp(i,4) = -temp(i)
	pp(i,5) = -pp(i,5)
	enddo
#endif
***]fix momenta

	end
***********************************************************************
      double precision function dot(p1,p2) 
      implicit none
      double precision p1(0:3),p2(0:3)
      dot=p1(0)*p2(0)-p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)
      end
***********************************************************************
      double precision function dot_p2(p) 
! calculate p^2:
      implicit none
      double precision p(0:3)
      dot_p2=sqrt(p(1)**2+p(2)**2+p(3)**2)
      dot_p2=(p(0)-dot_p2)*(p(0)+dot_p2)
      end
***********************************************************************
      double precision function inv_mass2(p1,p2) 
! calculate squared invariant mass
      implicit none
      double precision p1(0:3),p2(0:3),p(0:3)
      p(:) = p1(:) + p2(:)
      inv_mass2=sqrt(p(1)**2+p(2)**2+p(3)**2)
      inv_mass2=(p(0)-inv_mass2)*(p(0)+inv_mass2)
      end
***********************************************************************
      double precision function inv_mass(p1,p2) 
! calculate invariant mass
      implicit none
      double precision p1(0:3),p2(0:3),p(0:3)
      p(:) = p1(:) + p2(:)
      inv_mass=sqrt(p(1)**2+p(2)**2+p(3)**2)
      inv_mass=(p(0)-inv_mass)*(p(0)+inv_mass)
      inv_mass=sqrt(inv_mass)
      end
***********************************************************************
      double precision function inv_mass_p(p) 
! calculate sqrt(p^2):
      implicit none
      double precision p(0:3)
      inv_mass_p=sqrt(p(1)**2+p(2)**2+p(3)**2)
      inv_mass_p=(p(0)-inv_mass_p)*(p(0)+inv_mass_p)
      inv_mass_p=sqrt(inv_mass_p)
      end
***********************************************************************
      double precision function inv_mass_p2(p) 
! calculate p^2:
      implicit none
      double precision p(0:3)
      inv_mass_p2=sqrt(p(1)**2+p(2)**2+p(3)**2)
      inv_mass_p2=(p(0)-inv_mass_p2)*(p(0)+inv_mass_p2)
      end
***********************************************************************
      double precision function dot_space(p1,p2) 
      implicit none
      double precision p1(1:3),p2(1:3)
      dot_space=p1(1)*p2(1)+p1(2)*p2(2)+p1(3)*p2(3)
      end
***********************************************************************
	subroutine GetMomenta_BORN(qext)
	use global_vars
	use global_mom_xpdf
	implicit none
	integer i, j, i1
C Note: only qext(proc_nlegs - 1,4) is needed, but we keep this 
c for historical reason (same as the old code in ppWZ_vbf).
	double precision qext(proc_nlegs,4)

***[Exact momenta in the partonic CMS
        DO j = 1, 4
        do i = 1, 2
        qext(i,j) = pext(i,j)
        enddo

      if (lNLO_FIN_IND) then
        do i = 3, proc_nlegs-1
        qext(i,j) = pext(i+1,j)
        enddo
        qext(proc_nlegs,j) = pext(3,j)
      else
        do i =3, proc_nlegs
        qext(i,j) = pext(i,j)
        enddo
      endif
        ENDDO
***]Exact momenta in the partonic CMS

	end
***********************************************************************
	subroutine GetMomentaFC_OS_BORN_Tilde_obsolete(qext)
	use global_vars
	use global_mom_xpdf
	implicit none
	integer i, j
C Note: only qext(proc_nlegs - 1,4) is needed, but we keep this 
c for historical reason (same as the old code in ppWZ_vbf).
	double precision qext(proc_nlegs,4), kext ! important, array has to exactly match with other calls.


!!![matching util.h
#define SPEC_K 2
#define SPEC_E 3
#define SPEC_EX 10
#define SPEC_EY 11
#define SPEC_EZ 12
! do not remove vec, kcomp (for some reason, the code does not work without these declarations)
! the parameter nlegs_max = 20 has to match VecSet()
	double complex vec(2,2,8,0:20) ! FC6
	common /vectors/ vec
	double complex kcomp(32,0:20)
	equivalence (vec, kcomp)
	double precision momspec(16,0:20)
	common /momenta/ momspec
!!!]matching util.h

***[on-shell momenta
C These are momenta that enter the tilde (or reduced) amplitudes 
C needed for the subtraction terms.
        do i=1, 2
        kext=momspec(SPEC_K,i)
        qext(i,1)=momspec(SPEC_E,i)
        qext(i,2)=momspec(SPEC_EX,i)*kext
        qext(i,3)=momspec(SPEC_EY,i)*kext
        qext(i,4)=momspec(SPEC_EZ,i)*kext
	enddo

        do i=3, proc_nlegs-1
        j = i + 2 ! this shift is due to 2 gauge bosons
C ordering: q1 q2 > V3 V4 > L5 L6 L7 L8
        kext=momspec(SPEC_K,j)
        qext(i,1)=momspec(SPEC_E,j)
        qext(i,2)=momspec(SPEC_EX,j)*kext
        qext(i,3)=momspec(SPEC_EY,j)*kext
        qext(i,4)=momspec(SPEC_EZ,j)*kext
	enddo

C This setting to zero is important for cut and histogram routines
        do j = 1, 4
        qext(proc_nlegs,j) = 0
        enddo
***]on-shell momenta
! Ninh checking
c        do i = 1,proc_nlegs
c         print*,"tilde p =",i,qext(i,:)
c        enddo
c        stop

	end
***********************************************************************
	subroutine GetMomenta_NLO_FIN(qext)
	use global_vars
	use global_mom_xpdf
	implicit none
	integer i, j, i1
	double precision qext(proc_nlegs,4)

!***[Exact momenta in the partonic CMS
        DO j = 1, 4
        do i = 1, proc_nlegs
        qext(i,j) = pext(i,j)
        enddo
        ENDDO
!***]Exact momenta in the partonic CMS
	end
***********************************************************************
	subroutine GetMomenta_NLO_FIN2(qext)
	use global_vars
	use global_mom_xpdf	
	implicit none
	integer i, j, i1
	double precision qext(proc_nlegs,4)

!***[Exact momenta in the partonic CMS
        DO j = 1, 4
        do i = 1, 2
        qext(i,j) = pext(i,j)
        enddo

        do i=3, proc_nlegs - 1 ! 4 leptons
        qext(i,j) = pext(i+1,j)
        enddo

        qext(proc_nlegs,j) = pext(3,j) ! final-state quark
        ENDDO
!***]Exact momenta in the partonic CMS
	end
***********************************************************************
	subroutine GetMass(e, p, mass)
	implicit none
	double precision e, p, mass

	mass = (e - p)*(e + p)

	if (mass .ge. 0) then
	 mass = sqrt(mass)
	elseif (abs(mass) .lt. 1D-8*p**2) then
	 mass = 0
	else
c The mass cannot be too negative. 
c Stop here for safety reason. 
         print*,'Warning in GetMass(): negative mass set to 0'
	 print*,'squared mass =', mass

	 stop
	 mass = 0
	endif

	end
***********************************************************************
	subroutine PrintPSP()
	use global_vars
	use global_mom_xpdf	
	implicit none
	integer i, j
        double precision kext

        print*,'from PrintPSP: true momenta (pext):'
        print*,'number of legs =',proc_nlegs
        do i=1,proc_nlegs
        print*,'[E, px, py, pz] of particle',i
        print*,pext(i,1),pext(i,2),pext(i,3),pext(i,4)
        print*,'mass2=',pext(i,1)**2-pext(i,2)**2-pext(i,3)**2-pext(i,4)**2
	enddo

        do i=1,4
	kext = 0
	do j=1,2
	kext = kext + pext(j,i)
	enddo
	do j=3, proc_nlegs
	kext = kext - pext(j,i)
	enddo
	print*,'Check zero-sum EP: ',i,kext
	enddo
	
	print*,"x1_pdf =",x1_pdf
	print*,"x2_pdf =",x2_pdf
	print*,"y_coll =",y_coll
	end
***********************************************************************
	subroutine PrintPSP_OS()
	use global_vars
	use global_mom_xpdf	
	implicit none
	integer i, j
        double precision kext

        print*,'PrintPSP_OS: on-shell projected momenta (oext):'
        print*,'number of legs =',proc_nlegs_os
        do i=1,proc_nlegs_os
        print*,'[E, px, py, pz] of particle',i
        print*,oext(i,1),oext(i,2),oext(i,3),oext(i,4)
        print*,'mass2=',oext(i,1)**2-oext(i,2)**2-oext(i,3)**2-oext(i,4)**2
	enddo

        print*,'check for lepton production:'
        do i=1,4
	kext = 0
	do j=1,2
	kext = kext + oext(j,i)
	enddo
	do j=2 + 3, proc_nlegs_os
	kext = kext - oext(j,i)
	enddo
	print*,'Check zero-sum EP: ',i,kext
	enddo

	end
***********************************************************************
        subroutine PrintPSP_Decay(pp, qq)
        implicit none
        integer i, j

   

        double precision pp(4), qq(4,2), kext, dot

        print*,'[E, px, py, pz] of parent particle'
        print*,pp(1),pp(2),pp(3),pp(4)
        print*,"mass2 =",dot(pp(1),pp(1))
        do i=1,2
        print*,'[E, px, py, pz] of child particle',i
        print*,qq(1,i),qq(2,i),qq(3,i),qq(4,i)
        print*,"mass2 =",dot(qq(1,i),qq(1,i))
        enddo

        print*,'Check decay zero-sum EP for V > l1 l2'
        do i=1,4
        kext = 0
        kext = kext + pp(i)
        do j=1,2
        kext = kext - qq(i,j)
        enddo
        print*,'Check zero-sum EP: ',i,kext
        enddo

        end
********************************************************************
* Use BASES for phase space integration

	subroutine BASES2(ndim, result, error)
	use global_vars
	use global_xsection
	use global_dists
	implicit none

	!integer LEN_TRIM
	!external LEN_TRIM

	integer ndim
	double precision result(NCOMPO), error(NCOMPO)
        double precision FUNC
        EXTERNAL FUNC
	
	double precision prob(NCOMPO)
c        double precision HI,xdis
	
c        character*200 argi1,argi2,argix
        character*512 filen
        integer fail
        character*10 version	
	integer i
	
      Integer MXDIM,NDIMEN,NWILD,NCALL
      PARAMETER (MXDIM = 50)
      Integer IG(MXDIM),ITMX1,ITMX2
      
      double precision XL(MXDIM),XU(MXDIM)
      double precision ACC1,ACC2,ESTIM,SIGMA,CTIME
      Integer IT1, IT2
      Integer ISEED_FLAG, ISEED_IN
      COMMON /BPARM0/ ISEED_FLAG, ISEED_IN 
      COMMON /BPARM1/ XL,XU,NDIMEN,NWILD,IG,NCALL
      COMMON /BPARM2/ ACC1,ACC2,ITMX1,ITMX2
      
      integer lu, N_tot

C pass random-number seed to BASES
      ISEED_FLAG = 1 ! choose 1 to change the default Seed value, else otherwise
      ISEED_IN = seed ! read in from run argument obtained in main.F
      
*===> Initialization of BASES by calling BSINIT
           CALL BSINIT

           NDIMEN=ndim
           NWILD = NWILDin
           NCALL = 2**(N_stat)
           ITMX1 = ITMX1in
           ACC1  = ACC1in
           ITMX2 = ITMX2in
           ACC2  = ACC2in

           do i = 1, NDIMEN
             XL(i)= 0D0
             XU(i)= 1D0
           enddo

        if (lHist_Fill) then
***{Initialization of Histograms
      if (lWZ) then
       call hist_init_WZ()
      else if (lWW) then
       call hist_init_WW()
      else if (lZZ) then
       call hist_init_ZZ()
      else if (lWWW) then
       call hist_init_WWW()
      endif
***}Initialization of Histograms
        endif

        do icomp=NCOMPOmin, NCOMPO 
        call BASES(FUNC, ESTIM, SIGMA, CTIME, IT1, IT2 )
        result(icomp)=ESTIM
	error(icomp)=SIGMA
c	prob(icomp)=dabs(SIGMA/ESTIM)
	
	print*
	print*,'-------------------------------------------------------'
	write(*,14)'Energ, Xsect[pb], error : ',ecm,result(icomp),
     >	 error(icomp)
14      format (A, E9.3, '  ', E14.8, '  ', E14.8)     
        print*,"Cross_section_unit: pb"
        print*,"Cross_section_value:", result(icomp)
        print*,"Cross_section_error:", error(icomp)
	print*,'-------------------------------------------------------'
 
***{for histograms
        if (dist_case .eq. 1 .or. dist_case .eq. 3) then
        lu = 6 ! choose a logical unit value
        CALL BHPLOT( lu ) !print histograms in 0000001 output file
        endif

        if (dist_case .le. 2) then      
*===============================================*
* save the content of Histograms into data files
*===============================================*

      do i = 1, NHISMAX
cc      filen='data_dist/'//filem(i)
      filen=dir_save(1:LEN_TRIM(dir_save))//filem(i)
      filem(i)=filen(1:LEN_TRIM(filen))
      enddo

* write histograms to files
      call bshbok(filem(1))
        endif
***}for histograms

        enddo

      N_tot = NCALL*(ITMX1 + ITMX2)
      print*,'Number of before-cut points:',N_tot
      print*,'Number of after-cut points:',afcutcount,DBLE(afcutcount)/N_tot
      print*,'Number of NaN points:',nancount,DBLE(nancount)/N_tot
      print*,'Number of bad points:',badcount,DBLE(badcount)/N_tot
      print*,'Number of bad gram det. points (boxes):',gramcount4,DBLE(gramcount4)/N_tot

      if (lLHAPDF) then
        call GetLHAPDFVersion(version)
        print*,'LHAPDFVersion:',version
      endif

        return
	end
************************************************************************
	double precision function FUNC(XX)
	use global_vars
	use global_xsection
	implicit none



      Integer MXDIM,NDIMEN,NWILD,NCALL
      PARAMETER (MXDIM = 50)
      Integer IG(MXDIM)
      double precision XL(MXDIM),XU(MXDIM)

      COMMON /BPARM1/ XL,XU,NDIMEN,NWILD,IG,NCALL	

	integer i
	double precision x(NDIMEN), result(NCOMPO)
        double precision XX(MXDIM)

!	print*,"hello from FUNC() in xsection.F"      

        do i=1,NDIMEN
	 x(i)=XX(i)
	enddo
        Call DifferentialCS(NDIMEN, x, NCOMPO, result)
        FUNC=result(icomp)
	
	return
	end
***********************************************************************
      double precision function LambdaX(x,y,z)
      implicit none
      double precision x,y,z
      LambdaX=x**2+y**2+z**2-2*x*y-2*x*z-2*z*y
      end
********************* BOOST **************************
* Warning: rp0, rpz are changed.
        subroutine boost0(rgam, rbetgam, rp0, rpx, rpy, rpz)    
        
        implicit none
        double precision rpx,rpy,rpz,rp0,temp
        double precision rgam, rbetgam
        temp = rp0
        rp0 = rgam*rp0 - rbetgam*rpz
        rpz = -rbetgam*temp + rgam*rpz
        end
***********************************************************************
      subroutine GetLeptonMomProjOS(pV,pi1,pi2,po1,po2)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C input: pV(0:3) is the momentum of the gauge boson, 
c pi1 and pi2 are 2 reference momenta of OS decay leptons; 
c the asumption po1(mu) = c*pi1(mu) is used.

C output: po1 and po2 satisfying pV = po1 + po2, po1^2 = po2^2 = 0.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      double precision pV(0:3), pi1(0:3), pi2(0:3)
      double precision po1(0:3), po2(0:3), c, dot
      integer i

      c = dot(pV(0),pV(0))/(2*dot(pi1(0),pV(0)))
      do i = 0, 3
       po1(i) = c*pi1(i)
       po2(i) = pV(i) - po1(i)
      enddo

c      print*,"mass2 reference lepton 1 =",dot(pi1(0),pi1(0))/pi1(0)**2
c      print*,"mass2 reference lepton 2 =",dot(pi2(0),pi2(0))/pi2(0)**2
c      print*,"mass2 gauge boson =",dot(pV(0),pV(0))
c      print*,"mass2 lepton 1 =",dot(po1(0),po1(0))/po1(0)**2
c      print*,"mass2 lepton 2 =",dot(po2(0),po2(0))/po2(0)**2
      END
***********************************************************************
***[regular subtraction function for (n+1) term
      subroutine DoTildeAmpHis_OffS(tp,iflag,Gij,Jacb)
C tp(nleg0,4) are the exact qq > 4l momenta.
C with on-shell projection after tilde (reduced) kinematics
C The input momenta tp(nleg0,4) are defined in the same reference frame as 
C the momenta entering the (N+1) amplitudes. 
C Output: return new Gij with the amp multiplied, or set Gij = 0
      use global_vars
      use global_mom_xpdf 
      implicit none
#include "model.h"
        integer iflag
	integer i, j, i1, j1, j2, nlegs, cfail, nleg0
        parameter (nlegs = 8) ! number of legs for child q1 q2 > V3 V4 > l5 l6 l7 l8
        parameter (nleg0 = 6) ! number of legs for child q1 q2 > l3 l4 l5 l6

	double precision Gij, temp, temp1, temp2	
	double precision result(2)
        double precision tp(nleg0,4)
        double precision Vext(4,nleg0), qcm(4,nleg0)
        double precision tp1(4,nlegs) ! in the VV cms frame
        double precision tp1_pp(4,nlegs) ! in the Lab frame
        double precision qX(4), qX_inv(4)
        double precision MVint2(2) ! the off-shell invariant mass of the gauge bosons
!	double precision p, ex,ey,ez, psq
	double precision Jacb, MVV

        double precision rext(4,6), rm2(2)

	double precision qext(proc_nlegs,4) ! important, array has to exactly match with FillHIS routines

       do j1 = 1,nleg0 
       do j2 = 1,4
        Vext(j2,j1) = tp(j1,j2)
       enddo
       enddo
       
***[boost from partonic CMS to the X (= V1 + V2) rest frame
      do i = 1, 4
       qX(i) = tp(1,i) + tp(2,i) ! sum over initial particles
      enddo
      
!!![check input momenta
c      do i = 1, nleg0
c       print*,"tilde tp =",i,tp(i,:)
c      enddo
      ! VV cms frame ? NO
c      print*,"tilde tp12 =",tp(1,:)+tp(2,:)
!!!]check input momenta

! This cut has effect on final-initial terms
      MVV = sqrt(qX(1)**2 - qX(2)**2 - qX(3)**2 - qX(4)**2)
        if (MVV .le. (mass_VV(1) + mass_VV(2)) ) then 
         Gij = 0D0
         return
        endif

      qX_inv(1) = qX(1)
      do i = 2, 4
      qX_inv(i) = -qX(i)
      enddo
      do i = 1, nleg0
      call BOOSTN(Vext(1,i),qX_inv(1),qcm(1,i))
      enddo
***]boost from partonic CMS to the X (= V1 + V2) rest frame

!      print*,"from DoTildeAmpHis_OffS(): end boost 1"

***[momenta for cut and hist in the partonic CMS in special format
       do j1 = 1,nleg0 
       do j2 = 1,4
!!!        qext(j1,j2) = qcm(j2,j1) ! in the WZ rest frame
! Better use the original tilde mom, because we have to boost 
! back to the Lab frame using the (N+1) boost parameters.
! If we take qext() to be in the WZ rest frame, we have to boost 
! back in 2 steps, which is not done in DoCut_parCM(). 
        qext(j1,j2) = tp(j1,j2) ! the original tilde mom in the partonic CMS
       enddo
       enddo
C This setting to zero is important for cut and histogram routines
        do j = 1, 4
        qext(proc_nlegs,j) = 0
        enddo
        ! store qext(proc_nlegs,4) > pext_tilde(nlegs_max,4)
        pext_tilde = 0D0
        do i = 1,proc_nlegs
         pext_tilde(i,:)=qext(i,:)
        enddo
***]momenta for cut and hist in the partonic CMS in special format

***[tp1() is exact momenta in the VV CMS
! for child q1 q2 > V3 V4 > l5 l6 l7 l8
       rm2(1) = mass_VV(1)**2
       rm2(2) = mass_VV(2)**2
! initial particles
       do j1 = 1,2 
       do j2 = 1,4
        tp1(j2,j1) = qcm(j2,j1)
       enddo
       enddo
! final leptons
       do j1 = 3,nleg0
       j2 = j1 + 2
       do i = 1,4
        tp1(i,j2) = qcm(i,j1)
       enddo
       enddo
! intermediate gauge bosons
       do j2 = 1,4
        tp1(j2,3) = qcm(j2,3)+qcm(j2,4)
        tp1(j2,4) = qcm(j2,5)+qcm(j2,6)
       enddo
***]tp1() is exact momenta in the VV CMS

! This is Lorentz invariant   
       temp1 = sqrt(tp1(2,3)**2 + tp1(3,3)**2 + tp1(4,3)**2)
       MVint2(1) = (tp1(1,3)-temp1)*(tp1(1,3)+temp1)
       temp2 = sqrt(tp1(2,4)**2 + tp1(3,4)**2 + tp1(4,4)**2)
       MVint2(2) = (tp1(1,4)-temp2)*(tp1(1,4)+temp2)

***[checking momenta
#ifdef CHECK_MOMENTA
       if (MVint2(1).lt.0D0) then
c	call PrintPSP()
c        call PrintPSP_OS()
        print*,"exact qq > VV > 4l, tp1:"
        do i1 = 1, nlegs
        print*,i1,tp1(1,i1),tp1(2,i1),tp1(3,i1),tp1(4,i1)
        enddo
        print*,"MVint2(1) =",MVint2(1)
        stop
       endif
#endif
***]checking momenta

!      print*,"from DoTildeAmpHis_OffS(): begin OSProject"

***[on-shell mapping
      call OSProject(osmap_case,tp1(1,5),rm2(1),rext(1,1))
       do j1 = 3,nlegs 
       i1 = j1 - 2
       do j2 = 1,4
        tp1(j2,j1) = rext(j2,i1) ! tp1() is now OS momenta
       enddo
       enddo
***]on-shell mapping

!      print*,"from DoTildeAmpHis_OffS(): end OSProject"

      ! store tp1(4,nlegs) > oext_tilde_vv(nlegs_os_max,4)
      oext_tilde_vv = 0D0
      do i = 1, nlegs
       oext_tilde_vv(i,:) = tp1(:,i)
      enddo

!      print*,"from DoTildeAmpHis_OffS(): begin boost 2"

!!![get os momenta in the Lab frame
! for polarizations in the Lab frame:
! boost tp1(4,nlegs) back in 2 steps: 
! (1): boost back to the partonic CMS using qX()
! (2): boost back to the Lab frame using (rgam_boost, rbetgam_boost)
      do i = 1, nlegs
! boost back from VV CMS to parent partonic CMS
       call BOOSTN(tp1(1,i),qX(1),tp1_pp(1,i))
      enddo
      ! store tp1_pp(4,nlegs) > oext_tilde(nlegs_os_max,4)
      oext_tilde=0D0
      do i = 1, nlegs
       oext_tilde(i,:) = tp1_pp(:,i)
      enddo
      
      do i = 1, nlegs
! boost from partonic CMS to the LAB frame
      call boost0(rgam_boost, -rbetgam_boost,
     &  tp1_pp(1,i), tp1_pp(2,i), tp1_pp(3,i), tp1_pp(4,i))
      enddo
      ! store tp1_pp(4,nlegs) > oext_tilde_pp(nlegs_os_max,4)
      oext_tilde_pp=0D0
      do i = 1, nlegs
       oext_tilde_pp(i,:) = tp1_pp(:,i)
      enddo
!!!]get os momenta in the Lab frame

!      print*,"from DoTildeAmpHis_OffS(): end boost 2"

        call Do_VecSet_tilde()

***[impose cuts
C the momenta here has to be the tilde kinematics
        cfail = 1 ! initialization
c        print*,"qext =",qext ! Ninh check
	call DoCut_parCM(cfail, qext(1,1))
	if ( cfail .eq. 0 ) then
         Gij = 0D0
         return
        endif
***]impose cuts

!      print*,"from DoTildeAmpHis_OffS(): begin amp2"

	call P0SquaredME(result, MVint2(1), 1)
        Gij =  result(1)*Gij
       IF (lscaleFR_dyn) then ! dynamical scale
        call scaleFR(.true.) ! choose true for subtraction terms
        if (lNLO_QCD) then
         call GetAlfas_dyn()
        endif
        fac_pdf = 1d0
        call MultiplyPdfs(fac_pdf)
       ELSE ! fixed scale
        if (.not.lpdf_fixed) then
         fac_pdf = 1d0
         call MultiplyPdfs(fac_pdf)
         lpdf_fixed = .true.
        endif
       endif

      if (lNLO_QCD) then
       Gij = Alfas*Gij*fac_pdf
      else if (lNLO_EW) then
       Gij = Alfa*Gij*fac_pdf
      endif

!      print*,"from DoTildeAmpHis_OffS(): end amp2"

        if (lHist_Fill) then
C the momenta here has to be the tilde kinematics
        call FillHIS(Gij, Jacb, qext)
        endif

!      print*,"from DoTildeAmpHis_OffS(): end"

       END
***]regular subtraction function for (n+1) term
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!![for photon-photon contribution of WW case
      subroutine DoTildeAmpHis_OffS_aa(tp,iflag,Gij,Jacb,fac1,fac2,tpp)
C tp(nleg0,4) are the exact qq > 4l momenta.
C with on-shell projection after tilde (reduced) kinematics
C The input momenta tp(nleg0,4) are defined in the same reference frame as 
C the momenta entering the (N+1) amplitudes. 
C Output: return new Gij with the amp multiplied, or set Gij = 0
      use global_vars
      use global_mom_xpdf 
      implicit none
#include "model.h"
        integer iflag
	integer i, j, i1, j1, j2, nlegs, cfail, nleg0
        parameter (nlegs = 8) ! number of legs for child q1 q2 > V3 V4 > l5 l6 l7 l8
        parameter (nleg0 = 6) ! number of legs for child q1 q2 > l3 l4 l5 l6

	double precision Gij, temp, temp1, temp2	
	double precision result(2)
        double precision tp(nleg0,4)
        double precision Vext(4,nleg0), qcm(4,nleg0)
        double precision tp1(4,nlegs) ! in the VV cms frame
        double precision tp1_pp(4,nlegs) ! in the Lab frame
        double precision qX(4), qX_inv(4)
        double precision MVint2(2) ! the off-shell invariant mass of the gauge bosons
!	double precision p, ex,ey,ez, psq
	double precision Jacb, MVV

        double precision rext(4,6), rm2(2)

	double precision qext(proc_nlegs,4) ! important, array has to exactly match with FillHIS routines
	double complex epsp(-1:1,4), epspc(-1:1,4) ! polarization vectors
	double complex ctemp, amp(2)
	integer flags ! for the reduced aa amplitude
	double precision Hij, fac1, fac2
	double precision tpp(4) ! exact momenta in the partonic CMS
	double precision tpp_vv(4) ! exact momenta in the VV CMS
	double precision tpp_pp(4) ! exact momenta in the Lab frame
	double precision tqq(4)

       do j1 = 1,nleg0 
       do j2 = 1,4
        Vext(j2,j1) = tp(j1,j2)
       enddo
       enddo
       
***[boost from partonic CMS to the X (= V1 + V2) rest frame
      do i = 1, 4
       qX(i) = tp(1,i) + tp(2,i) ! sum over initial particles
      enddo
      
!!![check input momenta
c      do i = 1, nleg0
c       print*,"tilde tp =",i,tp(i,:)
c      enddo
      ! VV cms frame ? NO
c      print*,"tilde tp12 =",tp(1,:)+tp(2,:)
!!!]check input momenta

! This cut has effect on final-initial terms
      MVV = sqrt(qX(1)**2 - qX(2)**2 - qX(3)**2 - qX(4)**2)
        if (MVV .le. (mass_VV(1) + mass_VV(2)) ) then 
         Gij = 0D0
         return
        endif

      qX_inv(1) = qX(1)
      do i = 2, 4
      qX_inv(i) = -qX(i)
      enddo
      do i = 1, nleg0
      call BOOSTN(Vext(1,i),qX_inv(1),qcm(1,i))
      enddo
      ! boost tpp():
      call BOOSTN(tpp(1),qX_inv(1),tpp_vv(1))
***]boost from partonic CMS to the X (= V1 + V2) rest frame

***[momenta for cut and hist in the partonic CMS in special format
       do j1 = 1,nleg0 
       do j2 = 1,4
!!!        qext(j1,j2) = qcm(j2,j1) ! in the WZ rest frame
! Better use the original tilde mom, because we have to boost 
! back to the Lab frame using the (N+1) boost parameters.
! If we take qext() to be in the WZ rest frame, we have to boost 
! back in 2 steps, which is not done in DoCut_parCM(). 
        qext(j1,j2) = tp(j1,j2) ! the original tilde mom in the partonic CMS
       enddo
       enddo
C This setting to zero is important for cut and histogram routines
        do j = 1, 4
        qext(proc_nlegs,j) = 0
        enddo
        ! store qext(proc_nlegs,4) > pext_tilde(nlegs_max,4)
        pext_tilde = 0D0
        do i = 1,proc_nlegs
         pext_tilde(i,:)=qext(i,:)
        enddo
***]momenta for cut and hist in the partonic CMS in special format

***[tp1() is exact momenta in the VV CMS
! for child q1 q2 > V3 V4 > l5 l6 l7 l8
       rm2(1) = mass_VV(1)**2
       rm2(2) = mass_VV(2)**2
! initial particles
       do j1 = 1,2 
       do j2 = 1,4
        tp1(j2,j1) = qcm(j2,j1)
       enddo
       enddo
! final leptons
       do j1 = 3,nleg0
       j2 = j1 + 2
       do i = 1,4
        tp1(i,j2) = qcm(i,j1)
       enddo
       enddo
! intermediate gauge bosons
       do j2 = 1,4
        tp1(j2,3) = qcm(j2,3)+qcm(j2,4)
        tp1(j2,4) = qcm(j2,5)+qcm(j2,6)
       enddo
***]tp1() is exact momenta in the VV CMS

! This is Lorentz invariant   
       temp1 = sqrt(tp1(2,3)**2 + tp1(3,3)**2 + tp1(4,3)**2)
       MVint2(1) = (tp1(1,3)-temp1)*(tp1(1,3)+temp1)
       temp2 = sqrt(tp1(2,4)**2 + tp1(3,4)**2 + tp1(4,4)**2)
       MVint2(2) = (tp1(1,4)-temp2)*(tp1(1,4)+temp2)

***[checking momenta
#ifdef CHECK_MOMENTA
       if (MVint2(1).lt.0D0) then
c	call PrintPSP()
c        call PrintPSP_OS()
        print*,"exact qq > VV > 4l, tp1:"
        do i1 = 1, nlegs
        print*,i1,tp1(1,i1),tp1(2,i1),tp1(3,i1),tp1(4,i1)
        enddo
        print*,"MVint2(1) =",MVint2(1)
        stop
       endif
#endif
***]checking momenta

***[on-shell mapping
      call OSProject(osmap_case,tp1(1,5),rm2(1),rext(1,1))
       do j1 = 3,nlegs 
       i1 = j1 - 2
       do j2 = 1,4
        tp1(j2,j1) = rext(j2,i1) ! tp1() is now OS momenta
       enddo
       enddo
***]on-shell mapping

      ! store tp1(4,nlegs) > oext_tilde_vv(nlegs_os_max,4)
      oext_tilde_vv = 0D0
      do i = 1, nlegs
       oext_tilde_vv(i,:) = tp1(:,i)
      enddo

!!![get os momenta in the Lab frame
! for polarizations in the Lab frame:
! boost tp1(4,nlegs) back in 2 steps: 
! (1): boost back to the partonic CMS using qX()
! (2): boost back to the Lab frame using (rgam_boost, rbetgam_boost)
      do i = 1, nlegs
! boost back from VV CMS to parent partonic CMS
       call BOOSTN(tp1(1,i),qX(1),tp1_pp(1,i))
      enddo
      ! store tp1_pp(4,nlegs) > oext_tilde(nlegs_os_max,4)
      oext_tilde=0D0
      do i = 1, nlegs
       oext_tilde(i,:) = tp1_pp(:,i)
      enddo
      
      do i = 1, nlegs
! boost from partonic CMS to the LAB frame
      call boost0(rgam_boost, -rbetgam_boost,
     &  tp1_pp(1,i), tp1_pp(2,i), tp1_pp(3,i), tp1_pp(4,i))
      enddo
      ! store tp1_pp(4,nlegs) > oext_tilde_pp(nlegs_os_max,4)
      oext_tilde_pp=0D0
      do i = 1, nlegs
       oext_tilde_pp(i,:) = tp1_pp(:,i)
      enddo
      
! tpp(): boost from partonic CMS to the LAB frame
      tpp_pp = tpp
      call boost0(rgam_boost, -rbetgam_boost,
     &  tpp_pp(1), tpp_pp(2), tpp_pp(3), tpp_pp(4))
!!!]get os momenta in the Lab frame

      call Do_VecSet_tilde()
      call GetPolarVec(1,epsp,epspc) ! collinear photon at the first position

***[impose cuts
C the momenta here has to be the tilde kinematics
        cfail = 1 ! initialization
c        print*,"qext =",qext ! Ninh check
	call DoCut_parCM(cfail, qext(1,1))
	if ( cfail .eq. 0 ) then
         Gij = 0D0
         return
        endif
***]impose cuts

!!![diagonal term
	ctemp = 0
	do flags = 11, 22, 11
	call P1SquaredME(amp, tp, flags)
	ctemp = ctemp + amp(1)
	enddo
	
        Hij = DBLE(ctemp)*fac1
!!!]diagonal term

!!![(non-)diagonal term
        if (frame_polar.eq.1) then ! VV cms frame
         tqq = tpp_vv
        else if (frame_polar.eq.2) then ! lab frame
         tqq = tpp_pp
        else if (frame_polar.eq.3) then ! parent partonic frame
         tqq = tpp
        endif

	ctemp = 0
	do i1 = -1, 1, 2
	 do j1 = -1, 1, 2
         if (i1.eq.-1.and.j1.eq.-1) then
          flags = 11
         else if (i1.eq.-1.and.j1.eq.1) then
          flags = 12
         else if (i1.eq.1.and.j1.eq.-1) then
          flags = 21
         else if (i1.eq.1.and.j1.eq.1) then
          flags = 22
         endif
	call P1SquaredME(amp, tp, flags)

         ctemp = ctemp + amp(1)*
     &   (tqq(1)*epsp(i1,1) - tqq(2)*epsp(i1,2) 
     &   - tqq(3)*epsp(i1,3) - tqq(4)*epsp(i1,4))*
     &   (tqq(1)*epspc(j1,1) - tqq(2)*epspc(j1,2) 
     &   - tqq(3)*epspc(j1,3) - tqq(4)*epspc(j1,4))
	 enddo
	enddo
!!!](non-)diagonal term

        Hij = Hij + DBLE(ctemp)*fac1*fac2

c	if (abs(Hij).gt.0d0) then
c	print*,'Amp Hij =', Hij
c	stop
c	endif

        Gij = Hij
       IF (lscaleFR_dyn) then ! dynamical scales
        call scaleFR(.true.) ! choose true for subtraction terms
        if (lNLO_QCD) then
         call GetAlfas_dyn()
        endif
        fac_pdf = 1d0
        call MultiplyPdfs(fac_pdf)
       ELSE ! fixed scale
        if (.not.lpdf_fixed) then
         fac_pdf = 1d0
         call MultiplyPdfs(fac_pdf)
         lpdf_fixed = .true.
        endif
       ENDIF

      if (lNLO_QCD) then
       Gij = Alfas*Gij*fac_pdf
      else if (lNLO_EW) then
       Gij = Alfa*Gij*fac_pdf
      endif

        if (lHist_Fill) then
C the momenta here has to be the tilde kinematics
        call FillHIS(Gij, Jacb, qext)
        endif

       END
!!!]for photon-photon contribution of WW case
***********************************************************************
      subroutine DoTildeAmpHis_decayZ_OffS(tp,tp_off,iflag,Gij,Jacb)
C with on-shell projection after tilde (reduced) kinematics
        use global_vars
        use global_mom_xpdf
      implicit none
#include "model.h"
        integer iflag
	integer i1, j1, j2, i, j, nlegs, cfail
        parameter (nlegs = 8) ! for OS child q1 q2 > V3 V4 > l1 l2 l3 l4

	double precision Gij
	double precision result(2)
        double precision tp(2,4)
        double precision tp_off(2,4)
        double precision Lep(4,4)
        double precision MVint2(2) ! the off-shell invariant mass of the gauge bosons
	double precision Jacb

	double precision qext(proc_nlegs,4) ! important, array has to exactly match with FillHIS routines
	double precision tp1_pp(4,nlegs) ! in the Lab frame

! oext: q(1)q(2) > W(3)Z(4)gamma(5) > e(6)nu_e(7)mu+(8)mu-(9)
! pext: q(1)q(2) > e(3)nu_e(4)mu+(5)mu-(6)gamma(7)

***[for the on-shell gauge-boson decay particles
      if (lDECAYV2_NLO_EW) then
C first gauge boson (V1)
        Lep(:,1) = oext(6,:)
        Lep(:,2) = oext(7,:)
C second gauge boson (V2)
        Lep(:,3) = tp(1,:)
        Lep(:,4) = tp(2,:)
      else
C first gauge boson (V1)
        Lep(:,1) = tp(1,:)
        Lep(:,2) = tp(2,:)
C second gauge boson (V2)
        Lep(:,3) = oext(8,:)
        Lep(:,4) = oext(9,:)
      endif
***]for the on-shell gauge-boson decay particles

! store oext_tilde(nlegs_os_max,4): partonic CMS
      oext_tilde = 0D0
      do i = 1, 4
       oext_tilde(i,:) = oext(i,:)
       oext_tilde(i+4,:) = Lep(:,i)
      enddo

      oext_tilde_vv = oext_tilde
      
!!![get os momenta in the Lab frame
! for polarizations in the Lab frame:
! boost back to the Lab frame using (rgam_boost, rbetgam_boost)
      do i = 1, nlegs
       tp1_pp(:,i) = oext_tilde(i,:)
! boost from partonic CMS to the LAB frame
      call boost0(rgam_boost, -rbetgam_boost,
     &  tp1_pp(1,i), tp1_pp(2,i), tp1_pp(3,i), tp1_pp(4,i))
      enddo
      ! store tp1_pp(4,nlegs) > oext_tilde_pp(nlegs_os_max,4)
      oext_tilde_pp=0D0
      do i = 1, nlegs
       oext_tilde_pp(i,:) = tp1_pp(:,i)
      enddo
!!!]get os momenta in the Lab frame

      call Do_VecSet_tilde()
        
***[momenta for cut and hist in the partonic CMS in special format
! pext: q(1)q(2) > e(3)nu_e(4)mu+(5)mu-(6)gamma(7)
       do j1 = 1,2
        qext(j1,:) = pext(j1,:) ! off-shell mom
       enddo
      if (lDECAYV2_NLO_EW) then
        qext(3,:) = pext(3,:) ! off-shell mom
        qext(4,:) = pext(4,:) ! off-shell mom        
        qext(5,:) = tp_off(1,:)
        qext(6,:) = tp_off(2,:)
      else
        qext(3,:) = tp_off(1,:)
        qext(4,:) = tp_off(2,:)
        qext(5,:) = pext(5,:) ! off-shell mom
        qext(6,:) = pext(6,:) ! off-shell mom        
      endif
C This setting to zero is important for cut and histogram routines
        do j = 1, 4
        qext(proc_nlegs,j) = 0
        enddo

        ! store qext(proc_nlegs,4) > pext_tilde(nlegs_max,4)
        pext_tilde = 0D0
        do i = 1,proc_nlegs
         pext_tilde(i,:)=qext(i,:)
        enddo
***]momenta for cut and hist in the partonic CMS in special format

C Note: qext(6,4) are the off-shell tilde qq -> 4l momenta
	MVint2(1) = qext(3,1)*qext(4,1) - qext(3,2)*qext(4,2) 
     &  - qext(3,3)*qext(4,3) - qext(3,4)*qext(4,4)
        MVint2(1) = MVint2(1)*2
	MVint2(2) = qext(5,1)*qext(6,1) - qext(5,2)*qext(6,2) 
     &  - qext(5,3)*qext(6,3) - qext(5,4)*qext(6,4)
        MVint2(2) = MVint2(2)*2
!        print*,"DoTildeAmpHis: MVint2=",MVint2

***[impose cuts
C the momenta here has to be the tilde kinematics
        cfail = 1 ! initialization
	call DoCut_parCM(cfail, qext(1,1))
	if ( cfail .eq. 0 ) then
         Gij = 0D0
         return
        endif
***]impose cuts

	call P0SquaredME(result, MVint2(1), 1)
        Gij =  result(1)*Gij

       IF (lscaleFR_dyn) then ! dynamical scale
        call scaleFR(.true.) ! choose true for subtraction terms
        if (lNLO_QCD) then
         call GetAlfas_dyn()
        endif
        fac_pdf = 1d0
        call MultiplyPdfs(fac_pdf)
       ELSE ! fixed scale
        if (.not.lpdf_fixed) then
         fac_pdf = 1d0
         call MultiplyPdfs(fac_pdf)
         lpdf_fixed = .true.
        endif
       endif

      if (lNLO_QCD) then
       Gij = Alfas*Gij*fac_pdf
      else if (lNLO_EW) then
       Gij = Alfa*Gij*fac_pdf
      endif

        if (lHist_Fill) then
C the momenta here has to be the tilde kinematics
        call FillHIS(Gij, Jacb, qext) 
        endif

       END subroutine DoTildeAmpHis_decayZ_OffS
***********************************************************************
      subroutine DoTildeAmpHis_decayW_OffS(tp,tp_off,iflag,Gij,Jacb)
C with on-shell projection after tilde (reduced) kinematics
        use global_vars
        use global_mom_xpdf
      implicit none
#include "model.h"
        integer iflag
	integer i1, j1, j2, i, cfail
	integer nlegs
        parameter (nlegs = 8) ! for child q1 q2 > V3 V4 > l1 l2 l3 l4

	double precision Gij
	double precision result(2)
        double precision tp(2,4)
        double precision tp_off(2,4)
        double precision Lep(4,4)
        double precision MVint2(2) ! the off-shell invariant mass of the gauge bosons
	double precision Jacb

	double precision qext(proc_nlegs,4) ! important, array has to exactly match with FillHIS routines
	double precision tp1_pp(4,nlegs) ! in the Lab frame

! oext: q(1)q(2) > W(3)Z(4)gamma(5) > e(6)nu_e(7)mu+(8)mu-(9)
! pext: q(1)q(2) > e(3)nu_e(4)mu+(5)mu-(6)gamma(7), proc_nlegs = 7

***[for the on-shell gauge-boson decay particles
      if (lDECAYV1_NLO_EW) then
C first gauge boson
        Lep(:,1) = tp(1,:) ! electron
        Lep(:,2) = tp(2,:) ! neutrino
C second gauge boson
        Lep(:,3) = oext(8,:)
        Lep(:,4) = oext(9,:)
      else
! first gauge boson
        Lep(:,1) = oext(6,:)
        Lep(:,2) = oext(7,:)
! second gauge boson
        Lep(:,3) = tp(1,:) ! muon
        Lep(:,4) = tp(2,:) ! neutrino
      endif
***]for the on-shell gauge-boson decay particles

! store oext_tilde(nlegs_os_max,4): partonic CMS
      oext_tilde = 0D0
      do i = 1, 4
       oext_tilde(i,:) = oext(i,:)
       oext_tilde(i+4,:) = Lep(:,i)
      enddo

      oext_tilde_vv = oext_tilde
      
!!![get os momenta in the Lab frame
! for polarizations in the Lab frame:
! boost back to the Lab frame using (rgam_boost, rbetgam_boost)
      do i = 1, nlegs
       tp1_pp(:,i) = oext_tilde(i,:)
! boost from partonic CMS to the LAB frame
      call boost0(rgam_boost, -rbetgam_boost,
     &  tp1_pp(1,i), tp1_pp(2,i), tp1_pp(3,i), tp1_pp(4,i))
      enddo
      ! store tp1_pp(4,nlegs) > oext_tilde_pp(nlegs_os_max,4)
      oext_tilde_pp=0D0
      do i = 1, nlegs
       oext_tilde_pp(i,:) = tp1_pp(:,i)
      enddo
!!!]get os momenta in the Lab frame

      call Do_VecSet_tilde()

***[momenta for cut and hist in the partonic CMS in special format
! off-shell mom
! pext(nlegs_max,4): q(1)q(2) > e(3)nu_e(4)mu+(5)mu-(6)gamma(7)
       qext(1,:) = pext(1,:) 
       qext(2,:) = pext(2,:) 

      if (lDECAYV1_NLO_EW) then
       do j2 = 1,4
        qext(3,j2) = tp_off(1,j2) ! electron
        qext(4,j2) = tp_off(2,j2) ! neutrino
        qext(5,j2) = pext(5,j2) 
        qext(6,j2) = pext(6,j2) 
       enddo
      else
        qext(3,:) = pext(3,:)
        qext(4,:) = pext(4,:)
        qext(5,:) = tp_off(1,:) ! muon
        qext(6,:) = tp_off(2,:) ! neutrino
      endif
C This setting to zero is important for cut and histogram routines
        do j2 = 1, 4
        qext(proc_nlegs,j2) = 0
        enddo
        ! store qext(proc_nlegs,4) > pext_tilde(nlegs_max,4)
        pext_tilde = 0D0
        do i = 1,proc_nlegs
         pext_tilde(i,:)=qext(i,:)
        enddo
***]momenta for cut and hist in the partonic CMS in special format

!        print*,"E-p init=",qext(1,:)+qext(2,:)
!        print*,"E-p sum =",qext(1,:)+qext(2,:)-qext(3,:)-qext(4,:)-qext(5,:)-qext(6,:)

C Note: qext(6,4) are the off-shell tilde qq -> 4l momenta
	MVint2(1) = qext(3,1)*qext(4,1) - qext(3,2)*qext(4,2) 
     &  - qext(3,3)*qext(4,3) - qext(3,4)*qext(4,4)
        MVint2(1) = MVint2(1)*2
	MVint2(2) = qext(5,1)*qext(6,1) - qext(5,2)*qext(6,2) 
     &  - qext(5,3)*qext(6,3) - qext(5,4)*qext(6,4)
        MVint2(2) = MVint2(2)*2
!        print*,"DoTildeAmpHis: MVint2=",MVint2

***[impose cuts
C the momenta here has to be the tilde kinematics
        cfail = 1 ! initialization
	call DoCut_parCM(cfail, qext(1,1))
	if ( cfail .eq. 0 ) then
         Gij = 0D0
         return
        endif
***]impose cuts

	call P0SquaredME(result, MVint2(1), 1)
        Gij =  result(1)*Gij
        
       IF (lscaleFR_dyn) then ! dynamical scale
        call scaleFR(.true.) ! choose true for subtraction terms
        if (lNLO_QCD) then
         call GetAlfas_dyn()
        endif
        fac_pdf = 1d0
        call MultiplyPdfs(fac_pdf)
       ELSE ! fixed scale
        if (.not.lpdf_fixed) then
         fac_pdf = 1d0
         call MultiplyPdfs(fac_pdf)
         lpdf_fixed = .true.
        endif
       endif

      if (lNLO_QCD) then
       Gij = Alfas*Gij*fac_pdf
      else if (lNLO_EW) then
       Gij = Alfa*Gij*fac_pdf
      endif

        if (lHist_Fill) then
C the momenta here has to be the tilde kinematics
        call FillHIS(Gij, Jacb, qext) 
        endif

      End subroutine DoTildeAmpHis_decayW_OffS
***********************************************************************
      subroutine OSProject_2to2(Vext,rm2,rext,flag)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C input: Vext(0:3,2) are the momenta of 2 gauge bosons in 
C the rest frame of X = V1 + V2; rm2(2) are the physical 
C masses of the 2 gauge bosons. 
C output: rext(0:3,2) are the OS-projected momenta of the 
C two gauge bosons.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      logical flag
      double precision rext(0:3,2), Vext(0:3,2), rm2(2)
      double precision dot, a, b1, b2, b, Del, c
      double precision tmp1, tmp2, Vsum(0:3), SX, sqrtSX, temp(0:3)
      integer i

      flag = .true. ! initialization

      do i = 0, 3
       Vsum(i) = Vext(i,1)+Vext(i,2)
      enddo
      SX = Vsum(0)**2 - Vsum(1)**2 - Vsum(2)**2 - Vsum(3)**2
      sqrtSX = sqrt(SX)

c      print*,"sqrtSX=", sqrtSX
c      print*,"pW in=",Vext(0,1),Vext(1,1),Vext(2,1),Vext(3,1)
c      print*,"pZ in=",Vext(0,2),Vext(1,2),Vext(2,2),Vext(3,2)

***[first on-shell gauge boson
       tmp1 = (rm2(1) - rm2(2))/(2*sqrtSX) + sqrtSX/2 ! r10
       rext(0,1) = tmp1

       tmp2 = Vext(1,1)**2
     & + Vext(2,1)**2 + Vext(3,1)**2

       tmp2 = (rext(0,1)**2 - rm2(1))/tmp2
c      print*,"tmp2=", tmp2
       if ( tmp2 .lt. 0 ) then
        flag = .false.
        return
       endif
       tmp2 = sqrt(tmp2) ! b

       do i = 1, 3
       rext(i,1) = tmp2*Vext(i,1)
       enddo
***]first on-shell gauge boson

***[second on-shell gauge boson
       do i = 0, 3
       rext(i,2) = Vsum(i) - rext(i,1)
       enddo
***]second on-shell gauge boson
c      print*,"pW os=",rext(0,1),rext(1,1),rext(2,1),rext(3,1)
c      print*,"pZ os=",rext(0,2),rext(1,2),rext(2,2),rext(3,2)
      END
***********************************************************************
        subroutine Cal_mom_ID()
        use global_vars
        use global_mom_xpdf
        implicit none
        integer i

      vext_jend = 0d0 ! initialize all to zero

      IF (lNLO_FIN_IND) THEN
      ! move the jet from 3rd to the end
      do i = 1,2
       vext_jend(1:4,i) = pext(i,1:4)
       vext_jend(5,i) = leg(i)
      enddo
      do i = 3,proc_nlegs-1
       vext_jend(1:4,i) = pext(i+1,1:4)
       vext_jend(5,i) = leg(i+1)
      enddo
       vext_jend(1:4,proc_nlegs) = pext(3,1:4)
       vext_jend(5,proc_nlegs) = leg(3)
      ELSE
      ! no change in particle positions compared to pext()
      do i = 1,proc_nlegs
       vext_jend(1:4,i) = pext(i,1:4)
       vext_jend(5,i) = leg(i)
      enddo
      ENDIF

        end subroutine Cal_mom_ID
***********************************************************************
        subroutine translate_V2P(Vext,N)
        use global_vars
        use global_mom_xpdf
        implicit none
        integer N
        double precision Vext(0:3,N)
        integer i, j

!      Print*,"N =",N
!      do i = 1, N
!      print*,"Vext begin =",i,Vext(:,i),Vext(0,i)**2-Vext(1,i)**2-Vext(2,i)**2-Vext(3,i)**2
!      enddo

       pext = 0d0 ! initialize ALL to zero

***[ordering of particles (true kinematics)
      if (lNLO_FIN_IND) then
***[ for the induced process
C Note: pext(7,4) are the qq -> j + 4l momenta
       do i = 1,N
         do j = 1,4
          pext(i,j) = Vext(j-1,i)
         enddo
        enddo
***] for the induced process
      endif

      if (lNLO_FIN_RAD) then
***[ for the radiation process
C Note: pext(7,4) are the qq -> 4l + j momenta
       do i = 1,2
         do j = 1,4
          pext(i,j) = Vext(j-1,i)
         enddo
        enddo

       do i = 3,6
         do j = 1,4
          pext(i,j) = Vext(j-1,i+1)
         enddo
        enddo

       do j = 1,4
        pext(7,j) = Vext(j-1,3)
       enddo
***] for the radiation process
      endif
***]ordering of particles (true kinematics)

!          do i = 1, N
!          print*,"pext translated:",i,pext(i,:),pext(i,1)**2-pext(i,2)**2-pext(i,3)**2-pext(i,4)**2
!          enddo
!          print*,"check EP:",pext(1,:)+pext(2,:)-pext(3,:)-pext(4,:)-pext(5,:)-pext(6,:)-pext(7,:)


        end
***********************************************************************
        subroutine translate_P2V(Vext,N)
        use global_vars
        use global_mom_xpdf
        implicit none
        integer N
        double precision Vext(0:3,N)
        integer i, j

      if (lNLO_FIN_IND) then
***[ for the induced process
C Note: Vext(7,4) are the qq -> j + 4l momenta
       do i = 1,N
         do j = 1,4
          Vext(j-1,i) = pext(i,j)
         enddo
        enddo
***] for the induced process
      endif

      if (lNLO_FIN_RAD) then
***[ for the radiation process
C Note: pext(7,4) are the qq -> 4l + j momenta
C Note: Vext(7,4) are the qq -> j + 4l momenta
       do i = 1,2
         do j = 1,4
          Vext(j-1,i) = pext(i,j)
         enddo
        enddo

       do i = 3,6
         do j = 1,4
          Vext(j-1,i+1) = pext(i,j)
         enddo
        enddo

       do j = 1,4
        Vext(j-1,3) = pext(7,j)
       enddo
***] for the radiation process
      endif
        end
***********************************************************************
      subroutine OSProject(imap,Vext,rm2,rext)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C input: imap to choose mapping method,
C Vext(0:3,3:6) are the momenta of 4 leptons in 
C the rest frame of X = V1 + V2; rm2(2) are the physical 
C masses of the 2 gauge bosons. 
C output: rext(0:3,6) are the OS-projected momenta of the 
C two gauge bosons (first 2 entries) and of the 4 leptons, 
C in the rest frame of X.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer imap
      double precision rext(0:3,6), Vext(0:3,3:6), rm2(2)
      double precision tmp1, tmp2, Vsum(0:3), SX, sqrtSX
      integer i
      double precision rm(2)

      do i = 1, 2
       rm(i) = sqrt(rm2(i))
      enddo

      do i = 0, 3
       Vsum(i) = Vext(i,3)+Vext(i,4)+Vext(i,5)+Vext(i,6)
      enddo
      SX = Vsum(0)**2 - Vsum(1)**2 - Vsum(2)**2 - Vsum(3)**2
      sqrtSX = sqrt(SX)

!      print*,"SX, rm =", SX, rm(:)
!      do i = 3, 6
!      print*,"Vext i =",i, Vext(:,i), Vext(0,i)**2-Vext(1,i)**2-Vext(2,i)**2-Vext(3,i)**2
!      enddo

!      print*,"OSProject: sqrtSX child, imap=", sqrtSX, imap

***[first on-shell gauge boson
       tmp1 = (rm2(1) - rm2(2))/(2*sqrtSX) + sqrtSX/2 ! r10
       rext(0,1) = tmp1

       tmp2 = (Vext(1,3) + Vext(1,4))**2
     & + (Vext(2,3) + Vext(2,4))**2 + (Vext(3,3) + Vext(3,4))**2

       tmp2 = (rext(0,1)**2 - rm2(1))/tmp2
       tmp2 = sqrt(tmp2) ! b

c      print*,"from OSProject b=", tmp2, rext(0,1)**2, rm2(1)

       do i = 1, 3
       rext(i,1) = tmp2*(Vext(i,3) + Vext(i,4))
       enddo
***]first on-shell gauge boson

***[second on-shell gauge boson
       do i = 0, 3
       rext(i,2) = Vsum(i) - rext(i,1)
       enddo
***]second on-shell gauge boson

***[projected lepton momenta as Denner ZZ 2021
      if (imap .eq. 1) then
       call OSMap_decayV_LO(Vext(0,3),rm(1),rext(0,1),rext(0,3))
       call OSMap_decayV_LO(Vext(0,5),rm(2),rext(0,2),rext(0,5))
       return
      endif
***]projected lepton momenta as Denner ZZ 2021

***[projected lepton momenta as Denner RacoonWW 2000
      if (imap .eq. 2) then
       tmp1 = rext(0,1)*Vext(0,3) - rext(1,1)*Vext(1,3) 
     & - rext(2,1)*Vext(2,3) - rext(3,1)*Vext(3,3)
       tmp1 = rm2(1)/(2*tmp1)
       do i = 0, 3
        rext(i,3) = tmp1*Vext(i,3)
        rext(i,4) = rext(i,1) - rext(i,3)
       enddo

       tmp1 = rext(0,2)*Vext(0,5) - rext(1,2)*Vext(1,5)
     & - rext(2,2)*Vext(2,5) - rext(3,2)*Vext(3,5)
       tmp1 = rm2(2)/(2*tmp1)
       do i = 0, 3
        rext(i,5) = tmp1*Vext(i,5)
        rext(i,6) = rext(i,2) - rext(i,5)
       enddo
       return
      endif
***]projected lepton momenta as Denner RacoonWW 2000
      END subroutine OSProject
***********************************************************************
      subroutine OSMap_decayV_LO(p,MV,qV,q)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C input: p(0:3,1:2) are the momenta of 2 off-shell leptons, 
C MV is the physical mass of the gauge boson, 
C qV(0:3) is the OS-projected mom of the gauge boson 
C output: q(0:3,1:2) are the OS-projected momenta of the 
C two leptons in the original frame (where p are defined).
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      double precision p(0:3,2), q(0:3,2), MV, qV(0:3)
      double precision pV(0:3), pV_inv(0:3)
      double precision r(0:3) ! mom of the 1st particle in the off-shell V rest frame
      double precision nr(1:3) ! spatial directions of r
      double precision f(0:3,2) ! mom in the on-shell V rest frame
      integer i, j

!!![boost to the off-shell V-rest frame
      pV(:) = p(:,1) + p(:,2)
      pV_inv(0) = pV(0)
      do i = 1, 3
      pV_inv(i) = -pV(i)
      enddo
      call BOOSTN(p(0,1),pV_inv(0),r(0))
!!!]boost to the off-shell V-rest frame

!!![get the spatial directions
      if (r(0).gt.0D0) then
      do j = 1, 3
      nr(j) = r(j)/r(0) ! p_magnitude = Energy
      enddo
      else
      nr(:) = 0D0
      endif
!!!]get the spatial directions

!!![mom in the on-shell V-rest frame
      f(0,1) = MV/2
      f(0,2) = f(0,1)
      do i = 1, 3
      f(i,1) = nr(i)*f(0,1)
      f(i,2) = -f(i,1)
      enddo
!!!]mom in the on-shell V-rest frame

! boost from the on-shell V-rest frame to the original frame
      do i = 1, 2
      call BOOSTN(f(0,i),qV(0),q(0,i))
      enddo
      end subroutine OSMap_decayV_LO
***********************************************************************
      subroutine OSMap_decayV_NLO(p,MV,qV,q)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C input: p(0:3,1:3) are the momenta of 2 off-shell leptons + photon, 
C MV is the physical mass of the gauge boson, 
C qV(0:3) is the OS-projected mom of the gauge boson 
C output: q(0:3,1:3) are the OS-projected momenta of the 
C two leptons + photon in the original frame (where p are defined).
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      double precision p(0:3,3), q(0:3,3), MV, qV(0:3)
      double precision pV(0:3), pV_inv(0:3)
      double precision r(0:3,3) ! mom in the off-shell V rest frame
      double precision f(0:3,3) ! mom in the on-shell V rest frame
      double precision a, E(3), pV2
      double precision nr(1:3,3) ! spatial directions of r
      integer i, j

!!![boost to the off-shell V-rest frame
      pV(:) = p(:,1) + p(:,2) + p(:,3)
      pV_inv(0) = pV(0)
      do i = 1, 3
      pV_inv(i) = -pV(i)
      enddo
      do i = 1, 3
      call BOOSTN(p(0,i),pV_inv(0),r(0,i))
      enddo
!!!]boost to the off-shell V-rest frame

!!![get the spatial directions
      do i = 1, 3
      if (r(0,i).gt.0D0) then
      do j = 1, 3
      nr(j,i) = r(j,i)/r(0,i) ! p_magnitude = Energy 
      enddo
      else
      nr(:,i) = 0D0
      endif
      enddo
!!!]get the spatial directions

! rescale the energies
      pV2 = pV(0)**2-pV(1)**2-pV(2)**2-pV(3)**2
      a = MV/sqrt(pV2)
      do i = 1, 3
      f(0,i) = r(0,i)*a 
      enddo
!      print*,"sum s0 =",f(0,1)+f(0,2)+f(0,3),MV

!!![mom in the on-shell V-rest frame
      do j = 1, 3
      do i = 1, 3
      f(i,j) = nr(i,j)*f(0,j) ! p_magnitude = Energy
      enddo
      enddo
!!!]mom in the on-shell V-rest frame

! boost from the on-shell V-rest frame to the original frame
      do i = 1, 3
      call BOOSTN(f(0,i),qV(0),q(0,i))
      enddo
      end subroutine OSMap_decayV_NLO
***********************************************************************
      subroutine OSMap_AB(Vext,rm2,rext)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C input: Vext(0:3,2) are the momenta of A and B in 
C the rest frame of X = A + B; rm2(2) are the physical 
C masses of A and B. 
C output: rext(0:3,2) are the OS-projected momenta of the 
C two particles A and B in the rest frame of X.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      double precision rext(0:3,2), Vext(0:3,2), rm2(2)
      double precision tmp1, tmp2, Vsum(0:3), SX, sqrtSX
      integer i

      do i = 0, 3
       Vsum(i) = Vext(i,1)+Vext(i,2)
      enddo
      SX = Vsum(0)**2 - Vsum(1)**2 - Vsum(2)**2 - Vsum(3)**2
      sqrtSX = sqrt(SX)

c      print*,"sqrtSX child=", sqrtSX

***[first on-shell A
       tmp1 = (rm2(1) - rm2(2))/(2*sqrtSX) + sqrtSX/2 ! r10
       rext(0,1) = tmp1

       tmp2 = Vext(1,1)**2 + Vext(2,1)**2 + Vext(3,1)**2

       tmp2 = (rext(0,1)**2 - rm2(1))/tmp2
       tmp2 = sqrt(tmp2) ! b

       do i = 1, 3
       rext(i,1) = tmp2*Vext(i,1)
       enddo
***]first on-shell A

***[second on-shell B
       do i = 0, 3
       rext(i,2) = Vsum(i) - rext(i,1)
       enddo
***]second on-shell B
      END subroutine OSMap_AB
***********************************************************************
      subroutine OSMap_AB_frameX(fext,rm2,sext)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C input: fext(0:3,2) are the momenta of A and B in 
C the frame X; rm2(2) are the physical 
C masses of A and B. 
C output: sext(0:3,2) are the OS-projected momenta of the 
C two particles A and B in the frame X.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      double precision rext(0:3,2), Vext(0:3,2), rm2(2)
      double precision tmp1, tmp2, Vsum(0:3), SX, sqrtSX
      double precision Vsum_inv(0:3), Vsum_AB(0:3), sext(0:3,2), fext(0:3,2)
      double precision inv_mass_p
      integer i

      do i = 0, 3
       Vsum(i) = fext(i,1)+fext(i,2)
      enddo
      SX = Vsum(0)**2 - Vsum(1)**2 - Vsum(2)**2 - Vsum(3)**2
      sqrtSX = sqrt(SX)

c      print*,"sqrtSX child=", sqrtSX

!!![boost to the CMS_AB frame
      Vsum_inv(0) = Vsum(0)
      Vsum_inv(1:3) = -Vsum(1:3)
      do i = 1, 2
      call BOOSTN(fext(0,i),Vsum_inv(0),Vext(0,i))
      enddo
!!!]boost to the CMS_AB frame

      Vsum_AB(:) = Vext(:,1) + Vext(:,2)
!      print*,"Vsum_AB =",Vsum_AB(:)

***[first on-shell A
       tmp1 = (rm2(1) - rm2(2))/(2*sqrtSX) + sqrtSX/2 ! r10
       rext(0,1) = tmp1

       tmp2 = Vext(1,1)**2 + Vext(2,1)**2 + Vext(3,1)**2

       tmp2 = (rext(0,1)**2 - rm2(1))/tmp2
       tmp2 = sqrt(tmp2) ! b

       do i = 1, 3
       rext(i,1) = tmp2*Vext(i,1)
       enddo
***]first on-shell A

***[second on-shell B
       do i = 0, 3
       rext(i,2) = Vsum_AB(i) - rext(i,1)
       enddo
***]second on-shell B

!!![boost back to the X frame
      do i = 1, 2
      call BOOSTN(rext(0,i),Vsum(0),sext(0,i))
      enddo
!!!]boost back to the X frame

      END subroutine OSMap_AB_frameX
***********************************************************************
      subroutine OSMap_decayV_AB(p,MV,qV,q)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C decay: V -> AB;
C input: p(0:3,1:2) are the momenta of 2 off-shell decay products A and B; 
C MV(3) is the mass of the gauge boson V, A, and B; 
C qV(0:3) is the mom of V;
C output: q(0:3,1:2) are the OS-projected momenta of A and B 
C in the original frame (where p are defined).
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      double precision p(0:3,2), q(0:3,2), MV(3), qV(0:3)
      double precision pV(0:3), pV_inv(0:3)
      double precision r(0:3) ! mom of the 1st particle in the off-shell V rest frame
      double precision f(0:3,2) ! mom in the on-shell V rest frame
      double precision nr(1:3) ! normalized spatial direction
      double precision pp1 ! magnitude of the spatial vector
      integer i, j

!!![boost to the off-shell V-rest frame
      pV(:) = p(:,1) + p(:,2)
      pV_inv(0) = pV(0)
      do i = 1, 3
      pV_inv(i) = -pV(i)
      enddo
      call BOOSTN(p(0,1),pV_inv(0),r(0))
!!!]boost to the off-shell V-rest frame

!!![get the spatial directions
      pp1 = r(1)**2 + r(2)**2 + r(3)**2
      pp1 = sqrt(pp1)
      if (pp1.gt.0D0) then
      do j = 1, 3
      nr(j) = r(j)/pp1 ! normalized direction
      enddo
      else
      nr(:) = 0D0
      endif
!!!]get the spatial directions

!!![mom in the on-shell V-rest frame
      f(0,1) = 0.5D0*((MV(2)**2-MV(3)**2)/MV(1) + MV(1))
      f(0,2) = MV(1) - f(0,1)
      pp1 = f(0,1)**2 - MV(2)**2
! at this step here, the pp1 should be positive because 
! the on-shell condition of MV > (MA + MB) must be fullfiled:
      if (pp1 .ge. 0d0) then
       pp1 = sqrt(pp1)
      else
       if (abs(pp1)/MV(2)**2 .lt. 1D-10) then
        pp1 = 0d0
       else
        print*,"Error: something wrong with on-shell mapping, pp1 =",pp1
        stop
       endif
      endif
      do i = 1, 3
      f(i,1) = nr(i)*pp1
      f(i,2) = -f(i,1)
      enddo
!!!]mom in the on-shell V-rest frame

! boost from the on-shell V-rest frame to the original frame
      do i = 1, 2
      call BOOSTN(f(0,i),qV(0),q(0,i))
      enddo

      end subroutine OSMap_decayV_AB
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_amp2_hel_VVV(amp, result)
      use global_vars
      implicit none
      double complex amp(3,3,3), amp_all, amp_TTT, amp_TTL, amp_TLL, amp_LLL
      double precision result, amp2_all, amp2_TTT, amp2_TTL, amp2_TLL, amp2_LLL
      integer Hel3, Hel4, Hel5
      integer HelProd

      amp_all = (0d0, 0d0)
      amp_TTT = (0d0, 0d0)
      amp_TTL = (0d0, 0d0)
      amp_TLL = (0d0, 0d0)
      amp_LLL = (0d0, 0d0)
      do Hel3 = 1, 3
      do Hel4 = 1, 3
      do Hel5 = 1, 3
       amp_all = amp_all + amp(Hel3, Hel4, Hel5)
       HelProd = Hel3*Hel4*Hel5
       if (mod(HelProd,2).ne.0) then
        amp_TTT = amp_TTT + amp(Hel3, Hel4, Hel5)
       else if (mod(HelProd,4).ne.0) then
        amp_TTL = amp_TTL + amp(Hel3, Hel4, Hel5)
       else if (mod(HelProd,8).ne.0) then
        amp_TLL = amp_TLL + amp(Hel3, Hel4, Hel5)
       else
        amp_LLL = amp_LLL + amp(Hel3,Hel4,Hel5)
       endif
      enddo
      enddo
      enddo

      amp2_all = DBLE(DCONJG(amp_all)*amp_all)
      amp2_TTT = DBLE(DCONJG(amp_TTT)*amp_TTT)
      amp2_LLL = DBLE(DCONJG(amp_LLL)*amp_LLL)
      amp2_TTL = DBLE(DCONJG(amp_TTL)*amp_TTL)
      amp2_TLL = DBLE(DCONJG(amp_TLL)*amp_TLL)

      result = amp2_all - amp2_TTT - amp2_LLL - amp2_TTL - amp2_TLL

      end subroutine get_amp2_hel_VVV
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine check_mom(p, n, lpass)
      use, intrinsic :: ieee_arithmetic
      implicit none
      double precision p(0:3,n), q_all(0:3), q_sum
      integer n, i
      logical lpass
      lpass = .true.
      q_all(:) = p(:,1) + p(:,2) 
      do i = 3, n
       q_all(:) = q_all(:) - p(:,i) 
      enddo
      q_sum = q_all(0) + q_all(1) + q_all(2) + q_all(3)
      if (ieee_is_nan(q_sum) .or. .not. ieee_is_finite(q_sum)) then
       lpass = .false.
       return
      endif
      end subroutine check_mom
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine OSmap_solve_alpha2(n_res,mQ,m_sq,kk_sq,alpha2_k)
! target: find alpha2_k satisfying this equation
! sum_{i,1,n_res} sqrt{m_sq(i)^2 + alpha2_k*kk_sq(i)^2} - mQ = 0
! method: Newton-Raphson
      implicit none
      integer n_res, i
      integer :: i_max = 1000
      double precision :: eps_prec = 1D-14
      double precision mQ, alpha2_k, df, f, b
      double precision m_sq(1:n_res), kk_sq(1:n_res)
      ! first guess
      alpha2_k = 1D0
      do i = 1, i_max
      ! calculate f and its derivative:
      call OSmap_df_alpha2(n_res,mQ,m_sq,kk_sq,alpha2_k,df,f)
      ! check the value of f:
      if (abs(f/mQ).le.eps_prec) then
!       print*,"f/mQ, alpha2_k, i, eps_prec =", f/mQ, alpha2_k, i, eps_prec
       exit
      else
      ! calculate b:
      b = f - df*alpha2_k
      ! update alpha2_k
      alpha2_k = -b/df
      endif
      enddo

      end subroutine OSmap_solve_alpha2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine OSmap_df_alpha2(n_res,mQ,m_sq,kk_sq,alpha2_k,df,f)
! target: calculate f and its derivative wrt alpha2_k at alpha2_k, with
! f(alpha2_k) = sum_{i,1,n_res} sqrt{m_sq(i)^2 + alpha2_k*kk_sq(i)^2} - mQ
      implicit none
      integer n_res, i
      double precision mQ, alpha2_k, df, f, temp
      double precision m_sq(1:n_res), kk_sq(1:n_res)
      df = 0d0
      f = -mQ
      do i = 1, n_res
      temp = sqrt(m_sq(i) + kk_sq(i)*alpha2_k)
      df = df + kk_sq(i)/temp
      f = f + temp
      enddo
      df = df*0.5D0
      end subroutine OSmap_df_alpha2



