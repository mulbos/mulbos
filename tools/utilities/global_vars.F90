!*************************************************************************  
! This file helps to calculate global parameters.
!*************************************************************************  

module global_vars
      integer, parameter :: nlegs_max=9 ! TriBoson maximally (NLO)
      integer, parameter :: nlegs_os_max=12 ! added 3 bosons
      integer, parameter :: nlegs_res_max=3 ! 3 resonances maximally

      integer, parameter :: n_id_max=100 ! max number of process ids
      integer, parameter :: n_cont_max=30 ! max number of contributions per process

      integer, parameter :: n_scales_max=100 ! max number of scale choices

      integer :: nlegs_res ! number of resonances
      integer leg(nlegs_max) ! particle IDs
      double precision charge_qed(nlegs_max), charge_qcd(nlegs_max)
      integer proc_id, proc_cont, proc_nlegs, proc_nlegs_os
      !integer n_id ! actual number of idp (process id)
      integer n_cont(n_id_max) ! the maximal idc for each idp (contribution id, starting from zero)
      integer Nlepton, Nlep_charged, Nquark, Nphoton, Ngluon
      integer Nup, Ndown, Nbottom
      character*200 dir_name
      character*512 dir_save
      character*200 pdf_set
      integer seed, pdf_mem, N_stat
      integer polar_case, osmap_case, dist_case
      double precision ecm
      double precision coupling_GF
      double precision mass_W, mass_Z, mass_Top, mass_Bot, mass_H 
      double precision width_W, width_Z, width_Top, width_H
      double precision mass_Wos, mass_Zos
      double precision width_Wos, width_Zos
      double precision mass_VV(nlegs_res_max), width_VV(nlegs_res_max)
      integer scale_muf_case, scale_mur_case
      double precision muf_user, mur_user
      double precision scale_def ! the default scale
      double precision scale ! factorization scale
      double precision scaleren ! renormalization scale
      double precision scale_facF ! pre-factor for factorization scale  
      double precision scale_facR ! pre-factor for renormalization scale
      logical lscaleF_dyn, lscaleR_dyn, lscaleFR_dyn ! switches for dynamical scales
      logical lscaleFR_fixed ! switches for fixed scales
      double precision fac_identical, fac_color, fac_pdf
      logical lpdf_fixed
      
      integer hel_W, hel_Z
      integer hel_W_min, hel_W_max, hel_Z_min, hel_Z_max
      integer hel_W_step, hel_Z_step
      integer hel_V1, hel_V2, hel_V3
      integer hel_V1_step, hel_V2_step, hel_V3_step
      integer hel_V1_min, hel_V2_min, hel_V3_min
      integer hel_V1_max, hel_V2_max, hel_V3_max
      integer frame_polar
      
      integer NCOMPOmin, NCOMPO
      
      integer ITMX1in, ITMX2in, NWILDin, NWILD_born
      double precision ACC1in, ACC2in
      logical lWW, lZZ, lWZ
      logical lWWW
      logical lDPA, lFULL_CAL
      logical lLHAPDF
      logical lTREE, lVIRT, lNLO_COLL_RAD, lNLO_COLL_IND, lNLO_FIN_RAD, lNLO_FIN_IND
      logical lNLO_FIN_IND_TW
      logical lNLO_COLL, lNLO_FIN, lGLUGLU, lBOTBOT, lGAMGAM 
      logical lNLO_EW, lNLO_QCD, lDECAYV1_NLO_EW, lDECAYV2_NLO_EW, lDECAYV3_NLO_EW
      logical lMassReg, lDimReg
      logical lPDF_SCHEME_DIS, lPDF_SCHEME_EW_DIS
      integer input_scheme
      logical lHist_Fill
      logical lcomplex_mass_scheme, lFLV5
      logical lFC9 ! whether to use FormCalc-9X (X=8,10); FC-6.0 otherwise
      logical lBproc4FS, lTTbar_OS_4FS
      logical lDiBoson, lTriBoson
      logical lTTL_combined, lTLL_combined ! for combination of polarizations
      logical lPOLAR_INTERF ! for polarization interference
      
!!![technical cuts
! normalized 4-point Gram determinant cut: possible values 1D-4, 1D-5 (default), 1D-6
      double precision, parameter :: eps_GramDet4=1D-5
!!!]technical cuts

contains

      subroutine Read_Input(path_file)
      use readinput
      implicit none
      character*(*) path_file

        call loadfile(path_file)
        call read_int("PROC_ID",proc_id,1)
        call read_int("PROC_CONT",proc_cont,1)
        
        call read_string("OUTPUT_DIR",dir_name,"dummy")

        call read_int("HIST_WRITEOUT",dist_case,2)
        
        call read_int("RANDOM_SEED",seed,1)
        call read_int("RANDOM_POINTS",N_stat,20)
        
        call read_string("PDF_SET",pdf_set,"LUXqed17_plus_PDF4LHC15_nnlo_30.LHgrid")
        call read_int("PDF_MEMBER",pdf_mem,0)
        
        call read_int("OS_MAP",osmap_case,1)
        call read_int("POLARIZATION",polar_case,444)
        call read_int("FRAME_POLAR",frame_polar,1)
        
        call read_Real("ECM",ecm,13000d0)
        
        call read_int("SCALE_MUF",scale_muf_case,1)
        call read_int("SCALE_MUR",scale_mur_case,1)
        call read_Real("MUF_USER",muf_user,90D0)
        call read_Real("MUR_USER",mur_user,90D0)
        call read_Real("MUF_FAC",scale_facF,1D0)
        call read_Real("MUR_FAC",scale_facR,1D0)
        
        call read_logical("PDF_SCHEME_EW_DIS",lPDF_SCHEME_EW_DIS,.true.)
        
        call read_Real("COUPLING_GF",coupling_GF,1.16637D-5)
        
        call read_Real("WMASS",mass_Wos,80.385D0)
        call read_Real("ZMASS",mass_Zos,91.1876D0)
        call read_Real("TOPMASS",mass_Top,173.0D0)
        call read_Real("BOTMASS",mass_Bot,4.7D0)
        call read_Real("HMASS",mass_H,125D0)
        
        call read_Real("WWIDTH",width_Wos,2.085D0)
        call read_Real("ZWIDTH",width_Zos,2.4952D0)
        call read_Real("TOPWIDTH",width_Top,1.42D0)
        call read_Real("HWIDTH",width_H,0.0041D0)
        
        call closefile
      end subroutine
!**********************************************************************
      subroutine set_run_params()
! this code is to set important information for the run
      implicit none
      integer i
      integer :: processes(nlegs_max,0:n_cont_max,n_id_max)

      call load_data_processes(processes)

     lDiBoson = .false.
     lTriBoson = .false.
     mass_VV = 0d0
     width_VV = 0d0
     
     lWW = .false.      
     lZZ = .false.
     lWZ = .false.
     lWWW = .false.
     lFLV5 = .false.
     lFC9 = .false.

     if (proc_cont.eq.0) then
      lcomplex_mass_scheme = .true.
     else
      lcomplex_mass_scheme = .false.
     endif
     
     if (lcomplex_mass_scheme) then
     mass_Z = mass_Zos/sqrt(1+(width_Zos/mass_Zos)**2) 
     mass_W = mass_Wos/sqrt(1+(width_Wos/mass_Wos)**2)
     width_Z = width_Zos/sqrt(1+(width_Zos/mass_Zos)**2)
     width_W = width_Wos/sqrt(1+(width_Wos/mass_Wos)**2)
     else
     mass_Z = mass_Zos
     mass_W = mass_Wos
     width_Z = width_Zos
     width_W = width_Wos
     endif
     
     if (proc_id.eq.1.or.proc_id.eq.2) then
       mass_VV(1)=mass_W
       mass_VV(2)=mass_Z
       width_VV(1)=width_W
       width_VV(2)=width_Z
       lWZ=.true.
     else if ((proc_id-3)*(proc_id-7).le.0.or.(proc_id-13)*(proc_id-16).le.0) then
       mass_VV(1)=mass_W
       mass_VV(2)=mass_W
       width_VV(1)=width_W
       width_VV(2)=width_W
       lWW=.true.
     else if ((proc_id-8)*(proc_id-12).le.0) then! ZZ from 8 to 12
       mass_VV(1)=mass_Z
       mass_VV(2)=mass_Z
       width_VV(1)=width_Z
       width_VV(2)=width_Z
       lZZ=.true.
       lFLV5=.false. ! can be true for a fast and (approximate) results, not recommended
     else if ((proc_id-51)*(proc_id-52).le.0) then ! WWW from 51 to 52
       mass_VV(1)=mass_W
       mass_VV(2)=mass_W
       mass_VV(3)=mass_W
       width_VV(1)=width_W
       width_VV(2)=width_W
       width_VV(3)=width_W
       lWWW=.true.
     else
       print*,"invalid process_id:",proc_id
       stop
     endif

     lGLUGLU = .false.
     lBOTBOT = .false.
     lGAMGAM = .false.
     lBproc4FS = .false.
     lTTbar_OS_4FS = .false.
     if (proc_id.eq.5.or.proc_id.eq.10) then
      lGLUGLU = .true.
     else if (proc_id.eq.6.or.proc_id.eq.11) then
      lBOTBOT = .true.
     else if (proc_id.eq.7.or.proc_id.eq.12) then
      lGAMGAM = .true.
     else if (proc_id.ge.13.and.proc_id.le.16) then
      lBproc4FS = .true.
     endif
     
      leg(:) = processes(:,proc_cont,proc_id)
      proc_nlegs = 0
      Nlepton = 0
      Nquark = 0
      Nup = 0
      Ndown = 0
      Nbottom = 0
      Nphoton = 0
      Ngluon = 0
      Nlep_charged = 0
      do i = 1, nlegs_max
       if (leg(i).ne.0) then
        proc_nlegs = proc_nlegs + 1
        if (abs(leg(i)).ge.11.and.abs(leg(i)).le.16) then
         Nlepton = Nlepton + 1
         charge_qed(i) = -mod(leg(i),2)
         charge_qcd(i) = 0
         if (mod(abs(leg(i)),2).eq.1) then
          Nlep_charged = Nlep_charged + 1
         endif
        else if (abs(leg(i)).ge.1.and.abs(leg(i)).le.6) then
         Nquark = Nquark + 1
         charge_qcd(i) = leg(i)/abs(leg(i))         
         if (mod(leg(i),2).eq.0) then
          Nup = Nup + 1
          charge_qed(i) = 2D0/3*leg(i)/abs(leg(i))
         else
          Ndown = Ndown + 1
          charge_qed(i) = -1D0/3*leg(i)/abs(leg(i))
          if (abs(leg(i)) .eq. 5) Nbottom = Nbottom + 1
         endif
        else if (leg(i).eq.21) then
         Ngluon = Ngluon + 1
         charge_qcd(i) = 0         
        else if (leg(i).eq.22) then
         Nphoton = Nphoton + 1
         charge_qcd(i) = 0         
        endif
       endif
      enddo

      fac_identical = 1D0
      fac_color = 1D0
      fac_pdf = 1D0
      lpdf_fixed = .false.
      do i=1,2
      if (abs(leg(i)).ge.1.and.abs(leg(i)).le.6) then ! a quark
       fac_color = fac_color/3D0
      else if (leg(i).eq.21) then ! a gluon
       fac_color = fac_color/8D0
      endif
      enddo
      
      if (Nlepton .eq. 4) then
       nlegs_res = 2
       lDiBoson = .true.
      else if (Nlepton .eq. 6) then
       nlegs_res = 3
       lTriBoson = .true.
      else
       print*,"invalid process Nlepton =",Nlepton
       stop
      endif
      proc_nlegs_os = proc_nlegs + nlegs_res ! add intermediate bosons
      print*,"Nlepton =",Nlepton
      print*,"Nlep_charged =",Nlep_charged
      
      lDPA = .true.
      
      lFULL_CAL = .false.
      lTREE = .false.
      lVIRT = .false.
      lNLO_COLL_RAD = .false.
      lNLO_COLL_IND = .false.
      lNLO_FIN_RAD = .false.
      lNLO_FIN_IND = .false.
      lNLO_FIN_IND_TW = .false.
      lNLO_EW = .false.
      lNLO_QCD = .false.
      lDECAYV1_NLO_EW = .false.
      lDECAYV2_NLO_EW = .false.
      lDECAYV3_NLO_EW = .false.
      
      if (proc_cont.eq.0) then
       lFULL_CAL = .true.
       lDPA = .false.
       lTREE = .true.
      else if (proc_cont.eq.1) then
       lTREE = .true.
      else if (proc_cont.eq.2.or.proc_cont.eq.8) then
       lVIRT = .true.
      else if (proc_cont.eq.3.or.proc_cont.eq.9) then
       lNLO_COLL_RAD = .true.
      else if (proc_cont.eq.4.or.proc_cont.eq.10) then
       lNLO_COLL_IND = .true.
      else if (proc_cont.eq.5.or.proc_cont.eq.11) then
       lNLO_FIN_RAD = .true.
      else if (proc_cont.eq.6.or.proc_cont.eq.12) then
       lNLO_FIN_IND = .true.      
      else if (proc_cont.eq.7.or.proc_cont.eq.13) then
       lNLO_FIN_IND = .true.
      endif
      
      IF (lDiBoson) then
      if (proc_cont.eq.14) then
       lVIRT = .true.
       lDECAYV2_NLO_EW = .true.
      else if (proc_cont.eq.15) then
       lNLO_FIN_RAD = .true.
       lDECAYV2_NLO_EW = .true.
      else if (proc_cont.eq.16) then
       lVIRT = .true.
       lDECAYV1_NLO_EW = .true.
      else if (proc_cont.eq.17) then
       lNLO_FIN_RAD = .true.
       lDECAYV1_NLO_EW = .true.
      endif
      if (lWW) then
      if (proc_cont.ge.18.and.proc_cont.le.21) then
       lNLO_FIN_IND = .true.
       lNLO_FIN_IND_TW = .true.
      endif
      endif
      ENDIF
      
      IF (lTriBoson) then
      if (proc_cont.eq.14) then
       lVIRT = .true.
       lDECAYV3_NLO_EW = .true.
      else if (proc_cont.eq.15) then
       lNLO_FIN_RAD = .true.
       lDECAYV3_NLO_EW = .true.
      else if (proc_cont.eq.16) then
       lVIRT = .true.
       lDECAYV2_NLO_EW = .true.
      else if (proc_cont.eq.17) then
       lNLO_FIN_RAD = .true.
       lDECAYV2_NLO_EW = .true.
      else if (proc_cont.eq.18) then
       lVIRT = .true.
       lDECAYV1_NLO_EW = .true.
      else if (proc_cont.eq.19) then
       lNLO_FIN_RAD = .true.
       lDECAYV1_NLO_EW = .true.
      endif
      ENDIF
      
      if (lDiBoson) then
      if (proc_cont.ge.2.and.proc_cont.le.7) then
       lNLO_QCD = .true.
      else if (proc_cont.ge.8.and.proc_cont.le.17) then
       lNLO_EW = .true.
      endif
      if (lWW) then
      if (proc_cont.ge.18.and.proc_cont.le.19) then
       lNLO_QCD = .true.
      else if (proc_cont.ge.20.and.proc_cont.le.21) then
       lNLO_EW = .true.
      endif
      endif
      endif
      
      if (lTriBoson) then
      if (proc_cont.ge.2.and.proc_cont.le.7) then
       lNLO_QCD = .true.
      else if (proc_cont.ge.8.and.proc_cont.le.19) then
       lNLO_EW = .true.
      endif
      endif

!      if (lTriBoson) then
!       if (osmap_case .lt. 51) then 
!        osmap_case = 51 ! default choice
!       endif
!      endif
      
!      if (lDiBoson) then
!       if (osmap_case .gt. 40) then
!        osmap_case = 1 ! default choice
!       endif
!      endif
      
      lNLO_COLL = lNLO_COLL_RAD.or.lNLO_COLL_IND
      lNLO_FIN = lNLO_FIN_RAD.or.lNLO_FIN_IND
      if (lNLO_EW) then
       lPDF_SCHEME_DIS = lPDF_SCHEME_EW_DIS
      else
       lPDF_SCHEME_DIS = .false.
      endif
            
      ! Alpha_Gmu = 1 (recommended), Alpha(MZ) = 2, Alpha(0) = 3
      input_scheme = 1
      lLHAPDF = .true.
      lMassReg = .true. 
      lDimReg = .not.lMassReg
      
      if (lVIRT) then
       NCOMPOmin = 2
      else
       NCOMPOmin = 1
      endif
      NCOMPO = NCOMPOmin

! Bases parameters:      
      ITMX1in = 5
      ITMX2in = 5
      NWILD_born = (proc_nlegs-2)*3 + 2 - 5 ! remove the trivial-phi angle
      if (lTREE) then
       NWILDin = NWILD_born
       ACC1in = 0.002D0
       ACC2in = 0.001D0
      else if (lVIRT) then
       NWILDin = NWILD_born
       ACC1in = 0.02D0
       ACC2in = 0.01D0
      else if (lNLO_COLL) then
       NWILDin = NWILD_born + 1
       ACC1in = 0.02D0
       ACC2in = 0.01D0
      else if (lNLO_FIN) then
       NWILDin = NWILD_born
       ACC1in = 0.02D0
       ACC2in = 0.01D0
      else
       print*,"no valid component"
       stop
      endif
      
      if (NWILDin .gt. 15) then
       NWILDin = 15
      endif

      IF (dist_case .le. 3) THEN
        lHist_Fill = .true.
        ELSE IF (dist_case .eq. 4) THEN        
        lHist_Fill = .false.
        ELSE
        print*,"Invalid dist_case: dist_case =", dist_case
        stop
      ENDIF

      lscaleF_dyn = .false.
      lscaleR_dyn = .false.
      scale_def = 0d0
      do i = 1, nlegs_res
      scale_def = scale_def + mass_VV(i)
      enddo
      scale_def = scale_def/2d0 ! default choice

      scale = scale_def ! factorization scale
      scaleren = scale ! renormalization scale
      
! factorization scale:
      if (scale_muf_case.eq.1) then
       scale = scale_def
      else if (scale_muf_case.eq.0) then
       scale = muf_user
      else
       lscaleF_dyn = .true.
      endif

! renormalization scale:      
      if (scale_mur_case.eq.1) then
       scaleren = scale_def ! central scale
      else if (scale_mur_case.eq.0) then
       scaleren = mur_user
      else
       lscaleR_dyn = .true.
      endif
      lscaleFR_dyn = lscaleF_dyn.or.lscaleR_dyn
      lscaleFR_fixed = .not.lscaleFR_dyn
      
!!![fixed scales known here
      if (.not. lscaleF_dyn) then
       scale = scale_facF*scale ! factorization scale (appears in pdf)
      endif
      if (.not. lscaleR_dyn) then
       scaleren = scale_facR*scaleren ! renormalization scale
      endif
!!!]fixed scales known here
        
!***[polarization states
      IF (lDPA) THEN
! set initial values:
      hel_V1_min = -1
      hel_V1_max = 1
      hel_V1_step = 1

      hel_V2_min = -1
      hel_V2_max = 1
      hel_V2_step = 1

      hel_V3_min = -1
      hel_V3_max = 1
      hel_V3_step = 1

      lTTL_combined = .false.
      lTLL_combined = .false.
      lPOLAR_INTERF = .false.

! hel_V = 1 (transverse left), 2 (longitudinal), 3 (transverse right), 4 (unpolarized), 
! 5 (transverse left and right and L-R interference)
      if (lDiBoson) then
       if (polar_case.eq.444) then
        polar_case = 44 ! change default value for DiBoson case
       endif
      endif
      
      if (polar_case.lt.0) then
      ! for combined polarization cases:
       if (polar_case.eq.-552) then
        lTTL_combined = .true.
       else if (polar_case.eq.-522) then
        lTLL_combined = .true.
       else if (polar_case.eq.-444) then
        lPOLAR_INTERF = .true.
       endif
      else
      ! update hel ranges and step:
      call get_polar_indices()
      endif

      ENDIF
!***]polarization states

      end subroutine
!**********************************************************************
      subroutine get_polar_indices()
! this code is to set important information for the run
      implicit none
      integer i

      if (lTriBoson) then
!!![triboson
       hel_V1 = floor(polar_case/100d0)
       hel_V2 = floor((polar_case-hel_V1*100)/10d0)
       hel_V3 = floor((polar_case-hel_V1*100-hel_V2*10)/1d0)
       call get_hel_bounds(hel_V1,hel_V1_step,hel_V1_min,hel_V1_max)
       call get_hel_bounds(hel_V2,hel_V2_step,hel_V2_min,hel_V2_max)
       call get_hel_bounds(hel_V3,hel_V3_step,hel_V3_min,hel_V3_max)
!!!]triboson
      endif
!!!!!!!!!!!!!!!!!!!!!!!!      
      if (lDiBoson) then
!!![diboson
       hel_V1 = floor(polar_case/10d0)
       hel_V2 = floor((polar_case-hel_V1*10)/1d0)
       call get_hel_bounds(hel_V1,hel_V1_step,hel_V1_min,hel_V1_max)
       call get_hel_bounds(hel_V2,hel_V2_step,hel_V2_min,hel_V2_max)
       hel_W = hel_V1
       hel_Z = hel_V2
       hel_W_step = hel_V1_step
       hel_W_min = hel_V1_min
       hel_W_max = hel_V1_max
       hel_Z_step = hel_V2_step
       hel_Z_min = hel_V2_min
       hel_Z_max = hel_V2_max
!!!]diboson
      return
      endif
      end subroutine get_polar_indices
!**********************************************************************
      subroutine get_hel_bounds(hel_V,hel_V_step,hel_V_min,hel_V_max)
! this code is to set important information for the run
      implicit none
      integer hel_V,hel_V_step,hel_V_min,hel_V_max 

      hel_V_step = 1 ! default value

        if (hel_V .eq. 1) then
          hel_V_min = -1
          hel_V_max = hel_V_min
        elseif (hel_V .eq. 2) then
          hel_V_min = 0
          hel_V_max = hel_V_min
        elseif (hel_V .eq. 3) then
          hel_V_min = 1
          hel_V_max = hel_V_min
        elseif (hel_V .eq. 4) then
          hel_V_min = -1
          hel_V_max = 1
        elseif (hel_V .eq. 5) then
          hel_V_min = -1
          hel_V_max = 1
          hel_V_step = 2
        else
         print*,"invalid input for V helicity: hel_V=", hel_V
         stop
        endif

      end subroutine get_hel_bounds
!**********************************************************************
      subroutine get_sm_para_ext()
! this code is to pass information to old-style common block.
! do not change this code.      
        implicit none
        double precision GF
        double precision MW, MZ, MT, MH, MBot 
        double precision GamW, GamZ, GamTq, GamH
        ! this common block has to match with the one in model.h
        common /sm_para_ext/ GF, MZ, MW, MT, MH, GamZ, GamW, GamTq, GamH, MBot
        GF = coupling_GF
        MW = mass_W
        MZ = mass_Z
        MT = mass_Top
        MBot = mass_Bot
        MH = mass_H
        GamW = width_W
        GamZ = width_Z
        GamTq = width_Top
        GamH = width_H
      end subroutine
!**********************************************************************
      subroutine load_data_processes(processes)
! this code is to define the list of processes and its components
      implicit none
      integer i, j, i_proc, i_cont
!!![List of processes
! use zero (0) for no particle
! follow MC particle numbering scheme      
      integer :: processes(nlegs_max,0:n_cont_max,n_id_max)
!!!]List of processes

      processes(:,:,:) = 0 ! intial values
!!![process 1: WpZ
      ! born and NLO QCD corrections:
      ! UDbarWpZ-born-Hahn
      do i = 0, 4 ! full-born (0), DPA-born(1), virt(2), coll_rad (3), coll_ind(4)
      processes(1,i,1) = 2 ! up
      processes(2,i,1) = -1 ! anti-down
      processes(3,i,1) = -11 ! e+
      processes(4,i,1) = 12 ! nu_e
      processes(5,i,1) = -13 ! mu+
      processes(6,i,1) = 13 ! mu-
      enddo
      ! qext(LEGS,4), qext_os(LEGS,4), pext(LEGS,4)
      ! dir[5]=UDbarWpZglu-fin-Hahn
      do i = 5, 5 ! fin_rad(5): gluon at the end
      processes(1,i,1) = 2
      processes(2,i,1) = -1
      processes(3,i,1) = -11
      processes(4,i,1) = 12
      processes(5,i,1) = -13
      processes(6,i,1) = 13
      processes(7,i,1) = 21 ! gluon
      enddo
      !dir[6]=DbargluUbarWpZ-fin-Hahn
      do i = 6, 6 ! fin_ind(6): jet at the 3rd place
      processes(1,i,1) = -1
      processes(2,i,1) = 21
      processes(3,i,1) = -2
      processes(4,i,1) = -11
      processes(5,i,1) = 12
      processes(6,i,1) = -13
      processes(7,i,1) = 13
      enddo
      ! dir[7]=UgluDWpZ-fin-Hahn
      do i = 7, 7 ! fin_ind(7): jet at the 3rd place
      processes(1,i,1) = 2
      processes(2,i,1) = 21
      processes(3,i,1) = 1
      processes(4,i,1) = -11
      processes(5,i,1) = 12
      processes(6,i,1) = -13
      processes(7,i,1) = 13
      enddo
      ! EW corrections:
      do i = 8, 10 ! virt(8), coll_rad (9), coll_ind(10)
      processes(:,i,1) = processes(:,1,1)
      enddo
      ! dir[11]=UDbarWpZgam-fin-Hahn
      do i = 11, 11 ! fin_rad(11): photon at the end
      processes(1,i,1) = 2
      processes(2,i,1) = -1
      processes(3,i,1) = -11
      processes(4,i,1) = 12
      processes(5,i,1) = -13
      processes(6,i,1) = 13
      processes(7,i,1) = 22 ! photon
      enddo
      ! dir[12]=DbargamUbarWpZ-fin-Hahn
      do i = 12, 12 ! fin_ind(12): jet at the 3rd place
      processes(1,i,1) = -1
      processes(2,i,1) = 22 ! photon
      processes(3,i,1) = -2
      processes(4,i,1) = -11
      processes(5,i,1) = 12
      processes(6,i,1) = -13
      processes(7,i,1) = 13
      enddo
      ! dir[13]=UgamDWpZ-fin-Hahn
      do i = 13, 13 ! fin_ind(13): jet at the 3rd place
      processes(1,i,1) = 2
      processes(2,i,1) = 22 ! photon
      processes(3,i,1) = 1
      processes(4,i,1) = -11
      processes(5,i,1) = 12
      processes(6,i,1) = -13
      processes(7,i,1) = 13
      enddo

      ! dir[14]=UDbarWpZ-virtEW-decayZ-Ninh
      processes(:,14,1) = processes(:,8,1)
      ! dir[15]=UDbarWpZgam-fin-decayZ-Hahn
      processes(:,15,1) = processes(:,11,1)
      ! dir[16]=UDbarWpZ-virtEW-decayW-Ninh
      processes(:,16,1) = processes(:,8,1)
      ! dir[17]=UDbarWpZgam-fin-decayW-Hahn
      processes(:,17,1) = processes(:,11,1)
      n_cont(1) = 17
!!!]process 1: WpZ
!!![process 2: WmZ
      processes(:,:,2) = processes(:,:,1)
      n_cont(2) = n_cont(1)
      do i = 0, n_cont(2) ! over contributions
       do j = 1, 7 ! over legs
        ! for quarks
        IF ((i-6)*(i-7).le.0.or.(i-12)*(i-13).le.0) THEN
        ! induced n+1 processes (q glu > q glu) 6,7; 12, 13
        if (processes(j,i,2).eq.-1) then ! anti-down
         processes(j,i,2) = -2 ! anti-up
        else if (processes(j,i,2).eq.-2) then ! anti-up
         processes(j,i,2) = -1 ! anti-down
        else if (processes(j,i,2).eq.1) then ! down
         processes(j,i,2) = 2 ! up   
        else if (processes(j,i,2).eq.2) then ! up
         processes(j,i,2) = 1 ! down
        endif        
        ELSE 
        ! radiation processes: UDbar > UbarD
        if (abs(processes(j,i,2)).ge.1.and.abs(processes(j,i,2)).le.6) then ! quark
         processes(j,i,2) = -processes(j,i,2)
        endif
        ENDIF
        ! for leptons
        if (abs(processes(j,i,2)).ge.11.and.abs(processes(j,i,2)).le.12) then ! nu_e, electron
         processes(j,i,2) = -processes(j,i,2)
        endif
       enddo
      enddo
!!!]process 2: WmZ
!!![process 3: WW_uu
      processes(:,:,3) = processes(:,:,1)
      n_cont(3) = n_cont(1)
      do i = 0, n_cont(3) ! over contributions
       do j = 1, 7 ! over legs
        ! for quarks
        IF ((i-6)*(i-7).le.0.or.(i-12)*(i-13).le.0) THEN
        ! induced n+1 processes (q glu > q glu) 6,7; 12, 13
        ! d > u
        if (processes(j,i,3).eq.-1) then ! anti-down
         processes(j,i,3) = -2 ! anti-up
        else if (processes(j,i,3).eq.1) then ! down
         processes(j,i,3) = 2 ! up
        endif
        ELSE
        ! radiation processes: UDbar > UbarU
        if (processes(j,i,3).eq.2) then ! up
         processes(j,i,3) = -2 ! anti-up
        else if (processes(j,i,3).eq.-1) then ! anti-down
         processes(j,i,3) = 2 ! up
        endif
        ENDIF
        ! for leptons
        if (processes(j,i,3).eq.-13) then ! mu+
         processes(j,i,3) = 13 ! mu-
        else if (processes(j,i,3).eq.13) then ! mu-
         processes(j,i,3) = -14 ! anti nu_mu
        endif
       enddo
      enddo
!!!]process 3: WW_uu
!!![process 4: WW_dd
      processes(:,:,4) = processes(:,:,3)
      n_cont(4) = n_cont(3)
      do i = 0, n_cont(4) ! over contributions
       do j = 1, 7 ! over legs
       ! for quarks
        if (processes(j,i,4).eq.2) then ! up
         processes(j,i,4) = 1 ! down
        else if (processes(j,i,4).eq.-2) then ! anti-up
         processes(j,i,4) = -1 ! anti-down
        endif
       enddo
      enddo
!!!]process 4: WW_dd
!!![process 5: WW_gg
      processes(:,:,5) = processes(:,:,4)
      n_cont(5) = 1
      do i = 0, n_cont(5) ! over contributions
       do j = 1, 7 ! over legs
       ! for quarks
        if (processes(j,i,5).eq.1) then ! down
         processes(j,i,5) = 21 ! gluon
        else if (processes(j,i,5).eq.-1) then ! anti-down
         processes(j,i,5) = 21 ! gluon
        endif
       enddo
      enddo
!!!]process 5: WW_gg
!!![process 6: WW_bb
      processes(:,:,6) = processes(:,:,4)
      n_cont(6) = n_cont(4) + 4
      do i = 0, n_cont(4) ! over contributions
       do j = 1, 7 ! over legs
       ! for quarks
        if (processes(j,i,6).eq.1) then ! down
         processes(j,i,6) = 5 ! bottom
        else if (processes(j,i,6).eq.-1) then ! anti-down
         processes(j,i,6) = -5 ! anti-bottom
        endif
       enddo
      enddo
      ! for TW processes
      do i = 18, 19
       processes(:,i,6) = processes(:,i-12,6)
      enddo
      do i = 20, 21
       processes(:,i,6) = processes(:,i-8,6)
      enddo
!!!]process 6: WW_bb
!!![process 7: WW_aa
      processes(:,:,7) = processes(:,:,5)
      n_cont(7) = 1
      do i = 0, n_cont(7) ! over contributions
       do j = 1, 7 ! over legs
       ! for gluon
        if (processes(j,i,7).eq.21) then ! gluon
         processes(j,i,7) = 22 ! photon
        endif
       enddo
      enddo
!!!]process 7: WW_aa
!!![process 8: ZZ_uu
      processes(:,:,8) = processes(:,:,3)
      n_cont(8) = n_cont(3)
      do i = 0, n_cont(8) ! over contributions
       do j = 1, 7 ! over legs
       ! W+W- > ZZ: e+ > e+, nu_e > e-, mu- > mu+, anti-nu_mu > mu-
        if (processes(j,i,8).eq.12) then ! nu_e
         processes(j,i,8) = 11 ! e-
        else if (processes(j,i,8).eq.13) then ! muon-
         processes(j,i,8) = -13 ! muon +
        else if (processes(j,i,8).eq.-14) then ! anti nu_mu
         processes(j,i,8) = 13 ! mu-
        endif
       enddo
      enddo
!!!]process 8: ZZ_uu
!!![process 9: ZZ_dd
      processes(:,:,9) = processes(:,:,8)
      n_cont(9) = n_cont(8)
      do i = 0, n_cont(9) ! over contributions
       do j = 1, 7 ! over legs
       ! for quarks
        if (processes(j,i,9).eq.2) then ! up
         processes(j,i,9) = 1 ! down
        else if (processes(j,i,9).eq.-2) then ! anti-up
         processes(j,i,9) = -1 ! anti-down
        endif
       enddo
      enddo
!!!]process 9: ZZ_dd
!!![process 10: ZZ_gg
      processes(:,:,10) = processes(:,:,9)
      n_cont(10) = 1
      do i = 0, n_cont(10) ! over contributions
       do j = 1, 7 ! over legs
       ! for quarks
        if (processes(j,i,10).eq.1) then ! down
         processes(j,i,10) = 21 ! gluon
        else if (processes(j,i,10).eq.-1) then ! anti-down
         processes(j,i,10) = 21 ! gluon
        endif
       enddo
      enddo
!!!]process 10: ZZ_gg
!!![process 11: ZZ_bb
      processes(:,:,11) = processes(:,:,9)
      n_cont(11) = n_cont(9)
      do i = 0, n_cont(11) ! over contributions
       do j = 1, 7 ! over legs
       ! for quarks
        if (processes(j,i,11).eq.1) then ! down
         processes(j,i,11) = 5 ! bottom
        else if (processes(j,i,11).eq.-1) then ! anti-down
         processes(j,i,11) = -5 ! anti-bottom
        endif
       enddo
      enddo
!!!]process 11: ZZ_bb
!!![process 12: ZZ_aa
      processes(:,:,12) = processes(:,:,10)
      n_cont(12) = 0 ! zero contribution to the polarized cross sections
      do i = 0, n_cont(12) ! over contributions
       do j = 1, 7 ! over legs
       ! for gluon
        if (processes(j,i,12).eq.21) then ! gluon
         processes(j,i,12) = 22 ! photon
        endif
       enddo
      enddo
!!!]process 12: ZZ_aa
!!![process 13: WW_ggbb_4FS
      processes(:,:,13) = processes(:,:,5) ! load from gg > W+W-
      n_cont(13) = 1
      do i = 0, n_cont(13) ! over contributions
       ! additional legs
       processes(7,i,13) = -5 ! anti-bottom
       processes(8,i,13) = -processes(7,i,13) ! bottom
      enddo
!!!]process 13: WW_ggbb_4FS
!!![process 14: WW_uubb_4FS
      processes(:,:,14) = processes(:,:,13)
      n_cont(14) = 1
      do i = 0, n_cont(14) ! over contributions
       processes(1,i,14) = -2 ! anti-up
       processes(2,i,14) = 2
      enddo
!!!]process 14: WW_uubb_4FS
!!![process 15: WW_ddbb_4FS
      processes(:,:,15) = processes(:,:,14)
      n_cont(15) = 1
      do i = 0, n_cont(15)
       processes(1,i,15) = -1 ! anti-down
       processes(2,i,15) = 1
      enddo
!!!]process 15: WW_ddbb_4FS
!!![process 16: WW_aabb_4FS
      processes(:,:,16) = processes(:,:,15)
      n_cont(16) = 1
      do i = 0, n_cont(16)
       processes(1,i,16) = 22 ! photon
       processes(2,i,16) = 22
      enddo
!!!]process 16: WW_aabb_4FS
!!![process 51: WpWmWm
      i_proc = 51 ! process id
      ! born and NLO QCD corrections:
      ! UbarDWpWmWm-born-Hahn
      do i_cont = 0, 1 ! full-born (0), DPA-born(1)
      processes(1,i_cont,i_proc) = -2 ! anti-up
      processes(2,i_cont,i_proc) = 1 ! down
      processes(3,i_cont,i_proc) = -11 ! e+
      processes(4,i_cont,i_proc) = 12 ! nu_e
      processes(5,i_cont,i_proc) = 13 ! mu-
      processes(6,i_cont,i_proc) = -14 ! anti-nu_mu
      processes(7,i_cont,i_proc) = 15 ! tau-
      processes(8,i_cont,i_proc) = -16 ! anti-nu_tau
      enddo
      ! virt_QCD(2), coll_rad (3), coll_ind(4), ...
      call add_contributions_NLOQCD(processes,i_proc)
      ! virt_EW(8), ...
      call add_contributions_NLOEW(processes,i_proc)
!!!]process 51: WpWmWm
!!![process 52: WmWpWp
      i_proc = 52 ! process id
      ! born and NLO QCD corrections:
      ! DbarUWmWpWp-born-Hahn
      do i_cont = 0, 1 ! full-born (0), DPA-born(1)
      processes(1,i_cont,i_proc) = -1 ! anti-down
      processes(2,i_cont,i_proc) = 2 ! up
      processes(3,i_cont,i_proc) = 11 ! e-
      processes(4,i_cont,i_proc) = -12 ! anti-nu_e
      processes(5,i_cont,i_proc) = -13 ! mu+
      processes(6,i_cont,i_proc) = 14 ! nu_mu
      processes(7,i_cont,i_proc) = -15 ! tau+
      processes(8,i_cont,i_proc) = 16 ! nu_tau
      enddo
      ! virt_QCD(2), coll_rad (3), coll_ind(4), ...
      call add_contributions_NLOQCD(processes,i_proc)
      ! virt_EW(8), ...
      call add_contributions_NLOEW(processes,i_proc)
!!!]process 52: WmWpWp
      end subroutine load_data_processes

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine add_contributions_NLOQCD(processes,i_proc)
! add NLOQCD contributions: idc = 2,3,...,7
      implicit none
      integer i_proc, i_cont, i_leg, nlegs_born, j_leg
!!![List of processes
! use zero (0) for no particle
! follow MC particle numbering scheme (21 = gluon, 22 = photon, ...)      
      integer :: processes(nlegs_max,0:n_cont_max,n_id_max)
!!!]List of processes

! born-like kinematics:
      do i_cont = 2, 4 ! virt (2), rad-coll(3), ind-coll(4)
      processes(:,i_cont,i_proc) = processes(:,1,i_proc)
      enddo
! NLO gluon radiation (idc = 5):
      nlegs_born = 0
      do i_leg = 1, nlegs_max
       if (processes(i_leg,1,i_proc).ne.0) then
        nlegs_born = nlegs_born + 1
       endif
      enddo
      processes(:,5,i_proc) = processes(:,1,i_proc)
      processes(nlegs_born+1,5,i_proc) = 21 ! gluon

! NLO gluon induced process (Q_bar glu, idc = 6):
      processes(1,6,i_proc) = processes(1,1,i_proc)
      processes(3,6,i_proc) = -processes(2,1,i_proc) ! crossing quark to 3rd position
      processes(2,6,i_proc) = 21 ! gluon at 2nd position
      do i_leg = 4, nlegs_born+1
       j_leg = i_leg - 1
       processes(i_leg,6,i_proc) = processes(j_leg,1,i_proc)
      enddo

! NLO gluon induced process (Q glu, idc = 7):
      processes(1,7,i_proc) = processes(2,1,i_proc)
      processes(3,7,i_proc) = -processes(1,1,i_proc) ! crossing anti-quark to 3rd position
      processes(2,7,i_proc) = 21 ! gluon at the 2nd position
      do i_leg = 4, nlegs_born+1
       j_leg = i_leg - 1
       processes(i_leg,7,i_proc) = processes(j_leg,1,i_proc)
      enddo
      
      end subroutine add_contributions_NLOQCD
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine add_contributions_NLOEW(processes,i_proc)
! add NLOEW contributions: idc = 8,9,...
      implicit none
      integer i_proc, i_cont, i_leg, nlegs_born, j_leg, i_res, j_cont
!!![List of processes
! use zero (0) for no particle
! follow MC particle numbering scheme (21 = gluon, 22 = photon, ...)      
      integer :: processes(nlegs_max,0:n_cont_max,n_id_max)
!!!]List of processes

! born-like kinematics:
      do i_cont = 8, 10 ! virt-EW (8), rad-coll(9), ind-coll(10)
      processes(:,i_cont,i_proc) = processes(:,1,i_proc)
      enddo
! NLO photon radiation (idc = 11):
      nlegs_born = 0
      do i_leg = 1, nlegs_max
       if (processes(i_leg,1,i_proc).ne.0) then
        nlegs_born = nlegs_born + 1
       endif
      enddo
      processes(:,11,i_proc) = processes(:,1,i_proc)
      processes(nlegs_born+1,11,i_proc) = 22 ! photon

! NLO photon induced process (Q_bar photon, idc = 12):
      processes(1,12,i_proc) = processes(1,1,i_proc)
      processes(3,12,i_proc) = -processes(2,1,i_proc) ! crossing quark to 3rd position
      processes(2,12,i_proc) = 22 ! photon at 2nd position
      do i_leg = 4, nlegs_born+1
       j_leg = i_leg - 1
       processes(i_leg,12,i_proc) = processes(j_leg,1,i_proc)
      enddo

! NLO photon induced process (Q glu, idc = 7):
      i_cont = 13
      processes(1,i_cont,i_proc) = processes(2,1,i_proc)
      processes(3,i_cont,i_proc) = -processes(1,1,i_proc) ! crossing anti-quark to 3rd position
      processes(2,i_cont,i_proc) = 22 ! photon at the 2nd position
      do i_leg = 4, nlegs_born+1
       j_leg = i_leg - 1
       processes(i_leg,i_cont,i_proc) = processes(j_leg,1,i_proc)
      enddo

! NLO EW decays:
      do i_res = 1, nlegs_res 
      i_cont = 13 + 2*i_res - 1 ! virtual corrections
      processes(:,i_cont,i_proc) = processes(:,1,i_proc)
      j_cont = i_cont + 1 ! real corrections
      processes(:,j_cont,i_proc) = processes(:,5,i_proc)
      enddo
      
      end subroutine add_contributions_NLOEW
!**********************************************************************
end module
