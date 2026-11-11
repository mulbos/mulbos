module global_cuts
double precision minv_4l_min, minv_4l_max
double precision minv_4lep_min(3), minv_4lep_max(3)
double precision DRLA, py_phot_recomb, p_recomb
double precision py_jet_veto_max, pt_jet_veto_min
double precision minv_LL12_min, minv_LL12_max
double precision minv_LL34_min, minv_LL34_max
double precision pTmin_L(6), pTmin_Llead, pTmin_Llead2
double precision pymax_L(6)
double precision m_V1_shift, m_V2_shift
double precision mT_V1_min, mT_V2_min
double precision DRmin_LL_all
double precision DRmin_LL_12,DRmin_LL_13,DRmin_LL_14,DRmin_LL_23,DRmin_LL_24,DRmin_LL_34
double precision pT_4L_min,pT_4L_max
double precision pT_LL12_min,pT_LL12_max,pT_LL34_min,pT_LL34_max
double precision pT_LL13_min,pT_LL13_max
double precision pT_miss_min, pT_miss_max
double precision minv_LL13_min, minv_LL13_max

contains

      subroutine Read_Cuts(path_file)
      use global_vars
      use readinput
      implicit none
      character*(*) path_file

        call loadfile(path_file)
        ! photon threshold for lepton-photon recombination
        call read_Real("RECOMB_PY_A_MAX",py_phot_recomb,5d0)
        ! lepton-photon recombination parameter        
        call read_Real("RECOMB_DR_LA",DRLA,0.1d0)
        ! p = 0d0 using Delta_R distance
        ! p = 1d0 the kt distance (algorithm)
        ! p = -1d0 the anti-kt algorithm
        call read_Real("RECOMB_EXP",p_recomb,0d0)
        
        ! default: no jet veto
        ! pseudo-rapidity max for jet veto
        call read_Real("VETO_PY_J_MAX",py_jet_veto_max,0d0)
        ! pt min for jet veto
        call read_real("VETO_PT_J_MIN",pt_jet_veto_min,13d10)
        
        ! this is a technical cut to probe the MW resonance
        call read_Real("MINV_4L_MIN",minv_4l_min,0d0)
        call read_Real("MINV_4L_MAX",minv_4l_max,13D10)
        call read_Real("MINV_4L_1234_MIN",minv_4lep_min(3),0d0)
        call read_Real("MINV_4L_1234_MAX",minv_4lep_max(3),13D10)

        if (lTriBoson) then
        call read_Real("MINV_4L_1256_MIN",minv_4lep_min(2),0d0)
        call read_Real("MINV_4L_3456_MIN",minv_4lep_min(1),0d0)
        call read_Real("MINV_4L_1256_MAX",minv_4lep_max(2),13D10)
        call read_Real("MINV_4L_3456_MAX",minv_4lep_max(1),13D10)
        endif
        
        call read_Real("MINV_LL12_MIN",minv_LL12_min,0d0)
        call read_Real("MINV_LL12_MAX",minv_LL12_max,13D10)
        
        call read_Real("MINV_LL34_MIN",minv_LL34_min,0d0)
        call read_Real("MINV_LL34_MAX",minv_LL34_max,13D10)
        
        call read_Real("MINV_LL13_MIN",minv_LL13_min,0d0)
        call read_Real("MINV_LL13_MAX",minv_LL13_max,13D10)
        
        ! a cut around MV pole
        call read_Real("MINV_LL12_SHIFT",m_V1_shift,10D0)
        call read_Real("MINV_LL34_SHIFT",m_V2_shift,10D0)
        ! a cut on MT_V
        call read_Real("MT_LL12_MIN",mT_V1_min,0D0)
        call read_Real("MT_LL34_MIN",mT_V2_min,0D0)
        
        call read_Real("PT_L1_MIN",pTmin_L(1),20d0)
        call read_Real("PT_L2_MIN",pTmin_L(2),20d0)
        call read_Real("PT_L3_MIN",pTmin_L(3),20d0)
        call read_Real("PT_L4_MIN",pTmin_L(4),20d0)
        call read_Real("PT_LLEAD_MIN",pTmin_Llead,20d0)
        call read_Real("PT_LLEAD2_MIN",pTmin_Llead2,20d0)
        if (lTriBoson) then
        call read_Real("PT_L5_MIN",pTmin_L(5),20d0)
        call read_Real("PT_L6_MIN",pTmin_L(6),20d0)
        endif
        
        call read_Real("PY_L1_MAX",pymax_L(1),2.5d0)
        call read_Real("PY_L2_MAX",pymax_L(2),2.5d0)
        call read_Real("PY_L3_MAX",pymax_L(3),2.5d0)
        call read_Real("PY_L4_MAX",pymax_L(4),2.5d0)
        if (lTriBoson) then
        call read_Real("PY_L5_MAX",pymax_L(5),2.5d0)
        call read_Real("PY_L6_MAX",pymax_L(6),2.5d0)
        endif

        if (lTriBoson) then
        call read_Real("DR_LL_ALL_MIN",DRmin_LL_all,0.2d0) ! for all charged leptons
        endif
        call read_Real("DR_LL12_MIN",DRmin_LL_12,0.2d0)
        call read_Real("DR_LL13_MIN",DRmin_LL_13,0.2d0)
        call read_Real("DR_LL14_MIN",DRmin_LL_14,0.2d0)
        call read_Real("DR_LL23_MIN",DRmin_LL_23,0.2d0)
        call read_Real("DR_LL24_MIN",DRmin_LL_24,0.2d0)
        call read_Real("DR_LL34_MIN",DRmin_LL_34,0.2d0)
        
        ! PT of 4L system
        call read_Real("PT_4L_MIN",pT_4L_min,0d0)
        call read_Real("PT_4L_MAX",pT_4L_max,13d10)
        
        ! PT of intermediate gauge bosons
        call read_Real("PT_LL12_MIN",pT_LL12_min,0d0)
        call read_Real("PT_LL12_MAX",pT_LL12_max,13d10)
        call read_Real("PT_LL34_MIN",pT_LL34_min,0d0)
        call read_Real("PT_LL34_MAX",pT_LL34_max,13d10)
        call read_Real("PT_LL13_MIN",pT_LL13_min,0d0)
        call read_Real("PT_LL13_MAX",pT_LL13_max,13d10)
        
        ! PT miss
        call read_Real("PT_MISS_MIN",pT_miss_min,0D0)
        call read_Real("PT_MISS_MAX",pT_miss_max,13d10)

        call closefile
      end subroutine

end module
