module global_mom_xpdf
      use global_vars
      double precision tau_pdf, x_pdf, y_coll
      double precision x1_pdf, x2_pdf
      double precision rgam_boost, rbetgam_boost ! LAB to partonic CMS
      double precision pext(nlegs_max,4) ! exact momenta in the partonic CMS (for cut, hist), same particle order as processes() in global_vars
      double precision vext_jend(5,nlegs_max) ! same as pext() but put the jet/photon at the end, exchange row-column, add particle ID
      double precision oext(nlegs_os_max,4)! on-shell momenta in the partonic CMS
      double precision oext_vv(nlegs_os_max,4)! on-shell momenta in the VV CMS
      double precision oext_pp(nlegs_os_max,4)! on-shell momenta in the proton-proton CMS (Lab frame)
      !!!
      double precision pext_tilde(nlegs_max,4) ! exact momenta in the parent partonic CMS (for cut, hist), same particle order as LO
      double precision oext_tilde(nlegs_os_max,4)! on-shell momenta in the parent partonic CMS (different from the VV cms)
      double precision oext_tilde_vv(nlegs_os_max,4)! on-shell momenta in the VV CMS (= child partonic CMS)
      double precision oext_tilde_pp(nlegs_os_max,4)! on-shell momenta in the proton-proton CMS (Lab frame)
end module
