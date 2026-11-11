* basic quantities for external particles
* last modified 27 Oct 2025

      integer N_legs_max
      parameter (N_legs_max = 20)
      double complex vec_K(2,2,N_legs_max)
      double complex vec_E(2,2,-1:1,N_legs_max)
      double complex vec_EC(2,2,-1:1,N_legs_max)
      double complex spin(2,-1:1,6:7,N_legs_max)
      double complex dspin(2,-1:1,6:7,N_legs_max)
      common /K_Eps_Spin/ vec_K, vec_E, vec_EC, spin, dspin

