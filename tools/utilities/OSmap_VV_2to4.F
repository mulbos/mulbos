      subroutine OSmap_VV_2to4(fac)
      use global_vars
      use global_mom_xpdf
c General OS mapping for pp > V1 V2 > 4 leptons
c Input: 
c	fac: some initial weight
c Output:
c	fac: updated weight
c	all external momenta calculated, see common in xsection.h
      implicit none
      integer i, j
      double precision fac
c For 2 gauge bosons and 4 leptons:
      double precision rext(0:3,proc_nlegs_os), Vext(0:3,proc_nlegs)

      double precision RM2(2)
      
      rm2(1) = mass_VV(1)**2        ! mass squared input
      rm2(2) = mass_VV(2)**2        ! mass squared input
      
C Note: pext(6,4) are the true qq -> 4l momenta in the partonic CMS
       do i = 1,proc_nlegs
         do j = 1,4
           Vext(j-1,i) = pext(i,j)
         enddo
        enddo
        
       do i = 0, 3
       rext(i,1) = Vext(i,1)
       rext(i,2) = Vext(i,2)
       enddo

      call OSProject(osmap_case,Vext(0,3),rm2(1),rext(0,3))

! momenta in the partonic CMS
        do i = 1, proc_nlegs_os
C store oext(LEGSX,4) in the common block in xsection.h
         do j = 1,4
          oext(i,j) = rext(j-1,i)
         enddo
        enddo
! store os momenta in the VV CMS
      oext_vv = oext
! store os momenta in the Lab frame
      oext_pp = oext
      ! boost from qq partonic cms to lab frame
      do i = 1, proc_nlegs_os
       call boost0(rgam_boost,-rbetgam_boost,oext_pp(i,1),oext_pp(i,2),oext_pp(i,3),oext_pp(i,4))
      enddo
        
      return
      end subroutine
