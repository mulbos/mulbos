	subroutine OSmap_VV_2to5_decayV1(fac)
        use global_vars
        use global_mom_xpdf

	implicit none

	double precision fac
	integer i, j
c For 1 final jets and 2 gauge bosons (4 fermions)
        double precision rext(0:3,proc_nlegs_os), Vext(0:3,proc_nlegs)
        double precision rm2(2)

! proc_nlegs = 7 
! proc_nlegs_os = 9

! Note: pext(6,4) are the true qq -> 4l + j momenta
       do i = 1,7 
         do j = 1,4
          Vext(j-1,i) = pext(i,j)
         enddo
        enddo

!!![finalizing for OS projection
       do i = 0, 3
       rext(i,1) = Vext(i,1)
       rext(i,2) = Vext(i,2)
       enddo
       
        rm2(1) = mass_VV(1)**2 ! mass squared input
        rm2(2) = mass_VV(2)**2 ! mass squared input

! use Denner ZZ 2021 paper:
      call OSMap_decayV1_NLO(Vext(0,3),rm2(1),rext(0,3))

!!![re-ordering
! store oext(proc_nlegs_os,4) in the common block in xsection.h
        do j = 1,4
         do i = 1, 4 ! first 4 particles qq > V1 V2
          oext(i,j) = rext(j-1,i)
         enddo
          oext(5,j) = rext(j-1,9) ! the photon
         do i = 6, 9 ! the 4 leptons
          oext(i,j) = rext(j-1,i-1)
         enddo
        enddo
!!!]re-ordering

! store os momenta in the VV CMS
      oext_vv = oext

! store os momenta in the Lab frame
      oext_pp = oext
      ! boost from qq partonic cms to lab frame
      do i = 1, proc_nlegs_os
       call boost0(rgam_boost,-rbetgam_boost,oext_pp(i,1),oext_pp(i,2),oext_pp(i,3),oext_pp(i,4))
      enddo
!!!]finalizing for OS projection
        return
	end subroutine
***********************************************************************
      subroutine OSMap_decayV1_NLO(Vext,rm2,rext)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C input: Vext(0:3,3:7) are the momenta of 4 leptons + 1 photon in 
C the rest frame of X = V1 + V2; rm2(2) are the physical 
C masses of the 2 gauge bosons. 
C output: rext(0:3,7) are the OS-projected momenta of the 
C two gauge bosons (first 2 entries) and of the 4 leptons + 1 photon.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      double precision rext(0:3,7), Vext(0:3,3:7), rm2(2)
      double precision tmp1, tmp2, Vsum(0:3), SX, sqrtSX
      integer i
      double precision rm(2), p(0:3,3), Lep(0:3,3)

      do i = 1, 2
      rm(i) = sqrt(rm2(i))
      enddo

      do i = 0, 3
       Vsum(i) = Vext(i,3)+Vext(i,4)+Vext(i,5)+Vext(i,6)+Vext(i,7)
      enddo
      SX = Vsum(0)**2 - Vsum(1)**2 - Vsum(2)**2 - Vsum(3)**2
      sqrtSX = sqrt(SX)

c      print*,"sqrtSX child=", sqrtSX

***[first on-shell gauge boson
       tmp1 = (rm2(1) - rm2(2))/(2*sqrtSX) + sqrtSX/2 ! r10
       rext(0,1) = tmp1

       tmp2 = (Vext(1,3) + Vext(1,4) + Vext(1,7))**2
     & + (Vext(2,3) + Vext(2,4) + Vext(2,7))**2 
     & + (Vext(3,3) + Vext(3,4) + Vext(3,7))**2

       tmp2 = (rext(0,1)**2 - rm2(1))/tmp2
       tmp2 = sqrt(tmp2) ! b

       do i = 1, 3
       rext(i,1) = tmp2*(Vext(i,3) + Vext(i,4) + Vext(i,7))
       enddo
***]first on-shell gauge boson

***[second on-shell gauge boson
       do i = 0, 3
       rext(i,2) = Vsum(i) - rext(i,1)
       enddo
***]second on-shell gauge boson

***[projected lepton momenta (V1 > r3 + r4 + r7)
       p(:,1) = Vext(:,3)
       p(:,2) = Vext(:,4)
       p(:,3) = Vext(:,7)
       call OSMap_decayV_NLO(p(0,1),rm(1),rext(0,1),Lep(0,1))
       do i = 0, 3
        rext(i,3) = Lep(i,1)
        rext(i,4) = Lep(i,2)
        rext(i,7) = Lep(i,3)
       enddo
***]projected lepton momenta (V1 > r3 + r4 + r7)

***[projected lepton momenta (V2 > r5 + r6)
       call OSMap_decayV_LO(Vext(0,5),rm(2),rext(0,2),rext(0,5))
***]projected lepton momenta (V2 > r5 + r6)

      END subroutine OSMap_decayV1_NLO
