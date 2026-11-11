      subroutine OSmap_VVV_2to6(fac)
      use global_vars
      use global_mom_xpdf
c An OS mapping for pp > V1 V2 V3 > 6 leptons
c Input: 
c	fac: some initial weight
c Output:
c	fac: updated weight
c	all external momenta calculated, see global_mom_xpdf
      implicit none
      integer i, j
      double precision fac
c For all resonance gauge bosons:
      double precision rext(0:3,proc_nlegs_os), Vext(0:3,proc_nlegs)
      double precision dot_p2

C Note: pext(6,4) are the true qq -> 4l momenta in the partonic CMS
       do i = 1,proc_nlegs
         do j = 1,4
           Vext(j-1,i) = pext(i,j)
         enddo
        enddo

! initial-state particles:
       do i = 0, 3
       rext(i,1) = Vext(i,1)
       rext(i,2) = Vext(i,2)
       enddo

! resonance and final-state particles:
      call OSProject_VVV(Vext(0,3),rext(0,3),fac)
      if (fac .eq. 0d0) then
       return
      endif

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
      end subroutine OSmap_VVV_2to6
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine OSProject_VVV(Vext,rext,fac)
      use global_vars
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C input: Vext(0:3,3:8) are the momenta of 4 leptons in 
C the rest frame of X = V1 + V2 + V3. 
C output: rext(0:3,9) are the OS-projected momenta of the 
C three gauge bosons (first 3 entries) and of the 6 leptons, 
C in the rest frame of X.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      double precision fac
      double precision Vext(0:3,3:8) ! input lepton off-shell momenta
      double precision rext(0:3,9) ! output on-shell momenta (resonances + leptons) 
      double precision Vsum(0:3), SX, mVVV
      double precision m12,m13,m23
      double precision Vint(0:3,3) ! resonances
      double precision inv_mass
      double precision dot_p2

      Vint(:,1) = Vext(:,3)+Vext(:,4) 
      Vint(:,2) = Vext(:,5)+Vext(:,6)
      Vint(:,3) = Vext(:,7)+Vext(:,8)
      Vsum(:) = Vint(:,1)+Vint(:,2)+Vint(:,3)
      SX = dot_p2(Vsum(0))
      mVVV = sqrt(SX)
      m12 = inv_mass(Vint(0,1),Vint(0,2))
      m13 = inv_mass(Vint(0,1),Vint(0,3))
      m23 = inv_mass(Vint(0,2),Vint(0,3)) 

      call OSMap_Res_123(Vint(0,1),mVVV,m12,m13,m23,rext(0,1),fac)
      if (fac .eq. 0d0) then
       return
      endif

***[projected lepton momenta as Denner ZZ 2021
      call OSMap_decayV_LO(Vext(0,3),mass_VV(1),rext(0,1),rext(0,4))
      call OSMap_decayV_LO(Vext(0,5),mass_VV(2),rext(0,2),rext(0,6))
      call OSMap_decayV_LO(Vext(0,7),mass_VV(3),rext(0,3),rext(0,8))
***]projected lepton momenta as Denner ZZ 2021

!      do i = 1, 9
!       print*,"projected OS moms =",rext(:,i),sqrt(abs(dot_p2(rext(0,i))))
!      enddo
!      stop

      end subroutine OSProject_VVV
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine OSMap_Res_123(Vint,mVVV,m12,m13,m23,p_os,fac)
C input: Vint(0:3,3) are the momenta of 3 resonances in 
C the rest frame of X = V1 + V2 + V3. 
C mVVV: invariant mass of the 3 resonances
C output: p_os(0:3,3)
      use global_vars
      implicit none
      double precision fac
      double precision mVVV,m12,m13,m23
      double precision p_os(0:3,3)
      integer iA, iB1, iB2, iB, iC, i
      double precision mA, mB
      double precision pAB(0:3,2),rm2(2)
      double precision pBi(0:3,2) ! in the CMS_VVV
      double precision pABC(0:3,3), pA(0:3), pB(0:3)
      double precision MV(3), inv_mass_p
      double precision Vint(0:3,3) ! resonances
      double precision alpha_k ! re-scaling factor of the spatial vectors
      double precision alpha2_k, m_sq(3), kk_sq(3)
      integer :: n_res = 3

      IF (osmap_case.eq.1) THEN
!!![simultaneous on-shell projection
!!![find the re-scaling factor alpha_k
      do i = 1, n_res
       m_sq(i) = mass_VV(i)**2
       kk_sq(i) = Vint(1,i)**2+Vint(2,i)**2+Vint(3,i)**2
      enddo
!      print*,"n_res, osmap_case, mVVV =", n_res, osmap_case, mVVV
      call OSmap_solve_alpha2(n_res,mVVV,m_sq(1),kk_sq(1),alpha2_k)
!      print*,"alpha2_k =", alpha2_k
      alpha_k = sqrt(alpha2_k)
!!!]find the re-scaling factor alpha_k
!!![calculate the momenta in the on-shell VVV rest frame
! the on-shell VVV rest frame = the off-shell VVV rest frame
      do i = 1, n_res
       p_os(1:3,i) = Vint(1:3,i)*alpha_k
       p_os(0,i) = sqrt(m_sq(i) + kk_sq(i)*alpha2_k)
      enddo
!!!]calculate the momenta in the on-shell VVV rest frame
!      do i = 1, n_res
!      print*,"p_os(i) =",i,p_os(:,i),sqrt(p_os(0,i)**2-p_os(1,i)**2-p_os(2,i)**2-p_os(3,i)**2)
!      print*,"Vint(i) =",i,Vint(:,i),sqrt(Vint(0,i)**2-Vint(1,i)**2-Vint(2,i)**2-Vint(3,i)**2)
!      enddo
!      print*,"sum p_os =",p_os(:,1)+p_os(:,2)+p_os(:,3)
!      print*,"sum Vint =",Vint(:,1)+Vint(:,2)+Vint(:,3)
!      stop
      return
!!!]simultaneous on-shell projection
      ENDIF

      IF ((osmap_case-51)*(osmap_case-69).le.0) THEN
!!![choose A fist
! find the best pair for on-shell mapping:
      call Find_iA(mVVV,m12,m13,m23,iA,iB1,iB2,mA,mB,fac)
      if (fac .eq. 0d0) then
       return
      endif

***[on-shell Map AB
      pAB(:,1) = Vint(:,iA)
      pAB(:,2) = Vint(:,iB1) + Vint(:,iB2)
      rm2(1) = mA**2
      rm2(2) = mB**2       
      call OSMap_AB(pAB(0,1),rm2(1),pABC(0,1))
      pB(:) = pABC(:,2)
      pBi(:,1) = Vint(:,iB1)
      pBi(:,2) = Vint(:,iB2)
***]on-shell Map AB

***[decay the B into 2 massive on-shell particles
      MV(1) = mB
      MV(2) = mass_VV(iB1)
      MV(3) = mass_VV(iB2)
      call OSMap_decayV_AB(pBi(0,1),MV(1),pB(0),pABC(0,2))
***]decay the B into 2 massive on-shell particles

!!![re-order resonances
! from ABC to 123 (V1 to position 1, V2 to position 2, ...):
      call ReOrderRes(pABC(0,1),p_os(0,1),iA,iB1,iB2)
!!!]re-order resonances

      return
!!!]choose A fist
      ENDIF

      IF ((osmap_case-71)*(osmap_case-99).le.0) THEN
!!![choose AB fist
      call Find_AB(mVVV,m12,m13,m23,iA,iB,iC,fac)
      if (fac .eq. 0d0) then
       return
      endif

!      print*,"iA,iB,iC =",iA,iB,iC

      pABC(:,3) = Vint(:,iC)

!!![mapOS for AB
      pAB(:,1) = Vint(:,iA)
      pAB(:,2) = Vint(:,iB)
      rm2(1) = mass_VV(iA)**2
      rm2(2) = mass_VV(iB)**2
      call OSMap_AB_frameX(pAB(0,1),rm2(1),pABC(0,1))
!!!]mapOS for AB

!      do i = 1,3
!      print*,"check os pABC 1=",i,pABC(:,i),inv_mass_p(pABC(0,i))
!      enddo

!!![mapOS for C
      call OSMap_C_givenAB(pABC(0,1),mVVV,iA,iB,iC,fac)
      if (fac .eq. 0d0) then
       return
      endif
!!!]mapOS for C

!!![re-order resonances
! from ABC to 123 (V1 to position 1, V2 to position 2, ...):
      call ReOrderRes(pABC(0,1),p_os(0,1),iA,iB,iC)
!!!]re-order resonances

!!!]choose AB fist
      ENDIF

      end subroutine OSMap_Res_123
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine OSMap_C_givenAB(pABC,mVVV,iA,iB,iC,fac)
! input: pABC(1) is the mom of iA resonance, 
! pABC(2) is the mom of iB resonance, ...
      use global_vars
      implicit none
      double precision fac, pABC(0:3,3), mVVV, mXC, h_test, miC, miC_min
      double precision pXC(0:3,2), pXC_os(0:3,2), rm2(2), inv_mass
      integer i, imap(2), iA,iB,iC, iX

      iX = -1 ! initial value
      
      imap(1) = iA
      imap(2) = iB

      pXC(:,2) = pABC(:,3)
      rm2(2) = mass_VV(iC)**2

!!![mapping 71
! choose the pair with minimal invariant mass
      IF (osmap_case .eq. 71) THEN
      mXC = mVVV ! initial value
      
      do i = 1,2
       miC = inv_mass(pABC(0,i),pABC(0,3))
       miC_min = mass_VV(imap(i)) + mass_VV(iC)
       h_test = miC - miC_min

       if (h_test .gt. 0d0) then
        if (miC .lt. mXC) then
         mXC = miC ! update mXC
         iX = i
        endif
       endif
      enddo
      ENDIF
!!!]mapping 71

!!![mapping 72
! choose the pair with maximal invariant mass
      IF (osmap_case .eq. 72) THEN
      mXC = 0D0 ! initial value
      
      do i = 1,2
       miC = inv_mass(pABC(0,i),pXC(0,2))
       miC_min = mass_VV(imap(i)) + mass_VV(iC)
       h_test = miC - miC_min
       
       if (h_test .gt. 0d0) then
        if (miC .gt. mXC) then
         mXC = miC ! update mXC
         iX = i
        endif
       endif
      enddo
      ENDIF
!!!]mapping 72

      if (iX .lt. 0) then
       fac = 0
       return
      endif

! found a good iX, continue to check out:
!      print*,"iX =",iX, mXC

      pXC(:,1) = pABC(:,iX)
      rm2(1) = mass_VV(imap(iX))**2
      call OSMap_AB_frameX(pXC(0,1),rm2(1),pXC_os(0,1))
      pABC(:,3) = pXC_os(:,2)
      pABC(:,iX) = pXC_os(:,1)

      
      end subroutine OSMap_C_givenAB
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Find_AB(mVVV,m12,m13,m23,iA,iB,iC,fac)
      use global_vars
      implicit none
      double precision fac
      double precision mVVV,m12,m13,m23,m(3),mAB
      integer iA, iB, iC
      integer i, j_pair(2,3)
      double precision h_test, mjk, mjk_min

c      print*,"hello from Find_AB mVVV =",mVVV,m12,m13,m23
c      stop

      j_pair(1,1) = 2
      j_pair(2,1) = 3
      j_pair(1,2) = 1
      j_pair(2,2) = 3
      j_pair(1,3) = 1
      j_pair(2,3) = 2

      m(1) = m23
      m(2) = m13
      m(3) = m12

      iC = -1 ! initial value

!!![mapping 71
! choose the pair with minimal invariant mass
      IF (osmap_case .eq. 71) THEN
      mAB = mVVV ! initial value
      
      do i = 1,3
       mjk = m(i)
       mjk_min = mass_VV(j_pair(1,i)) + mass_VV(j_pair(2,i))
       h_test = mjk - mjk_min

       if (h_test .gt. 0d0) then
        if (mjk .lt. mAB) then
         mAB = mjk ! update mAB
         iC = i
        endif
       endif
      enddo
      ENDIF
!!!]mapping 71

!!![mapping 72
! choose the pair with maximal invariant mass
      IF (osmap_case .eq. 72) THEN
      mAB = 0D0 ! initial value
      
      do i = 1,3
       mjk = m(i)
       mjk_min = mass_VV(j_pair(1,i)) + mass_VV(j_pair(2,i))
       h_test = mjk - mjk_min
       
       if (h_test .gt. 0d0) then
        if (mjk .gt. mAB) then
         mAB = mjk ! update mAB
         iC = i
        endif
       endif
      enddo
      ENDIF
!!!]mapping 72

      if (iC .lt. 0) then
       fac = 0
       return
      endif

      iA = j_pair(1,iC)
      iB = j_pair(2,iC)


      end subroutine Find_AB
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine ReOrderRes(pABC,p123,iA,iB,iC)
      implicit none
      double precision pABC(0:3,3), p123(0:3,3)
      integer iA, iB, iC
      if (iA .eq. 1) then
       p123(:,1) = pABC(:,1)
       if (iB .eq. 2) then
        p123(:,2) = pABC(:,2)
        p123(:,3) = pABC(:,3)
       else
        p123(:,3) = pABC(:,2)
        p123(:,2) = pABC(:,3)
       endif
      else if (iA .eq. 2) then
       p123(:,2) = pABC(:,1)
       if (iB .eq. 1) then
        p123(:,1) = pABC(:,2)
        p123(:,3) = pABC(:,3)
       else
        p123(:,3) = pABC(:,2)
        p123(:,1) = pABC(:,3)
       endif
      else ! iA = 3
       p123(:,3) = pABC(:,1)
       if (iB .eq. 1) then
        p123(:,1) = pABC(:,2)
        p123(:,2) = pABC(:,3)
       else
        p123(:,2) = pABC(:,2)
        p123(:,1) = pABC(:,3)
       endif
      endif
      
      end subroutine ReOrderRes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Find_iA(mVVV,m12,m13,m23,iA,iB1,iB2,mA,mB,fac)
      use global_vars
      implicit none
      double precision fac
      double precision mVVV,m12,m13,m23,m(3)
      double precision h, h_test, mjk, mjk_min, mjk_max
      integer i, j_pair(2,3), iA, iB1, iB2
      double precision mA, mB

      j_pair(1,1) = 2
      j_pair(2,1) = 3
      j_pair(1,2) = 1
      j_pair(2,2) = 3
      j_pair(1,3) = 1
      j_pair(2,3) = 2

      m(1) = m23
      m(2) = m13
      m(3) = m12

      iA = -1

!!![mapping 51
! choose the pair with minimal invariant mass
      IF (osmap_case .eq. 51) THEN
      mB = mVVV ! initial value
      
      do i = 1,3
       mjk = m(i)
       mjk_min = mass_VV(j_pair(1,i)) + mass_VV(j_pair(2,i))
       mjk_max = mVVV - mass_VV(i)
       h_test = (mjk - mjk_min)*(mjk_max - mjk)

       if (h_test .gt. 0d0) then
        if (mjk .lt. mB) then
         mB = mjk ! update mB
         iA = i
        endif
       endif
      enddo
      ENDIF

!!!]mapping 51

!!![mapping 52
! choose the pair with maximal invariant mass
      IF (osmap_case .eq. 52) THEN
      mB = 0d0 ! initial value
      
      do i = 1,3
       mjk = m(i)
       mjk_min = mass_VV(j_pair(1,i)) + mass_VV(j_pair(2,i))
       mjk_max = mVVV - mass_VV(i)
       h_test = (mjk - mjk_min)*(mjk_max - mjk)
       if (h_test .gt. 0d0) then
        if (mjk .gt. mB) then
         mB = mjk ! update mB
         iA = i
        endif
       endif
      enddo
      ENDIF

!!!]mapping 52

!!![mapping 53
! choose the pair with maximal h
      IF (osmap_case .eq. 53) THEN
      h = -1d0

      do i = 1, 3
       mjk = m(i)
       mjk_min = mass_VV(j_pair(1,i)) + mass_VV(j_pair(2,i))
       mjk_max = mVVV - mass_VV(i)
       h_test = (mjk - mjk_min)*(mjk_max - mjk)
       if (h_test .gt. 0d0) then
        if (h_test .gt. h) then
         h = h_test ! update the h value
         iA = i
        endif
       endif
      enddo
      ENDIF
!!!]mapping 53

!!![mapping 69 (pair 23)
      IF (osmap_case .eq. 69) THEN
      i = 1
       mjk = m(i)
       mjk_min = mass_VV(j_pair(1,i)) + mass_VV(j_pair(2,i))
       mjk_max = mVVV - mass_VV(i)
       h_test = (mjk - mjk_min)*(mjk_max - mjk)
       if (h_test .gt. 0d0) then
        iA = i
       endif
      ENDIF
!!!]mapping 69 (pair 23)

!      print*,"end find_iA iA =",iA
!      stop

      if (iA .lt. 0) then
       fac = 0
       return
      endif

! output: iA,iB1,iB2,mA,mB

      iB1 = j_pair(1,iA)
      iB2 = j_pair(2,iA)
      
      mA = mass_VV(iA)
      mB = m(iA) 

!      print*,"final values Find_iA:",mB,iA
!      stop

      end subroutine Find_iA
      




