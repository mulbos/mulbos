	subroutine OSmap_TW_2to5(fac)
        use global_vars
        use global_mom_xpdf

	implicit none

	double precision fac
	integer i, j
c For 1 final jets and 2 gauge bosons (4 fermions)
        double precision rext(0:3,proc_nlegs_os), Vext(0:3,proc_nlegs)
        double precision sext(0:3,proc_nlegs_os)
        double precision qX(0:3), qX_inv(0:3)
        double precision rm2(2)
        logical lW1fromT

! proc_nlegs = 7 
! proc_nlegs_os = 9

! Note: pext(7,4) are the qq -> j(3) + 4l(4,5;6,7) momenta
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
       
        rm2(1) = mass_Top**2 ! mass squared input
        rm2(2) = mass_VV(2)**2 ! mass squared input

      if (leg(1).eq.-5) then ! b_bar gluon process
       lW1fromT = .false.
      else ! b gluon process
       lW1fromT = .true.
      endif
      call OSMap_TW(lW1fromT,Vext(0,3),rm2(1),rext(0,3))
      
!      do i = 1,proc_nlegs_os
!      print*,"OSMap: os parCM i =",i,rext(:,i),
!     & sqrt(abs(rext(0,i)**2-rext(1,i)**2-rext(2,i)**2-rext(3,i)**2)) 
!      enddo

! the ordering: q1 q2 > g3 V4(W+) V5(W-) > l6 l7 l8 l9:

!!![store os momenta in the partonic CMS
! oext(nlegs_os_max,4)! on-shell momenta in the partonic CMS
      do i = 1, proc_nlegs_os
      do j = 0, 3
       oext(i,j+1) = rext(j,i)
      enddo
      enddo
!!!]store os momenta in the partonic CMS

!!![store os momenta in the VV CMS
***[ boost to the X (= V1 + V2) rest frame
      do i = 0, 3
       qX(i) = rext(i,1) + rext(i,2) - rext(i,3)
      enddo
      qX_inv(0) = qX(0)
      do i = 1, 3
      qX_inv(i) = -qX(i)
      enddo
      do i = 1, proc_nlegs_os
      call BOOSTN(rext(0,i),qX_inv(0),sext(0,i))
      enddo
***] boost to the X (= V1 + V2) rest frame

! oext_vv(nlegs_os_max,4)! on-shell momenta in the VV CMS
      do i = 1, proc_nlegs_os
      do j = 0, 3
       oext_vv(i,j+1) = sext(j,i)
      enddo
!      print*,"OSMap: os VV i =",i,oext_vv(i,:),
!     & sqrt(abs(oext_vv(i,1)**2-oext_vv(i,2)**2-oext_vv(i,3)**2-oext_vv(i,4)**2))
      enddo
!      print*,"OSMap VV sum =",oext_vv(1,:)+oext_vv(2,:)
!     & -oext_vv(3,:)-oext_vv(6,:)-oext_vv(7,:)-oext_vv(8,:)-oext_vv(9,:)
!!!]store os momenta in the VV CMS

!!![store os momenta in the Lab frame
      do i = 1, proc_nlegs_os
! boost from partonic CMS to the LAB frame
      call boost0(rgam_boost, -rbetgam_boost,
     &  rext(0,i), rext(1,i), rext(2,i), rext(3,i))
      enddo
! oext_pp(nlegs_os_max,4)! on-shell momenta in the pronton-proton CMS (Lab frame)
      do i = 1, proc_nlegs_os
      do j = 0, 3
       oext_pp(i,j+1) = rext(j,i)
      enddo
      enddo
!!!]store os momenta in the Lab frame

!!!]finalizing for OS projection
        return
	end subroutine
***********************************************************************
      subroutine OSMap_TW(lW1fromT,Vext,rm2,rext)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C input: Vext(0:3,3:7) are the momenta of b + 4 leptons in 
C the rest frame of X = T + W; rm2(2) are the physical 
C masses of T and W. 
C output: rext(0:3,7) are the OS-projected momenta of the b (first entry), 
C two gauge bosons (next 2 entries) and of the 4 leptons (WpWm ordering).
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      double precision rext(0:3,7), Vext(0:3,3:7), rm2(2)
      double precision Vab(0:3,2), Vab_os(0:3,2)
      integer i
      integer t1, t2, t3 ! top-quark decays
      integer w1, w2 ! W decays
      double precision rm(2), p(0:3,3), Lep(0:3,5)
      double precision Lep_os(0:3,2), r(0:3,3), pT_inv(0:3), pT(0:3)
      double precision rm2x(2)
      logical lW1fromT ! W1 = Wp

      do i = 1, 2
      rm(i) = sqrt(rm2(i))
      enddo

C Note: Vext(7,4) are the qq -> j(3) + 4l (4,5;6,7) momenta
      t1 = 3 ! bottom
      if (lW1fromT) then
      t2 = 4 ! lepton 1
      t3 = 5 ! lepton 2
      w1 = 6 ! lepton 3
      w2 = 7 ! lepton 4
      else
      t2 = 6
      t3 = 7
      w1 = 4
      w2 = 5
      endif
!!![map TW on-shell
!      print*,"OSmap_TW t1,t2,t3 =",t1,t2,t3      
!      print*,"OSmap_TW w1, w2 =",w1,w2
            
      Vab(:,1) = Vext(:,t1)+Vext(:,t2)+Vext(:,t3)
      Vab(:,2) = Vext(:,w1) + Vext(:,w2)
      call OSMap_AB(Vab,rm2,Vab_os)
!      print*,"Vab_os T =",Vab_os(:,1)
!      print*,"Vab_os W =",Vab_os(:,2)
!      print*,"Vab_os psum =",Vab_os(:,1)+Vab_os(:,2)
!      print*,"Vab_os MT =",sqrt(Vab_os(0,1)**2-Vab_os(1,1)**2-Vab_os(2,1)**2-Vab_os(3,1)**2)
!      print*,"Vab_os MW =",sqrt(Vab_os(0,2)**2-Vab_os(1,2)**2-Vab_os(2,2)**2-Vab_os(3,2)**2)

***[projected lepton momenta (T > L1 + L2 + L3)
       p(:,1) = Vext(:,t1)
       p(:,2) = Vext(:,t2)
       p(:,3) = Vext(:,t3)
       call OSMap_decayV_NLO(p(0,1),rm(1),Vab_os(0,1),Lep(0,1))
!       print*,"Lep t1 =",Lep(:,1)
!       print*,"Lep t2 =",Lep(:,2)
!       print*,"Lep t3 =",Lep(:,3)
!       print*,"Lep sum =",Lep(:,1)+Lep(:,2)+Lep(:,3)
!       print*,"Vab_os T=",Vab_os(:,1)
***]projected lepton momenta (T > L1 + L2 + L3)

***[projected lepton momenta (V2 > L4 + L5)
       call OSMap_decayV_LO(Vext(0,w1),rm(2),Vab_os(0,2),Lep(0,4))
***]projected lepton momenta (V2 > L4 + L5)
!!!]map TW on-shell

!!![make the W from T-decay on-shell
!!![boost to the top-quark rest frame
      pT(:) = Lep(:,1) + Lep(:,2) + Lep(:,3)
      pT_inv(0) = pT(0)
      do i = 1, 3
      pT_inv(i) = -pT(i)
      enddo
      do i = 1, 3
      call BOOSTN(Lep(0,i),pT_inv(0),r(0,i))
      enddo
!!!]boost to the top-quark rest frame

      Vab(:,2) = r(:,2) + r(:,3) ! off-shell W from Top-decay
      Vab(:,1) = r(:,1) ! the bottom
      rm2x(1) = 0d0
      rm2x(2) = rm2(2)
      
      call OSMap_AB(Vab,rm2x,Vab_os)
      Lep(:,1) = Vab_os(:,1) ! bottom quark
!      print*,"Vab_os B  =",Vab_os(:,1)
!      print*,"Vab_os Wm =",Vab_os(:,2)
!      print*,"Vab_os psum =",Vab_os(:,1)+Vab_os(:,2)
!      print*,"Vab_os MB  =",sqrt(abs(Vab_os(0,1)**2-Vab_os(1,1)**2-Vab_os(2,1)**2-Vab_os(3,1)**2))
!      print*,"Vab_os MWm =",sqrt(Vab_os(0,2)**2-Vab_os(1,2)**2-Vab_os(2,2)**2-Vab_os(3,2)**2)

***[projected lepton momenta (V2 > L2 + L3)
      call OSMap_decayV_LO(r(0,2),rm(2),Vab_os(0,2),Lep_os(0,1))
      Lep(:,2) = Lep_os(:,1)
      Lep(:,3) = Lep_os(:,2)
!      print*,"Vab_os Wm     =",Vab_os(:,2)
!      print*,"Lep Wm os sum =",Lep(:,2) + Lep(:,3)
***]projected lepton momenta (V2 > L2 + L3)
!!!]make the W from T-decay on-shell
!!![boost back to the TW CM frame
      do i = 1, 3
       call BOOSTN(Lep(0,i),pT(0),r(0,i))
       Lep(:,i) = r(:,i)
      enddo
!      print*,"pT =",Lep(:,1)+Lep(:,2)+Lep(:,3)
!      print*,"pW =",Lep(:,4)+Lep(:,5)
!!!]boost back to the TW CM frame
!!![re-ordering particles: W+W-
! the ordering: q1 q2 > g3 V4 V5 > l6 l7 l8 l9:
      rext(:,1) = Lep(:,1) ! bottom
      if (lW1fromT) then ! Wp from Top
      rext(:,2) = Lep(:,2) + Lep(:,3)
      rext(:,3) = Lep(:,4) + Lep(:,5)
      rext(:,4) = Lep(:,2)
      rext(:,5) = Lep(:,3)
      rext(:,6) = Lep(:,4)
      rext(:,7) = Lep(:,5)
      else
      rext(:,2) = Lep(:,4) + Lep(:,5) 
      rext(:,3) = Lep(:,2) + Lep(:,3)
      rext(:,4) = Lep(:,4)
      rext(:,5) = Lep(:,5)
      rext(:,6) = Lep(:,2)
      rext(:,7) = Lep(:,3)
      endif
!!!]re-ordering particles: W+W-
!      do i = 1, 7
!      print*,"OSMap_TW rext =",i, rext(:,i)
!      enddo
      END subroutine OSMap_TW
