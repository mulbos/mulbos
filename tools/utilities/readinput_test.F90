!*************************************************************************  
!     This file is taken from the VBFNLO-2.7.1 program.
!     VBFNLO: https://www.itp.kit.edu/vbfnlo/wiki/doku.php
!     This file contains subroutines for reading input files. 
!     The input files are in the format: key = value 
!*************************************************************************  


      subroutine loadfile(filename)
!*************************************************************************  
!     Opens a textfile and reads the data in an internal buffer.
!     End-of-line comments, starting with ! or //, are removed.
!     Do only process one single file at a time.
!     DON"T FORGET TO CALL CLOSEFILE AFTERWARDS.
!*************************************************************************  
          !use globalvars, only: ldoblha
      implicit none
      
      integer, parameter :: maxlines=200
      character*250 lines(maxlines)
      integer n_lines
      integer iunit
      character*250 fname
      logical usedefaults
      
      character*250 inputpath, pdfpath


      character*(*) filename
     
!#include "global.inc"
!#include "BLHAhelper.inc"
      
      integer i,le,k,j
      character*250 line

      !character*250 GetInputPath
      character*250 path
      !external GetInputPath

      iunit=10
      usedefaults=.false.
      
! add input file path to filename
      !path=GetInputPath()
      !path = trim(path)//"/"//filename
      path=filename
      
      print*,"path =",path
! Open file
      open (unit=iunit, &
     &      file=trim(path), &
     &      status="old",err=77)

      fname = filename
      
! read text lines
      n_lines=0
 50   read(iunit,"(a250)",end=99,err=98) lines(n_lines+1)
      n_lines=n_lines+1
      if (n_lines.lt.maxlines) GOTO 50
 98   CONTINUE
      print *," Read error on file",filename
      print *," Using default values."
      usedefaults=.true.
 99   CONTINUE
      n_lines=n_lines+1

! remove comments
      do i=1,n_lines
         line=lines(i)
         if (INDEX(line,"!").ne.0) then
            line=line(1:INDEX(line,"!")-1)
         endif
         if (INDEX(line,"//").ne.0) then
            line=line(1:INDEX(line,"//")-1)
         endif
! upper case everything in front of a "="
         le = index(line,"=")
         if (le.gt.0) then
            do j=1,le-1
               k=ichar(line(j:j))
               if(k.ge.97.and.k.le.122) then
                  k=ichar(line(j:j))-32   
                  line(j:j)=char(k)
               endif
            enddo
         endif
         lines(i)=line
      enddo

      return
 77   CONTINUE
      print *," Open error on file ",filename
      print *," Using default values."
      usedefaults=.true.
      end subroutine

