!*************************************************************************  
!     This file is taken from the VBFNLO-2.7.1 program.
!     VBFNLO: https://www.itp.kit.edu/vbfnlo/wiki/doku.php
!     This file contains subroutines for reading input files. 
!     The input files are in the format: key = value 
!*************************************************************************  

module readinput
      integer, parameter :: maxlines=200
      character*250 lines(maxlines)
      integer n_lines
      integer iunit
      character*250 fname
      logical usedefaults
      
      character*250 inputpath, pdfpath

contains

      subroutine loadfile(filename)
!*************************************************************************  
!     Opens a textfile and reads the data in an internal buffer.
!     End-of-line comments, starting with ! or //, are removed.
!     Do only process one single file at a time.
!     DON"T FORGET TO CALL CLOSEFILE AFTERWARDS.
!*************************************************************************  
          !use globalvars, only: ldoblha
      implicit none

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
      print*,"input_file:",trim(path)
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
!*************************************************************************  
      subroutine closefile
!*************************************************************************  
!     Closes the current file.
!*************************************************************************  
      implicit none

      close(iunit)
      end subroutine
      
!*************************************************************************  
      SUBROUTINE read_Logical(key, value, default_value)
!*************************************************************************  
      !use globalvars, only: lglobalprint
      implicit none

      character*(*) key
      Logical value, default_value
      character*250 line
      integer i,le,lk,lt

      value=default_value
      if (usedefaults) RETURN
      do i=1,n_lines
         line = lines(i)
         lk=INDEX(line,key)
         if (lk.ne.0) then
            le=INDEX(line,"=")
            lt = lk+LEN(key) 
            if (((le.ne.0).and.(le.gt.lk)) &
     &           .and.verifyWhitespaces(line(1:lk-1)) &
     &           .and.(line(lt:lt).eq." ".or.line(lt:lt).eq."=".or.line(lt:lt).eq.achar(9))) then
               line=" "//line(le+1:len(line))
               read(line,*,ERR=10,END=20) value
               !if (showinput) write(*,"(T4,A,T25,A,L8)") key," = ",value
               RETURN
            endif
         endif
      enddo

 10   CONTINUE
      !if (lglobalprint) then
          print *,"Warning: cannot read value for key = ",key, &
     &      ". Using default value : ",default_value
      !endif
      RETURN

 20   CONTINUE
      print *,"Input Error in line ",i," of ",fname
      print *,""
      GOTO 10

      RETURN
      END subroutine
      
!*************************************************************************  
      SUBROUTINE read_String(key, value, default_value)
!*************************************************************************  
      !use globalvars, only: lglobalprint
      implicit none
     
      character*(*) key
      character*(*) value, default_value
      character*250 line
      integer i,le,lk,lt
     
      value=default_value
      if (usedefaults) RETURN
      do i=1,n_lines
         line = lines(i)
         lk=INDEX(line,key)
         if (lk.ne.0) then
            le=INDEX(line,"=")
            lt = lk+LEN(key) 
            if (((le.ne.0).and.(le.gt.lk)) &
     &           .and.verifyWhitespaces(line(1:lk-1)) &
     &           .and.(line(lt:lt).eq." ".or.line(lt:lt).eq."=".or.line(lt:lt).eq.achar(9))) then           
               line=" "//line(le+1:len(line))
               read(line,*,ERR=20,END=20) value
               !if (showinput) write(*,"(T4,A,T25,A,A)") key," = ",value
               RETURN
            endif
         endif
      enddo

 10   CONTINUE
      !if (lglobalprint) then
          print *,"Warning: cannot read value for key = ",key, &
     &      ". Using default value : ",default_value
      !endif
      RETURN

 20   CONTINUE
      print *,"Input Error in line ",i," of ",fname
      print *,""
      GOTO 10

      RETURN
      END subroutine
      
!*************************************************************************  
      SUBROUTINE read_Int(key, value, default_value)
!*************************************************************************  
      !use globalvars, only: lglobalprint
      implicit none
     
      character*(*) key
      
      integer value, default_value
      character*250 line
      integer i,le,lk,lt

      value=default_value
      if (usedefaults) RETURN
      do i=1,n_lines
         line = lines(i)
         lk=INDEX(line,key)
         if (lk.ne.0) then
            le=INDEX(line,"=")
            lt = lk+LEN(key) 
            if (((le.ne.0).and.(le.gt.lk)) &
     &           .and.verifyWhitespaces(line(1:lk-1)) &
     &           .and.(line(lt:lt).eq." ".or.line(lt:lt).eq."=".or.line(lt:lt).eq.achar(9))) then
               line=" "//line(le+1:len(line))
               read(line,*,ERR=20,END=20) value
               !if (showinput) write(*,"(T4,A,T25,A,I8)") key," = ",value
               RETURN
            endif
         endif
      enddo

 10   CONTINUE
      !if (lglobalprint) then
          print *,"Warning: cannot read value for key = ",key, &
     &      ". Using default value : ",default_value
      !endif
      RETURN

 20   CONTINUE
      print *,"Input Error in line ",i," of ",fname
      print *,""
      GOTO 10

      RETURN
      END subroutine
      
!*************************************************************************  
      SUBROUTINE read_Real(key, value, default_value)
!*************************************************************************  
      !use globalvars, only: lglobalprint
      implicit none
     
      character*(*) key
      real*8 value, default_value
      character*250 line
      integer i,le,lk,lt

      value=default_value
      if (usedefaults) RETURN
      do i=1,n_lines
         line = lines(i)
         lk=INDEX(line,key)
         if (lk.ne.0) then
            le=INDEX(line,"=")
            lt = lk+LEN(key) 
            if (((le.ne.0).and.(le.gt.lk)) &
     &       .and.verifyWhitespaces(line(1:lk-1)) &
     &       .and.(line(lt:lt).eq." ".or.line(lt:lt).eq."=".or.line(lt:lt).eq.achar(9) )) then
               line=" "//line(le+1:len(line))
               read(line,*,ERR=20,END=20) value       
               !if (showinput) write(*,"(T4,A,T25,A,G12.6)") key," = ",value
               RETURN
            endif
         endif
      enddo

 10   CONTINUE
      print *,"Warning: cannot read value for key = ",key, &
     &          ". Using default value : ",default_value
      RETURN

 20   CONTINUE
      print *,"Input Error in line ",i," of ",fname
      print *,""
      GOTO 10

      END subroutine
      
!*************************************************************************  
      logical FUNCTION verifyWhitespaces(s)
!*************************************************************************  
      implicit none
      integer i
      character*(*) s
      do i=1,len_trim(s)
        if ((s(i:i).ne." ").and.(s(i:i).ne.achar(9))) then
           verifyWhitespaces = .false.
           return
        endif
      enddo
      verifyWhitespaces = .true.
      END function


end module
