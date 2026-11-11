#include "defs.h"

*#define DEBUG_fS3n

      subroutine CYPfS3n(cy,cy1,raa,cbb,ccc,signc,signy0,cS3)
************************************************************************
* Def. fS3n = \int_0^1 \fr{ln(raa*y^2+cbb*y+ccc+i*signc)}{y-cy} dy
* input: cy, cy1 = 1-cy, raa (real), cbb, ccc.
* signc = sign(img(ccc)), signy0 = sign(img(cy)) in case they are real.
* output: cS3 = fS3n.
* remark: ieps is needed.
************************************************************************
      IMPLICIT NONE

      CVAR cy,cy1,cbb,ccc,cS3
      CVAR ctemp
      CVAR czi(2)
      RVAR signzi(2)
      CVAR cza1,cdyza,czb1,cdyzb
      CVAR csigny1,csigny2
     
      CVAR ceta
      RVAR raa,signc,signy0,reps
      RVAR temp
      
#include "lt.h"

      cS3 = 0

#ifdef DEBUG_fS3n
      print*,'fS3n.F hello'
      print*,'fS3n.F input =',cy,cy1,raa,cbb,ccc,signc,signy0
#endif
      
***{check for end-point sing.
* some cases are very regular when y=0,1
* more work needed

c      if (abs(cy).lt.tiny.or.abs(cy1).lt.tiny) then
c       print*,'warning in fS3n.F: poles are close to end-points'
c       print*,'fS3n.F: cy =',cy
c       print*,'raa,cbb,ccc =',raa,cbb,ccc
c      endif      
***}check for end-point sing.

      if ( raa .eq. 0 ) then 
        call CYPfS3n_lin(cy,cy1,cbb,ccc,signc,signy0,cS3)
	return
      endif
      
      if (abs(cbb).eq.0 .and. abs(ccc).eq.0) then
      cS3 = raa
      cS3 = CYPLOG(cS3, signc)*CYPLOG(-cy1/cy, signy0)
       
      call CYPfS3n_lin(cy,cy1,DCMPLX(1D0),ccc,signc,signy0,ctemp)
      cS3=cS3+ctemp
      call CYPfS3n_lin(cy,cy1,DCMPLX(1D0),ccc,-signc,signy0,ctemp)
      cS3=cS3+ctemp            
      
c      print*,'new cS3=',cS3
      return
      endif 

***{reps is infinitesimal, having opposite sign to the 
* imaginary part of the argument
      reps = -DIMAG(ccc)
      if (reps .eq. 0) reps=-signc

***}reps is infinitesimal, having opposite sign to the 
* imaginary part of the argument

***{calculate the 2 roots
* eq.: y**2+(b/a)y+(c/a)=0

      call CYPRoot2(cone, cbb/raa, ccc/raa, -reps/raa, czi, signzi)

***}calculate the 2 roots

***{calculate intermediate variables
      cza1 = 1 - czi(1)
      cdyza = cy - czi(1)

      czb1 = 1 - czi(2)
      cdyzb = cy - czi(2)
***}calculate intermediate variables

***{calculate the eta function
* the sign of the imaginary parts

       temp = DIMAG(czi(1))
       if (temp .eq. 0) then
       csigny1 = DCMPLX(0D0, signzi(1)) 
       else
       csigny1 = DCMPLX(0D0, temp)
       endif 
       
       temp = DIMAG(czi(2))
       if (temp .eq. 0) then
       csigny2 = DCMPLX(0D0,signzi(2))
       else
       csigny2=DCMPLX(0D0, temp)
       endif       
       
       temp = DIMAG(ccc/raa)
       if (temp .eq. 0) then
       ctemp=DCMPLX(0D0, signc/raa)
       else
       ctemp=DCMPLX(0D0,temp)       
       endif
       
       ceta = CYPetax(-csigny1, -csigny2, ctemp)        
***}calculate the eta function

***{calculate the log
      cS3=raa
      cS3=CYPLOG(cS3, -reps)
      
      ctemp=CYPLOG(-cy1/cy, signy0)
***}calculate the log

       cS3=(ceta + cS3)*ctemp
       
***{calculate the 2 R-functions
      call CYPfRn(cy,cy1,czi(1),cza1,cdyza,signzi(1),signy0,ctemp)
      cS3=cS3+ctemp      
      
      call CYPfRn(cy,cy1,czi(2),czb1,cdyzb,signzi(2),signy0,ctemp)
      
      cS3=cS3+ctemp 
***}calculate the 2 R-functions        

#ifdef DEBUG_fS3n
      print*,'fS3n.F End'
#endif
       
      return
      end
      
************************************************************************
      subroutine CYPfS3n_lin(cy,cy1,cbb,ccc,signc,signy0,cS3)
************************************************************************
* Def. fS3n_lin = \int_0^1 \fr{ln(cbb*y+ccc+i*signc)}{y-cy} dy
* input: cy, cy1 = 1 - cy, cbb, ccc
* signc=sign(img(ccc)) in case ccc is real.
* output: cS3
* remark: ieps is needed.
************************************************************************
      IMPLICIT NONE

      CVAR cy,cy1,cbb,ccc,cS3
      CVAR ctemp
      CVAR cz,cz1,cdyz
      CVAR zfflog
      CVAR ceta
      
      RVAR signc,signy0,reps,signz,temp
      
#include "lt.h"

      cS3 = 0
      
      if (abs(cbb) .le. 0) then
        if (abs(ccc) .le. 0) then 
	  print*,'Error in CYPfS3n_lin(): log(0) encountered, return 0'
	  return
	endif  
	
	cS3 = CYPLOG(ccc, signc)*CYPLOG(-cy1/cy, signy0)
	return
      endif

***{reps is infinitesimal, having opposite sign to the 
* imaginary part of the argument
      reps = -DIMAG(ccc)
      if ( reps .eq. 0 ) reps=-signc
      
***}reps is infinitesimal, having opposite sign to the 

***{calculate the root: cz=-c/b
* eq.: by+c=0
      cz=-ccc/cbb
      
* calculate the sign of img(cz) related to ieps
      signz=-signc*DBLE(cbb)
      if (signz.eq.0) signz=-signc

***}calculate the root: cz=-c/b

      cz1 = 1 - cz
      cdyz = cy - cz   

***{calculate the eta functions

* for real b
      ceta = 0

* for complex b
       if ( DIMAG(cbb) .ne. 0 ) then
       ctemp=DCMPLX(0D0,-reps)
       temp = -DIMAG(cz)
       if (temp .eq. 0) temp = -signz  
       ceta = CYPetax(cbb, DCMPLX(0D0,temp), ctemp)
       endif
***}calculate the eta functions

      cS3=(CYPLOG(cbb, signc) + ceta)*CYPLOG(-cy1/cy, signy0)
             
***{calculate the R-function
      call CYPfRn(cy,cy1,cz,cz1,cdyz,signz,signy0,ctemp)
      
      cS3=cS3+ctemp
***}calculate the R-function 

      end

************************************************************************
      subroutine CYPfS2(raa,cbb,ccc,signc,cS2)
************************************************************************
* Def. S2=\int_0^1 dy \ln(raa*y^2+cbb*y+ccc), 
* where raa is real; cbb and ccc are complex.
* input: raa, cbb, ccc
* signc=sign(img(ccc)) in case ccc is real.
* output: cS2
* remark: ieps is needed.
************************************************************************
      IMPLICIT NONE

      CVAR cbb,ccc,cS2

      CVAR ctemp
      CVAR czi(2)
      CVAR csigny1,csigny2
      CVAR ceta
      RVAR signc,raa,signzi(2)
      RVAR temp
      CVAR caa

#include "lt.h"

c	print*,'fS2.F hello'

      cS2 = 0
     
      if ( abs(raa) .le. 0 ) then 
        call CYPfS2_lin(cbb,ccc,signc,cS2)
	return
      endif 	

***{calculate the 2 roots
* eq.: y**2+(b/a)y+(c/a)=0

      call CYPRoot2(cone, cbb/raa, ccc/raa, signc/raa, czi, signzi)

***}calculate the 2 roots

***{calculate the eta function
       temp=DIMAG(czi(1))
       csigny1=DCMPLX(0D0,temp)
       if (temp.eq.0) csigny1=DCMPLX(0D0,signzi(1)) 
       
       temp=DIMAG(czi(2))
       csigny2=DCMPLX(0D0,temp)
       if (temp.eq.0) csigny2=DCMPLX(0D0,signzi(2))
       
       temp=DIMAG(ccc/raa)
       ctemp=DCMPLX(0D0,temp)
       if (temp.eq.0) ctemp=DCMPLX(0D0,signc/raa)
       
       ceta=CYPetax(-csigny1,-csigny2,ctemp)
***}calculate the eta function

***{calculate the log
       temp=DIMAG(ccc)
       if (temp.eq.0) temp=signc
       caa = raa
       ctemp=CYPLOG(caa, temp)
       
       cS2 = ctemp + ceta
***}calculate the log      

      call CYPfS2_linr(czi(1),signzi(1),ctemp)
      cS2=cS2+ctemp
      call CYPfS2_linr(czi(2),signzi(2),ctemp)
      cS2=cS2+ctemp

      end
      
************************************************************************
      subroutine CYPfS2_lin(cbb,ccc,signc,cS2)
************************************************************************     
* Def. fS2_lin = \int_0^1 dy \ln(cbb*y + ccc), 
* where cbb and ccc are complex.
* input: cbb, ccc.
* signc=sign(img(ccc)) in case ccc is real.
* output: cS2
* remarks: ieps is needed.
************************************************************************
      IMPLICIT NONE

      CVAR cbb,ccc,cS2

      CVAR csignz
      CVAR cz
      RVAR temp,signc,signz,reps
      
#include "lt.h"

      cS2 = 0
      
      if (abs(cbb) .le. 0) then
        if (abs(ccc) .le. 0) then 
	  print*,'Error in CYPfS2_lin(): log(0) encountered, return 0'
	  return
	else
          cS2 = CYPLOG(ccc, signc)
	  return     
        endif
      endif 	

      reps = -DIMAG(ccc)
      if (reps.eq.0) reps = -signc
      
***{ compute cz
      cz = -ccc/cbb
      signz = -signc*DBLE(cbb)
      if (signz .eq. 0) signz = -signc
***} compute cz      
      
***{ for real-b
      if (DIMAG(cbb) .eq. 0) then 
       call CYPfS2_linr(cz,signz,cS2)
    
       cS2 = cS2 + CYPLOG(cbb, -reps)
       return
      endif        
***} for real-b

***{ for complex-b      
      call CYPfS2_linr(cz,signz,cS2)

      temp = DIMAG(cz)
      if (temp .eq. 0) temp = signz
      csignz = DCMPLX(0D0,temp)
      
      cS2 = cS2 + CYPetax(cbb, -csignz, DCMPLX(0D0, -reps)) + log(cbb)
***} for complex-b      

      end     
      
************************************************************************
      subroutine CYPfS2_linr(cz,signz,cS2)
************************************************************************     
* Def. fS2_linr = \int_0^1 dy ln(y - cz), 
* input: cz, signz=sign(img(cz)) in case cz is real.
* output: cS2
* remark: ieps is needed.
************************************************************************
      IMPLICIT NONE

      CVAR cz,cS2
      RVAR signz
      
#include "lt.h"

      cS2 = 0
      
***{ special cases
      if (abs(cz) .le. 0) then 
        cS2 = -1D0
	return
      endif
       
      if (abs(cz - 1) .le. 0) then
        cS2 = CYPLOG(-cone, -signz) - 1
	return
      endif
***} special cases

      cS2 = CYPLOG(1-cz, -signz)*(1-cz) + cz*CYPLOG(-cz, -signz) - 1
      end

************************************************************************
      subroutine CYPfRn(cy,cy1,cza,cza1,cdyza,signza,signy,cR3)
************************************************************************
* Def. Rn = \int_0^1 \fr{x-cza-i*signza}{x-cy-i*signy} dx
* Input: cy,cy1=1-cy,cza,cza1=1-cza,cdyza=cy-cza,signza,signy
* Output: cR3
* i*sign = -i*eps is needed.
************************************************************************
      IMPLICIT NONE

      CVAR cy,cy1,cza,cza1,cdyza,cR3
      CVAR ctemp
      RVAR signza,signy

#include "lt.h"
       
      cR3 = 0

***{cy=cza and check for singularities
*   be careful with log(0) singularity.
      if (abs(cdyza) .le. 0) then
      
      if (DIMAG(cy).eq.0.and.signy*signza.lt.0
     > .and.DBLE(cy).ge.0) then
     
        if (DBLE(cy).le.1) then
        print*,'Error in CYPfRn(): 
     &	pinch singularity encountered, return 0'
c          print*,'cy,cza =',cy,cza
	return

***{3pt-Landau sing. formulae: Eq.(B.18) of arXiv:0902.0325	  
c          ctemp=-signy
c          cR3=DCMPLX(0D0,-4D0*pi*signza/abs(signza))
c          cR3=cR3*(CYPLOG(-DCMPLX(dzero),1,ctemp, ier)
c     >	  -CYPLOG(-cy,1,ctemp,ier))
	  
c          cR3=cR3+CYPLOG(-cy1/cy,1,-ctemp,ier)
c     >       *(CYPLOG(cy1,1,ctemp,ier)+CYPLOG(-cy,1,ctemp,ier))
     
c          cR3=cR3/2D0
***}3pt-Landau sing. formulae: Eq.(B.18) of arXiv:0902.0325

	else
          cR3=DCMPLX(0D0,-4*pi*signza/abs(signza))
          
	  cR3=CYPLOG(-cy1/cy, signy)
     >    *(CYPLOG(cy1, -signy) + CYPLOG(-cy, -signy)+cR3)/2
	endif
	
      else
       cR3=CYPLOG(-cza1/cza, signza)
     >  *(CYPLOG(cza1, -signza)+CYPLOG(-cza, -signza))/2
      endif
      return
      endif
***}cy=cza and check for singularities

      cR3 = CYPLOG(cdyza, signy)
      if ( abs(cR3) .gt. 0 ) cR3 = cR3*CYPLOG(-cy1/cy, signy)
      
      call CYPfR(cy,cy1,cza,cza1,cdyza,signza,signy,ctemp)
      
      cR3 = cR3 + ctemp
      end

************************************************************************
      subroutine CYPfR(cy,cy1,cza,cza1,cdyza,signza,signy,cR3)
************************************************************************      
* Def. fR = \int_0^1 {ln(y - cza) - ln(cy - cza)}/(y - cy)
* input: cy, cy1=1-cy, cza, cza1=1-cza, cdyza=cy-cza. 
* signza=sign(img(cza)) in case it is real.
* output: cR3
************************************************************************
      IMPLICIT NONE

      CVAR cy,cy1,cza,cza1,cdyza,cR3
      CVAR csigny1,csigndy01,csignab1,csignab2
      CVAR ceta1,ceta2

      RVAR signza,signy
      RVAR temp

#include "lt.h"
       
      cR3 = 0
      
***{Special cases
      if (abs(cza) .le. 0) then   
      
c      print*,'fR.F: special cases'

      cR3 = -CYPLOG(cy, signy)*CYPLOG(-cy1/cy, signy)
      cR3 = cR3 + CYPSp(1/cy, -signy)

c      print*,'new new cR3=',cR3,cy
      return
      endif
      
      
***}Special cases      

***{calculate the SIGN of imaginary parts and eta functions
* we do not need the ieps for y0
* if img(y0)=img(y1) we may need the ieps for the logs

       temp=DIMAG(cza)

       csigny1=DCMPLX(0D0,temp)
       if (temp.eq.0) csigny1=DCMPLX(0D0,signza) 
       
       temp=DIMAG(cdyza)
      
       csigndy01=DCMPLX(0D0,temp)
       if (temp.eq.0) csigndy01=DCMPLX(0D0,signy)

       temp=DIMAG(-cza/cdyza)
       if (temp.eq.0) temp = DBLE(cza)*signy
*** this step: not checked but same as below ***
* choose +signy since this ieps is relevant if cza in (0,1)      
       if (temp.eq.0) temp = signy       
       csignab1=DCMPLX(0D0,temp)
       
      ceta1 = CYPetax(-csigny1,-csigndy01,csignab1)
       
      temp=DIMAG(cza1/cdyza)

       if (temp.eq.0) temp = -DBLE(cza1)*signy

*** this step: checked and worked ***
* choose -signy since this ieps is relevant if cza in (0,1)                    
       if (temp.eq.0) temp = -signy  	
       
       csignab2 = DCMPLX(0D0,temp)
       
      ceta2 = CYPetax(-csigny1,-csigndy01,csignab2)       
***}calculate the SIGN of imaginary parts and eta functions

***{calculate R-func from Sp-func
* Def. R(y0,y1)=Sp(y0/(y0-y1))-Sp((y0-1)/(y0-y1))+ln(y0/(y0-y1))*eta(-y1,1/(y0-y1))
*              -ln((1-y0)/(y1-y0))*eta(1-y1,1/(y0-y1))

      if (abs(ceta1) .ne. 0) cR3 = CYPLOG(cy/cdyza, -DIMAG(csignab1))*ceta1
      
      if (abs(ceta2) .ne. 0) cR3=cR3
     > - CYPLOG(-cy1/cdyza, -DIMAG(csignab2))*ceta2
      
      cR3 = cR3 + CYPSp(cy/cdyza, -DIMAG(csignab1))
     & - CYPSp(-cy1/cdyza, -DIMAG(csignab2))

***}calculate R-func from Sp-func

      end
      
