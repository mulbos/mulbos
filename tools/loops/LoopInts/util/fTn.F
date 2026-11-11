#include "defs.h"

*#define DEBUG_fTn

      subroutine CYPfT(raa,rbb,rcc,rgg,rhh,cdd,cee,cff,cjj,
     > signf,signj,csum)
************************************************************************
* Def. fT = \int_0^1dx\int_0^xdy 
*  \fr{1}{(rgg*x+rhh*y+cjj+i*signj)(raa*x^2+rbb*y^2+rcc*xy+cdd*x+cee*y+cff+i*signf)}
* with signf=-eps, {raa,rbb,rcc,rgg,rhh} are real, {cdd,cee,cff,cjj} are complex.
* Using the Barycentric technique.
* output: csum
* important: variables "signX" is the sign of img(X) in case X becomes real.
************************************************************************
      IMPLICIT NONE      

      RVAR raa,rbb,rcc,rgg,rhh,signf,signj
      CVAR cdd,cee,cff,cjj
      CVAR csum,cresd,ctv
      CVAR cbeta1,cbeta2,cbeta
      CVAR ckk,cll,cnn,cyi(2)
      CVAR ctemp
      CVAR cbbj(6),cccj(6),ccj(6)
      CVAR cyij,cyij1,cbj(6)
      RVAR raaj(6)

      RVAR sj,reps,rsj
      RVAR signnn,signyi(2)
      RVAR signyij,signccj,signtv
      
      integer checketa(2),nyi
      RVAR ru(2),rv(2)
      
      integer i,j,k
      
#include "lt.h"

	csum = 0

#ifdef DEBUG_fTn
      print*,'fTn.F hello'
      print*,'raa =',raa
      print*,'rbb =',rbb
      print*,'rcc =',rcc

      print*,'cdd =',cdd
      print*,'cee =',cee
      print*,'cff =',cff
      print*,'signf =',signf

      print*,'rgg =',rgg
      print*,'rhh =',rhh
      print*,'cjj =',cjj
      print*,'signj =',signj
#endif

      if ( abs(rbb) .le. 0 ) then 
	call CYPfT_lin(raa,rcc,rgg,rhh,cdd,cee,cff,cjj,
     >  signf,signj,csum)
     
	return
      endif 
      
      if ( abs(raa) .le. 0 ) then
* change the integration variables to get rbb=0 as above
        
	cjj=cjj+rgg+rhh
	cff=rbb+rcc+cdd+cee+cff
	ctemp=cdd
	cdd=-cdd-cee-2D0*(rbb+rcc)
	cee=rcc+ctemp
	
	raa=rbb+rcc
	rbb=rgg
	rgg=-(rgg+rhh)
	rhh=rbb
	rcc=-rcc
	       
	call CYPfT_lin(raa,rcc,rgg,rhh,cdd,cee,cff,cjj,
     >  signf,signj,csum)
     
	return
      endif 	      	
      
***[calculate the sign of img(J) and ieps
      rsj = DIMAG(cjj)
      if ( rsj .eq. 0) rsj=signj
      rsj=rsj/abs(rsj) 
      
      reps=-DIMAG(cff)
      if (reps.eq.0) reps=-signf
      reps=reps/abs(reps)
***]calculate the sign of img(J) and ieps

***[change the sign of G,H,J
* sj=1 or -1
       sj=rsj*reps
       rgg=-rgg*sj
       rhh=-rhh*sj
       cjj=-cjj*sj
***]change the sign of G,H,J

***[calculate beta
* beta is one root of: B*beta^2+C*beta+A=0
* we do not need the ieps for beta

      ctemp = rcc**2 - 4*rbb*raa
      cbeta1=(-rcc-sqrt(ctemp))/2/rbb
      cbeta2=(-rcc+sqrt(ctemp))/2/rbb
      if (abs(cbeta1).gt.abs(cbeta2)) then 
        cbeta2=raa/rbb/cbeta1
      else
	cbeta1=cbeta2      
        cbeta2=raa/rbb/cbeta2
      endif
***]calculate beta

#ifdef DEBUG_fTn
      print*,'cbeta1, cbeta2 =',cbeta1, cbeta2
#endif
      
*** cbeta(i) cannot be 0 or 1 here ***
* this is a special case considered above.
c      if (abs(cbeta2).lt.tiny .or. abs(cbeta1).lt.tiny 
c     > .or. abs(1D0-cbeta2).lt.tiny .or. abs(1D0-cbeta1).lt.tiny ) then
c        print*,'warning from fTn.F: cbeta close to 0,1'
c	print*,'cbeta1, cbeta2 =',cbeta1,cbeta2
c	print*,'something wrong with rotation?'
c      endif
      
     
*** which one for beta ? ***
      if (abs(log(abs(cbeta1))).gt.abs(log(abs(cbeta2)))) then
        ctemp=cbeta1
	cbeta1 = cbeta2
	cbeta2=ctemp
      endif 	            	  

***[look at the prefactor 1/(S*V-T*U)
* eq. (S*V-T*U)=K*y^2+L*y+N=0
* The two roots are calculated.
* To check if the roots are inside the beta triangle.
* Leading Landau Sing. can occur if y1=y2 and eps->0

      cbeta=cbeta1

      ckk=rbb*(rgg+rhh*cbeta) - rhh*(rcc+2*rbb*cbeta)

      cll=(rgg+rhh*cbeta)*cee-rhh*(cdd+cee*cbeta)
     > -cjj*(rcc+2*rbb*cbeta)
      cnn=(rgg+rhh*cbeta)*cff-cjj*(cdd+cee*cbeta)
      
* the ieps for cnn
      signnn=signf*DBLE(rgg+rhh*cbeta)
      if (signnn.eq.0) signnn=signf      

      checketa(1)=0
      checketa(2)=0

      if (abs(ckk).eq.0) then
       if (abs(cll).ne.0) then
        nyi=1
	cyi(1)=-cnn/cll
* ieps for this pole
        signyi(1)=-signnn*DBLE(cll)
	if (signyi(1).eq.0) signyi(1)=signf 	
	
       else
        nyi=0
	if (abs(cnn).eq.0) then 
	if ( runmode .ge. 2 ) then
	 print*,'Error in CYPfT(): the prefactor is 1/0, return 0'
	endif
	 return
	endif
	
	goto 100
       endif
       
      else
       nyi=2
       call CYPRoot2(ckk, cll, cnn, signnn, cyi, signyi)
      endif

      if (DIMAG(cbeta) .eq. 0) goto 100

***{calculate the parameters to say if the poles are inside the triangle
	do i=1,nyi
	 call CYPtriangle(cbeta,cyi(i),ru(i),rv(i))
	enddo  
***}calculate the parameters to say if the poles are inside the triangle

* to check if the poles are inside the beta triangle? 
      do i=1,nyi
	if (ru(i).ge.0.and.rv(i).ge.0.and.(ru(i)+rv(i)).le.1) then 
	  checketa(i)=1
	else
	  checketa(i)=0
	endif
      enddo
***]look at the prefactor 1/(S*V-T*U)     

***{ we need: checketa(i),cbeta,nyi,cyi,signyi,ckk,cll,cnn}

 100  continue

#ifdef DEBUG_fTn
      print*,'cbeta1, cbeta2 =',cbeta1, cbeta2
      print*,'cbeta =',cbeta

      print*,'ckk, cll, cnn =',ckk, cll, cnn
      print*,'Number of roots: nyi =',nyi
      print*,'Roots: cyi =',cyi

      print*,'Beta-triangle check: checketa =',checketa
#endif
 
***{ the coefficients of the 6 log arguments
      raaj(1)=0D0
      cbbj(1)=rhh
      cccj(1)=rgg+cjj
      
      raaj(2)=0D0
      cbbj(2)=rgg+rhh
      cccj(2)=cjj
      
      raaj(3)=0D0
      cbbj(3)=rgg
      cccj(3)=cjj
      
      raaj(4)=rbb
      cbbj(4)=rcc+cee
      cccj(4)=raa+cdd+cff
      
      raaj(5)=raa+rbb+rcc
      cbbj(5)=cee+cdd
      cccj(5)=cff
      
      raaj(6)=raa
      cbbj(6)=cdd
      cccj(6)=cff
      
      
* the ieps for the log arguments
      signccj=-reps      

***} the coefficients of the 6 log arguments

***{ the cj(6)-coefficients before the logs
      ccj(1) = 1
      ccj(2) = -1 + cbeta
      ccj(3) = -cbeta
      ccj(4) = -1
      ccj(5) = 1-cbeta
      ccj(6) = cbeta
***} the cj(6)-coefficients before the logs

***{ the case nyi=0, (SV-TU)=N=constant
* no extra term is needed

      if (nyi.eq.0) then

      do 1 j=1,6
      
      if (abs(ccj(j)) .ne. 0) then 
        call CYPfS2(raaj(j),cbbj(j),cccj(j),signccj,ctemp)
        csum=csum+ctemp*ccj(j)
      endif
       	
 1    continue
* cnn cannot be zero here     
      csum=-csum/cnn/sj
      return
      endif

***} the case nyi=0, (SV-TU)=N=constant

***{ cbj(6)-coefficients of cj/(aj*y-bj-yi)
      cbj(1)=cbeta
      cbj(2)=czero
      cbj(3)=czero
      cbj(4)=cbeta
      cbj(5)=czero
      cbj(6)=czero
***} cbj(6)-coefficients of cj/(aj*y-bj-yi)

***{ the case nyi=1, (SV-TU)=L*y+N
* extra term is needed

      if (nyi.eq.1) then
      i=1
      
      if (checketa(i).eq.0) then
* no extra term      

      do j=1,6
      
      if (abs(ccj(j)) .ne. 0) then 
        cyij=-(-1)**j*(cyi(i)+cbj(j))/ccj(j)
        cyij1=1-cyij           
        signyij=-(-1)**j*signyi(i)*DBLE(ccj(j))
	if (signyij.eq.0) signyij=signyi(i)
	
        call CYPfS3n(cyij,cyij1,raaj(j),cbbj(j),cccj(j),
     > 	signccj,signyij,ctemp)
     
        csum=csum-ctemp*(-1)**j
      endif
      
      enddo
      
      else
* extra term needed 

***{calculate the residue
      cresd = 0
      
* the denominator was checked above in fS3.F therefore the (V/T)_pole 
* should be safe now:	
	ctv=rbb*cyi(i)**2+cee*cyi(i)+cff
	ctv=(rhh*cyi(i)+cjj)/ctv
	
	ctemp=(rcc+2*cbeta*rbb)*cyi(i)+cdd+cee*cbeta
	ctemp=(rgg+cbeta*rhh)/ctemp
	
        if (abs(DIMAG(ctemp)).gt.abs(DIMAG(ctv))) ctv=ctemp		

* if im(ctv)=0 then take the ieps from T/V
	signtv=-signf*DBLE(rhh*cyi(i)+cjj)
	if (signtv.eq.0) signtv=-signf	

        ctv=CYPLOG(ctv, signtv)
	
      if (abs(ctv) .ne. 0) then
      
      do k=1,3
            
      if (abs(ccj(k)) .ne. 0) then
        cyij=-(-1)**k*(cyi(i)+cbj(k))/ccj(k)
        cyij1=1-cyij
	signyij=-(-1)**k*signyi(i)*DBLE(ccj(k))
	if (signyij.eq.0) signyij=signyi(i)
	  
        ctemp=CYPLOG(-cyij1/cyij, signyij)
	  
	cresd=cresd-(-1)**(k+i)*ctemp
      endif
      	  	     	
      enddo
      
      cresd=cresd*ctv
      
      endif 
***}calculate the residue      
      
***{calculate the main part      
      do j=1,6

      if (abs(ccj(j)) .ne. 0) then 
        cyij=-(-1)**j*(cyi(i)+cbj(j))/ccj(j)
        cyij1=1-cyij           
        signyij=-(-1)**j*signyi(i)*DBLE(ccj(j))
	if (signyij.eq.0) signyij=signyi(i)
	
        call CYPfS3n(cyij,cyij1,raaj(j),cbbj(j),cccj(j),
     > 	signccj,signyij,ctemp)
        csum=csum-ctemp*(-1)**j
      endif
      
      enddo 
***}calculate the main part            
      
      csum=csum-cresd
      
      endif     
* cll cannot be zero here
      csum=-csum/cll/sj     
      return
      endif
***} the case nyi=1, (SV-TU)=L*y+N

***{ the case nyi=2, (SV-TU)=K*y^2+L*y+N
* extra term is needed

***{if 21
      if (nyi .eq. 2) then
      
      do 5 i=1,2
      
***{if 23
      if (checketa(i).eq.0) then      
* no extra term (residue should be 0)

      do j=1,6
      if (abs(ccj(j)) .ne. 0) then 
        cyij=-(-1)**j*(cyi(i)+cbj(j))/ccj(j)
        cyij1=1-cyij           
        signyij=-(-1)**j*signyi(i)*DBLE(ccj(j))
	if (signyij.eq.0) signyij=signyi(i)
	
        call CYPfS3n(cyij,cyij1,raaj(j),cbbj(j),cccj(j),
     > 	signccj,signyij,ctemp)
     
        csum=csum-ctemp*(-1)**(j+i)
      endif
      
      enddo
      
*** if 23      
      else
* extra term needed 

***{calculate the residue
      cresd = 0
      
* the denominator was checked above in fS3.F therefore the (V/T)_pole 
* should be safe now:	
	ctv=rbb*cyi(i)**2+cee*cyi(i)+cff
	ctv=(rhh*cyi(i)+cjj)/ctv
	
	ctemp=(rcc+2*cbeta*rbb)*cyi(i)+cdd+cee*cbeta
	ctemp=(rgg+cbeta*rhh)/ctemp
	
        if (abs(DIMAG(ctemp)).gt.abs(DIMAG(ctv))) ctv=ctemp		

* if im(ctv)=0 then take the ieps from T/V	
	signtv=-signf*DBLE(rhh*cyi(i)+cjj)
	if (signtv.eq.0) signtv=-signf	

        ctv=CYPLOG(ctv, signtv)
	
      if (abs(ctv) .ne. 0) then
      
      do k=1,3
            
      if (abs(ccj(k)) .ne. 0) then
        cyij=-(-1D0)**k*(cyi(i)+cbj(k))/ccj(k)
        cyij1=1D0-cyij
	signyij=-(-1D0)**k*signyi(i)*DBLE(ccj(k))
	if (signyij.eq.0) signyij=signyi(i)
	  
        ctemp=CYPLOG(-cyij1/cyij, signyij)
	  
	cresd=cresd-(-1D0)**(k+i)*ctemp
      endif
      	  	     	
      enddo
      
      cresd=cresd*ctv
      
      endif 
***}calculate the residue 

***{calculate the main part
      do j=1,6
      
      if (abs(ccj(j)) .ne. 0) then 
        cyij=-(-1D0)**j*(cyi(i)+cbj(j))/ccj(j)
        cyij1=1D0-cyij           
        signyij=-(-1D0)**j*signyi(i)*DBLE(ccj(j))
	if (signyij.eq.0) signyij=signyi(i)
	
        call CYPfS3n(cyij,cyij1,raaj(j),cbbj(j),cccj(j),
     > 	signccj,signyij,ctemp)
        csum=csum-ctemp*(-1D0)**(j+i)
      endif
      
      enddo      
***}calculate the main part      

      csum=csum-cresd 
      
***}if 23
      endif
      
 5    continue
 
* the prefactor of Landau det.
      ctemp=ckk*(cyi(2)-cyi(1))     

      if ( abs(ctemp) .le. 0 ) then 
      if ( runmode .ge. 2 ) then
      print*,'Error in CYPfT(): zero Landau det, return 0'
      endif
      return
      endif

      csum=-csum/sj/ctemp

      endif
***}if 21      
***} the case nyi=2, (SV-TU)=K*y^2+L*y+N

#ifdef DEBUG_fTn
      print*,'fTn.F End'
#endif

      end
     
************************************************************************
      subroutine CYPtriangle(cbeta,cy,ru,rv)
************************************************************************
* to calculate the necessary parameters to tell if the 2 roots are 
* in the triangle [0,-cbeta,1-cbeta],
* based on barycentric technique.
************************************************************************
      IMPLICIT NONE

      CVAR cbeta,cy
      RVAR pab(2),pac(2),pay(2)
      RVAR ru,rv,dotyc,dotyb,dotbc,dotbb,dotcc,abc2
      integer i       
      
      
      pay(1)=DBLE(cy)+DBLE(cbeta)
      pay(2)=DIMAG(cy)+DIMAG(cbeta)
      
      pab(1)=1D0
      pab(2)=0D0
      
      pac(1)=DBLE(cbeta)
      pac(2)=DIMAG(cbeta)
      
      dotyc=0D0
      dotyb=0D0
      dotbc=0D0
      dotcc=0D0
      dotbb=1D0
      
      do i=1,2
      dotyc=dotyc+pay(i)*pac(i)
      
      dotyb=dotyb+pay(i)*pab(i)
      
      dotbc=dotbc+pab(i)*pac(i)
      
      dotcc=dotcc+pac(i)*pac(i)
      enddo
      
      abc2=DIMAG(cbeta)**2
      
      ru=(dotyc*dotbb-dotbc*dotyb)/abc2
      
      rv=(dotcc*dotyb-dotyc*dotbc)/abc2
      
      return
      end

************************************************************************
      subroutine CYPfT_lin(raa,rcc,rgg,rhh,cdd,cee,cff,cjj,
     > signf,signj,csum)
************************************************************************
* Def. fT_lin = \int_0^1dx\int_0^xdy
* \fr{1}{(rgg*x+rhh*y+cjj)(raa*x^2+rcc*xy+cdd*x+cee*y+cff+i*signf)} 
* with signf=-eps, {raa,rcc,rgg,rhh} are real, {cdd,cee,cff,cjj} are complex.
* output: csum = fT_lin
* important: variables "signX" is the sign of img(X) in case X becomes real.
************************************************************************
      IMPLICIT NONE

      CVAR cdd,cee,cff,cjj,csum
      CVAR ckk,cll,cnn,cyi(2)
      CVAR ctemp
      CVAR cbbj(4),cccj(4)
      CVAR cyij,cyij1
      RVAR raaj(4),signyi(2)

      integer i,j,nyi
      RVAR raa,rcc,rgg,rhh,signf,signj,signnn
      RVAR temp
      RVAR rsj,reps,sj,signccj
      
#include "lt.h"

      csum = 0

***{ special cases for 3-point integrals
* use Mathematica
      if ( abs(raa + rcc) .le. 0 .and. abs(raa + cdd) .le. 0 
     & .and. abs(raa - cee) .le. 0 .and. abs(cff) .gt. 0 
     & .and. abs(rgg) + abs(rhh) .le. 0 
     & .and. abs(cjj) .gt. 0 ) then

c	print*,'fTn.F: Ninh working on 3-point special case'

      if ( abs(raa) .le. 0 ) then
	csum = 1D0/(2*cff*cjj)
      else
       csum = sqrt(1-4*cff/raa)
       if ( DBLE(csum) .eq. 0 ) then
         temp = 1D0
       else
         temp = -2/DBLE(csum)*signf/raa
       endif
       csum = (CYPSp(-2D0/(-1 + csum), temp) 
     & + CYPSp(-2D0/(-1 - csum), -temp))/raa/cjj
      endif
      return
      endif

      if ( abs(rcc + cee) .le. 0 .and. abs(cff) .gt. 0
     & .and. abs(rgg) + abs(rhh) .le. 0
     & .and. abs(cjj) .gt. 0 ) then

c       print*,'fTn.F: Ninh working on 3-point special case'
      csum = sqrt(cdd**2 - 4*raa*cff)
      if ( DBLE(csum) .eq. 0 ) then
       temp = 1D0
      else
       temp = -2/DBLE(csum)*raa*signf
      endif

      ctemp = -CYPSp(2*raa/(2*raa+cdd+csum),-2*raa*temp)
     & - CYPSp(2*raa/(2*raa+cdd-csum),2*raa*temp)
     & - cI*pi*CYPLOG((cdd-csum)/(2*raa+cdd-csum),-2*raa*temp)
     & - cI*pi*CYPLOG((cdd+csum)/(2*raa+cdd+csum),2*raa*temp)

c       print*,'ctemp 1 =',ctemp

      csum = sqrt(rcc**2 + cdd**2 - 4*raa*cff - 2*rcc*(cdd + 2*cff))
      if ( DBLE(csum) .eq. 0 ) then
       temp = 1D0
      else
       temp = -2/DBLE(csum)*(raa+rcc)*signf
      endif

      csum = CYPSp(2*(raa+rcc)/(2*raa+rcc+cdd-csum),2*(raa+rcc)*temp)
     & + CYPSp(2*(raa+rcc)/(2*raa+rcc+cdd+csum),-2*(raa+rcc)*temp)
     & + cI*pi*CYPLOG((rcc-cdd+csum)/(-2*raa-rcc-cdd+csum),
     & -2*(raa+rcc)*temp)
     & + cI*pi*CYPLOG((rcc-cdd-csum)/(-2*raa-rcc-cdd-csum),
     & 2*(raa+rcc)*temp)

c       print*,'ctemp 2 =',csum
      csum = (csum + ctemp)/(rcc*cjj)
      return

      endif
***} special cases for 3-point integrals

***{ calculate the sign of img(J) and ieps
      rsj=DIMAG(cjj)
      if (rsj.eq.0) rsj=signj
      rsj=rsj/abs(rsj) 
      
      reps=-DIMAG(cff)
      if (reps.eq.0) reps=-signf
      reps=reps/abs(reps)
***} calculate the sign of img(J) and ieps

***{ change the sign of G,H,J
* sj=1 or -1
       sj=rsj*reps
       rgg=-rgg*sj
       rhh=-rhh*sj
       cjj=-cjj*sj
***} change the sign of G,H,J

***{ the prefactor 1/(S*V-T*U)
* eq. (S*V-T*U)=K*y^2+L*y+N=0
* Leading Landau Sing. can occur if y1=y2 and eps->0
* the ieps is needed for the roots

      ckk=rhh*raa-rcc*rgg
      cll=rhh*cdd-rcc*cjj-cee*rgg
      cnn=rhh*cff-cee*cjj

* the ieps for cnn      
      signnn=-reps*rhh
      if (signnn.eq.0) signnn=-reps      
      
      if (abs(ckk).eq.0) then
       if (abs(cll).ne.0) then
        nyi=1
	cyi(1)=-cnn/cll
* ieps for this pole
        signyi(1)=-signnn*DBLE(cll)
	if (signyi(1).eq.0) signyi(1)=signf 		
	
       else
        nyi=0
	if (abs(cnn).eq.0) then 
	if ( runmode .ge. 2 ) then
	 print*,'Error in CYPfTn_lin(): the prefactor is 1/0, return 0'
	endif
	 return
	endif  
       endif
       
      else
       nyi=2 
       call CYPRoot2(ckk, cll, cnn, signnn, cyi, signyi)
      endif        
       
***} the prefactor 1/(S*V-T*U)

***{ the coefficients of the 4 log arguments
      raaj(1)=rcc+raa
      cbbj(1)=cee+cdd
      cccj(1)=cff
      
      raaj(2)=0D0
      cbbj(2)=rhh+rgg
      cccj(2)=cjj
      
      raaj(3)=0D0
      cbbj(3)=rgg
      cccj(3)=cjj
      
      raaj(4)=raa
      cbbj(4)=cdd
      cccj(4)=cff
      
* the ieps is the same for all
      signccj=-reps      

***} the coefficients of the 4 log arguments

***{ the case nyi=0, (SV-TU)=N=constant

      if (nyi .eq. 0) then

      do 1 j=1,4
      
        call CYPfS2(raaj(j),cbbj(j),cccj(j),signccj,ctemp)
        csum=csum+ctemp*(-1D0)**j
       	
 1    continue     
* cnn cannot be zero here      
      csum=-csum/cnn/sj
      return
      endif

***} the case nyi=0, (SV-TU)=N=constant

***{ the case nyi=1, (SV-TU)=L*y+N

      if (nyi.eq.1) then
      i=1

      cyij=cyi(i)
      cyij1=1D0-cyij                 
      
      do 2 j=1,4
	
        call CYPfS3n(cyij,cyij1,raaj(j),cbbj(j),cccj(j),
     > 	signccj,signyi(i),ctemp)
        csum=csum+ctemp*(-1D0)**j
      
 2    continue       
* cll cannot be zero here      
      csum=-csum/cll/sj
           
      return
      endif

***} the case nyi=1, (SV-TU)=L*y+N

***{ the case nyi=2, (SV-TU)=K*y^2+L*y+N
      if (nyi.eq.2) then
      
      do 5 i=1,2
       cyij=cyi(i)
       cyij1=1D0-cyij           
      
      do 4 j=1,4

        call CYPfS3n(cyij,cyij1,raaj(j),cbbj(j),cccj(j),
     > 	signccj,signyi(i),ctemp)
     
        csum=csum+ctemp*(-1D0)**(j+i)
      
 4    continue       
 5    continue
 
* the prefactor of Landau det.
      ctemp=ckk*(cyi(2)-cyi(1))        

      if ( abs(ctemp) .le. 0 ) then
      if ( runmode .ge. 2 ) then
      print*,'Error in CYPfTn_lin(): zero Landau det, return 0'
      endif
      return
      endif

      csum=-csum/sj/ctemp
      
      return      
      endif
***} the case nyi=2, (SV-TU)=K*y^2+L*y+N

      end     

************************************************************************
      subroutine CYPfT13(raa,rcc,rgg,rhh,cdd,cee,cff,cii,
     > signf,csum)
************************************************************************
* Def. fT13 = \int_0^1dx\int_0^xdy 
* \fr{y}{(rgg*y^2+rhh*xy+cdd*x+cii*y+cff+i*signf)(raa*y^2+rcc*xy+cdd*x+cee*y+cff+i*signf)} 
* with signf=-eps, {raa,rcc,rgg,rhh} are real, {cdd,cee,cff,cii} are complex.
* output: csum = fT13
* important: variables "signX" is the sign of img(X) in case X becomes real.
************************************************************************
      IMPLICIT NONE

      CVAR cdd,cee,cff,cii,cjj,csum
      CVAR ckk,cll,cnn,cyi(2)
      CVAR ctemp
      CVAR cbbj(4),cccj(4)
      CVAR cyij,cyij1
      RVAR raaj(4),signyi(2)

      integer nyi,j,i
      RVAR raa,rcc,rgg,rhh,signf,signnn,signccj
      RVAR temp
      
#include "lt.h"

      csum = 0  

***{ the prefactor 1/(S*V-T*U)
* eq. (S*V-T*U)=K*y^2+L*y+N=0
* Leading Landau Sing. can occur if y1=y2 and eps->0
* the ieps is needed for the roots

      ckk=rhh*raa-rcc*rgg
      cll=raa*cdd+rhh*cee-cdd*rgg-rcc*cii
      cnn=(rhh-rcc)*cff+cdd*cee-cdd*cii
      
c      print*,'from fT13.F: ckk,cll,cnn=',ckk,cll,cnn
     
* the ieps for cnn
      signnn = signf*(rhh-rcc)
* if (rhh-rcc)=0 then we are at the boundary of phase space 
* and the signnn can be + or -
      if (signnn .eq. 0) signnn = signf
      
      if (abs(ckk).eq.0) then
       if (abs(cll).ne.0) then
        nyi=1
	cyi(1)=-cnn/cll
* ieps for this root
        signyi(1)=-signnn*DBLE(cll)
	if (signyi(1).eq.0) signyi(1)=signf
	
       else
        nyi=0
	if (abs(cnn).eq.0) then 
	if ( runmode .ge. 2 ) then
	 print*,'Error in CYPfT13(): the prefactor is 1/0, return 0'
	endif
	 return
	endif  
       endif
       
      else
       nyi=2
       call CYPRoot2(ckk, cll, cnn, signnn, cyi, signyi)
      endif        
***} the prefactor 1/(S*V-T*U)

***{ the coefficients of the 4 log arguments
      raaj(1)=raa
      cbbj(1)=cee+rcc
      cccj(1)=cff+cdd
      
      raaj(2)=rgg
      cbbj(2)=cii+rhh
      cccj(2)=cdd+cff
      
      raaj(3)=rgg+rhh
      cbbj(3)=cdd+cii
      cccj(3)=cff
      
      raaj(4)=raa+rcc
      cbbj(4)=cee+cdd
      cccj(4)=cff
      
* the ieps is the same for all
      signccj=signf      

***} the coefficients of the 6 log arguments

***{ the case nyi=0, (SV-TU)=N=constant
      if (nyi .eq. 0) then

      do 1 j=1,4
      
        call CYPfS2(raaj(j),cbbj(j),cccj(j),signccj,ctemp)
        csum = csum + ctemp*(-1)**j
       	
 1    continue     
* cnn cannot be zero here      
      csum=csum/cnn
      return
      endif

***} the case nyi=0, (SV-TU)=N=constant

***{ the case nyi=1, (SV-TU)=L*y+N
      if (nyi .eq. 1) then
      i=1
      
      cyij = cyi(i)
      cyij1 = 1 - cyij                 

      do 2 j=1,4
      
        call CYPfS3n(cyij,cyij1,raaj(j),cbbj(j),cccj(j),
     > 	signccj,signyi(i),ctemp)
     
        csum = csum + ctemp*(-1)**j
      
 2    continue       
* cll cannot be zero here      
      csum=csum/cll     
      return
      endif

***} the case nyi=1, (SV-TU)=L*y+N

***{ the case nyi=2, (SV-TU)=K*y^2+L*y+N
      if (nyi .eq. 2) then
      do 5 i=1,2
       cyij = cyi(i)
       cyij1 = 1 - cyij                 

      do 4 j=1,4  
      
        call CYPfS3n(cyij,cyij1,raaj(j),cbbj(j),cccj(j),
     > 	signccj,signyi(i),ctemp)
     
        csum = csum + ctemp*(-1)**(j+i)
      
 4    continue
 5    continue
 
* the prefactor of Landau det.
      ctemp = ckk*(cyi(2)-cyi(1))

      if ( abs(ctemp) .le. 0 ) then
      if ( runmode .ge. 2 ) then
      print*,'Error in CYPfT13(): zero Landau det, return 0'
      endif
      return
      endif

      csum = csum/ctemp
      
      endif
***} the case nyi=2, (SV-TU)=K*y^2+L*y+N

      return
      end
