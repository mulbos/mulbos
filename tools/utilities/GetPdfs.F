* This code links to the LHAPDF library.
* leg(1), leg(2) are the PDG parton codes

************************************************************************
	subroutine InitPdfs
	use global_vars
	implicit none

	integer partonsign
	common /hadron/ partonsign

cc for p anti-p colliders (like the Tevatron)
c	  partonsign = -1

cc for p p colliders (like the LHC)
	  partonsign = 1

      if (lLHAPDF) then
	call InitPDFsetByName(pdf_set)
	call InitPDF(pdf_mem)
      endif

	end
************************************************************************
	subroutine MultiplyPdfs(fac)
	use global_vars
	use global_mom_xpdf
	implicit none
	double precision fac, pdfs

	double precision xPDF(-6:6), tPDF(-6:6), phot
	double precision tau, x, xrange
	integer fail, isgn(2), i
	integer parton1, parton2
	integer partonsign
	common /hadron/ partonsign
        logical has_photon

	tau = tau_pdf**2
	xrange = 1 - tau
	x = tau + xrange*x_pdf

	if( x .lt. 1D-10 .or. x .gt. (1 - 1D-10) ) then
	  fac = 0
	  return
	endif
	
        if (has_photon()) then
        call evolvePDFphoton(x, scale, xPDF, phot)
        xPDF(-6) = phot
        xPDF(6) = phot
        call evolvePDFphoton(tau/x, scale, tPDF, phot)
        tPDF(-6) = phot
        tPDF(6) = phot
        else
        call evolvePDF(x, scale, xPDF)
        call evolvePDF(tau/x, scale, tPDF)
C Ninh check
c        do i = -6, 6
c         xPDF(i) = x
c         tPDF(i) = tau/x
c        enddo
        endif

***[convolute with PDfs, sum over light quarks
        parton1 = leg(1)
        parton2 = leg(2)
        ! gluon
        if (parton1.eq.21) then
         parton1 = 0
        endif
        if (parton2.eq.21) then
         parton2 = 0
        endif
        ! photon
        if (parton1.eq.22) then
         parton1 = 6
        endif
        if (parton2.eq.22) then
         parton2 = 6
        endif
        
	if ( parton1 .eq. 0 .or. abs(parton1) .eq. 6 ) then 
	 isgn(1) = 0
	elseif ( abs(parton1) .lt. 3 ) then
	 isgn(1) = sign(1, parton1)
	else
	 isgn(1) = 10
	endif 

	if ( parton2 .eq. 0 .or. abs(parton2) .eq. 6 ) then 
	 isgn(2) = 0
	elseif ( abs(parton2) .lt. 3 ) then
	 isgn(2) = sign(1, parton2)
	else
	 isgn(2) = 10
	endif 

	fac = fac*xrange/(tau_pdf*x)
	if ( isgn(1) .ne. 10 .and. isgn(2) .ne. 10 
     &	 .and. abs(parton1) + abs(parton2) .ne. 0 
     &	 .and. abs(parton1) + abs(parton2) .ne. 12 
     &	 .and. abs(parton1) + abs(parton2) .ne. 6 ) then 
         pdfs = tPDF(parton1*partonsign)*xPDF(parton2) +
     &      xPDF(parton2*partonsign)*tPDF(parton1) +
     &      tPDF((parton1 + 2*isgn(1))*partonsign)*
     &      xPDF(parton2 + 2*isgn(2)) +
     &      xPDF((parton2 + 2*isgn(2))*partonsign)*
     &      tPDF(parton1 + 2*isgn(1))
         if (lZZ.and.abs(parton1).eq.1.and.lFLV5) then
          pdfs = pdfs + tPDF((parton1 + 4*isgn(1))*partonsign)*
     &      xPDF(parton2 + 4*isgn(2)) +
     &      xPDF((parton2 + 4*isgn(2))*partonsign)*
     &      tPDF(parton1 + 4*isgn(1))
         endif
         fac = fac*pdfs
	else
         fac = fac*( tPDF(parton1*partonsign)*xPDF(parton2) +
     &      xPDF(parton2*partonsign)*tPDF(parton1) )
	endif
***]convolute with PDfs, sum over light quarks


      if (parton1 .ne. parton2) then
         fac = fac*2
      endif

      end subroutine
************************************************************************
	subroutine MultiplyPdfs_coll(fac1, fac2, fac3)
	use global_vars
	use global_mom_xpdf
	implicit none
	double precision fac1, fac2, fac3
	integer pcoll

	double precision xPDF(-6:6), tPDF(-6:6), phot, temp
	double precision tau, x, y, xrange 
	integer fail, isgn(2)
	integer parton1, parton2
	integer partonsign, i
	common /hadron/ partonsign
        logical has_photon

        fac1 = 1
	fac2 = 1
	fac3 = 1
	
	if (lNLO_QCD) then
	 pcoll = 0 ! gluon
	else
	 pcoll = 6 ! photon
	endif

	tau = tau_pdf**2
	xrange = 1 - tau
	x = tau + xrange*x_pdf

***[ momentum densities	  
        if (has_photon()) then
        call evolvePDFphoton(x, scale, xPDF, phot)
        xPDF(-6) = phot
        xPDF(6) = phot
        call evolvePDFphoton(tau/x, scale, tPDF, phot)
        tPDF(-6) = phot
        tPDF(6) = phot
        else
        call evolvePDF(x, scale, xPDF)
        call evolvePDF(tau/x, scale, tPDF)
        endif
***] momentum densities	  

***[convolute with PDfs, sum over light quarks
        parton1 = leg(1)
        parton2 = leg(2)
        ! gluon
        if (parton1.eq.21) then
         parton1 = 0
        endif
        if (parton2.eq.21) then
         parton2 = 0
        endif
        ! photon
        if (parton1.eq.22) then
         parton1 = 6
        endif
        if (parton2.eq.22) then
         parton2 = 6
        endif

	if ( parton1 .eq. 0 ) then 
	 isgn(1) = 0
	else
	 isgn(1) = sign(1, parton1)
	endif 

	if ( parton2 .eq. 0 ) then 
	 isgn(2) = 0
	else
	 isgn(2) = sign(1, parton2)
	endif 

* include the Jacobian
	temp = 2*xrange/(tau_pdf*x)
* include the densities
	if ( abs(parton1) .lt. 3 .and. abs(parton2) .lt. 3 
     &	 .and. abs(parton1) + abs(parton2) .gt. 0 ) then 
	fac1 = tPDF(parton1)*xPDF(parton2) + 
     &	 tPDF(parton1 + 2*isgn(1))*xPDF(parton2 + 2*isgn(2))
        if (lZZ.and.abs(parton1).eq.1.and.lFLV5) then
        fac1 = fac1 + tPDF(parton1 + 4*isgn(1))*xPDF(parton2 + 4*isgn(2))
        endif
        fac1 = fac1*temp
	else
	fac1 = temp*tPDF(parton1)*xPDF(parton2)	
	endif

      if (parton1 .ne. parton2) then
        fac1 = fac1*2
      endif	

* watch the order, this is qbar glu
* Warning: pcoll must be the gluon(0) or photon (6)
	if ( abs(parton1) .lt. 3 .and. abs(parton1) .gt. 0 ) then
	fac2 = tPDF(parton1)*xPDF(pcoll) + 
     &   tPDF(parton1 + 2*isgn(1))*xPDF(pcoll)
        if (lZZ.and.abs(parton1).eq.1.and.lFLV5) then
        fac2 = fac2 + tPDF(parton1 + 4*isgn(1))*xPDF(pcoll)
        endif
        fac2 = fac2*temp
	else
	fac2 = temp*tPDF(parton1)*xPDF(pcoll)
        endif
      if (parton1 .eq. pcoll) then
        print*,'error in lumi_hadron.F: invalid induced process'
	return
      else
	fac2 = fac2*2
      endif

* watch the order, this is glu q
* Warning: pcoll must be the gluon(0) or photon (6)
	if ( abs(parton2) .lt. 3 .and. abs(parton2) .gt. 0 ) then
	fac3 = tPDF(pcoll)*xPDF(parton2) + 
     &   xPDF(parton2 + 2*isgn(2))*tPDF(pcoll)
        if (lZZ.and.abs(parton2).eq.1.and.lFLV5) then
        fac3 = fac3 + xPDF(parton2 + 4*isgn(2))*tPDF(pcoll)
        endif
        fac3 = fac3*temp
	else
	fac3 = temp*tPDF(pcoll)*xPDF(parton2)
        endif
      if (parton2 .eq. pcoll) then
        print*,'error in lumi_hadron.F: invalid induced process'
	return
      else
	fac3 = fac3*2
      endif	
***]convolute with PDfs, sum over light quarks
	end

