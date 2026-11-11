**************************************************************************  
*     This file is taken from the VBFNLO-2.7.1 program.
*     VBFNLO: https://www.itp.kit.edu/vbfnlo/wiki/doku.php
****************************************************************************
*   This file contains useful functions and subroutines for building a 
*   phasespace generator.
c***************************************************************************
c   LIST OF ALL FUNCTIONS AND SUBROUTINES IN THIS FILE:
c
c  SUBROUTINE ThreeBodyDecay0(R, P_in, M2_in, P_1, P_2, P_3, weight)
c  FUNCTION Resonance(xm2, xmg, rm2min, rm2max, R, weight, qsq)
c  FUNCTION NonResonance(rm2min, rm2max, R, weight, qsq)
c  FUNCTION TwoBodyDecay(R1,R2,P_in,M2_in,M2_out1,M2_out2,P_out1,P_out2,weight)
c  FUNCTION TwoToJetsPlusX(N,R,RN,ECM,PTmin,XQ2,K1,K2,X1,X2,Q,P,weight)
c  SUBROUTINE TwoBodyDecay0(R1, R2, P_in, M2_in, P_out1, P_out2, weight)
c
c***************************************************************************

***************************************************************************
      LOGICAL FUNCTION ThreeBodyDecay(R, RN, P_in, m1sq, m2sq, m3sq, 
     &                  ql1, ql2, ql3, w)
****************************************************************************
*     Calculates the three-body decay P_in -> q1_out + q2_out + q3_out
*
*    INPUT
*       R(4), RN         : five random numbers between 0 and 1
*       P_in             : 4-momentum of mother particle
*       m1sq, m2sq, m3sq : mass squared of outgoing particles
*       w                : phasespace weight
*
*    OUTPUT
*       w                 : new phasespace weight
*       ql1, ql2, ql3     : 4-momenta of daughter particles
*       ThreeBodyDecay    : .true.,  if everything is ok 
*                           .false., if kinematics not ok -> weight = 0
****************************************************************************

      IMPLICIT NONE

      DOUBLE PRECISION RM2, RMG, RM2MIN, RM2MAX, S, PTC
      DOUBLE PRECISION RM2MIN12, RM2MAX12, RM2MIN23, RM2MAX23, M2MIN
      DOUBLE PRECISION RM2MAX1, RM2MAX2
      COMMON /TRIBOSBLIPSI/ RM2(3), RMG(3), RM2MIN(3), RM2MAX(3),
     1                RM2MIN12, RM2MAX12, RM2MIN23, RM2MAX23, M2MIN,
     2                S, PTC, RM2MAX1, RM2MAX2

      DOUBLE PRECISION CLR, XM2, XMG, B, VE, A
      COMMON /BKOPOU/ CLR(4,5,-1:1), XM2(6), XMG(6),
     1                B(6,6,6), VE(4,5), A(4,5)


      DOUBLE PRECISION m1sq, m2sq, m3sq, m12sq, m23sq
      DOUBLE PRECISION P_in(0:4), QL1(0:3), QL2(0:3), QL3(0:3)


      DOUBLE PRECISION M, E1, E2, E3, Abs1, Abs3, C, D
      DOUBLE PRECISION alpha, gamma, cosa, cosb, cosg, sina, sinb, sing

      DOUBLE PRECISION R(4), RN
      DOUBLE PRECISION w
      DOUBLE PRECISION twopi
      PARAMETER (twopi = 2d0*3.141592653589793d0)

c data to be useed in histogram routine
      double precision mZ1sq, mZ2sq, mZ3sq
      Common /histdata1/ m12sq, m23sq, mZ1sq, mZ2sq, mZ3sq, alpha, cosb, gamma

      DOUBLE PRECISION xxwidth, xxm
      PARAMETER ( xxwidth = 10d1, xxm = 45d1 )
      DOUBLE PRECISION x1width, x1m
      PARAMETER ( x1width = 70.0d0, x1m = 90.0d0 )
      DOUBLE PRECISION x2width, x2m
      PARAMETER ( x2width = 70.0d0, x2m = 90.0d0 )

      LOGICAL RESONANCE, NonResonance_KK, NonResonance
      EXTERNAL RESONANCE, NonResonance_KK, NonResonance

      ThreeBodyDecay = .false.

      M = sqrt(p_in(4))

c generate two invariant masses (qi + qj)^2  = mij^2 first **********

c generation of m12sq ---------------------------------
      if (.not. NonResonance(1.0d-5, p_in(4), r(1), W, m12sq)) return

c generation of m23sq next ---------------------------------
      if (.not. NonResonance(1.0d-5, p_in(4), r(2), W, m23sq)) return

c caculate E1, E2, E3, Abs1, Abs2 from that 
      E1 = (p_in(4) + m1sq - m23sq)/2.0d0/M
      E3 = (p_in(4) + m3sq - m12sq)/2.0d0/M
      E2 = M - E1 - E3
      Abs1 = sqrt(E1**2 - m1sq)
      Abs3 = sqrt(E3**2 - m3sq)

c now generate Euler angles (rd(3), rd(4), rn)
      alpha = R(3) * twopi
      W = W*twopi
      cosb = -1.0d0 + R(4)*2.0d0
      W = W*2.0d0
      gamma = RN * twopi
      W = W*twopi

      cosa = cos(alpha)
      cosg = cos(gamma)
      sina = sin(alpha)
      sinb = sqrt(1 - cosb**2)
      sing = sin(gamma)

c calculate now C and D
      D = -(m1sq-m2sq+m3sq+p_in(4)-2.0d0*M*(E1+E3)+2.0d0*E1*E3)/2/Abs1/Abs3
      if (D**2.lt.1) then
         C = sqrt(1 - D**2)
      else
         W = 0.0d0
         return
      endif

c Momenta of V1
      ql1(0) = E1
      ql1(1) = Abs1 * (cosa * cosg - cosb * sina * sing)
      ql1(2) = Abs1 * (-cosb * cosg * sina - cosa * sing)
      ql1(3) = Abs1 * (sina * sinb)

c Momenta of V3
      ql3(0) = E3
      ql3(1) = Abs3 * C * (cosg * sina + cosa * cosb * sing) 
     1      - Abs3 * D * (cosa * cosg - cosb * sina * sing)
      ql3(2) = Abs3 * C * (cosa * cosb * cosg - sina * sing)
     1      - Abs3 * D * (-cosb * cosg * sina - cosa * sing)
      ql3(3) = Abs3 * C * (- cosa * sinb) - Abs3 * D * (sina * sinb)

      call boostn( ql1(0), P_in(0), ql1(0) )
      call boostn( ql3(0), P_in(0), ql3(0) )

c Momenta of V2
      ql2(0) = P_in(0) - ql1(0) - ql3(0)
      ql2(1) = P_in(1) - ql1(1) - ql3(1)
      ql2(2) = P_in(2) - ql1(2) - ql3(2)
      ql2(3) = P_in(3) - ql1(3) - ql3(3)

c Additional weight factors ******************************************
      W = W/twopi**3
      W = W/32.0d0/P_in(4)

      ThreeBodyDecay = .true.

      RETURN
      END

****************************************************************************
      SUBROUTINE ThreeBodyDecay0(R, P_in, M2_in, P_1, P_2, P_3, weight)
****************************************************************************
*     Calculates the three-body decay P_in -> P_1 + P_2 + p_3
*     for MASSLESS decay particles via Euler angles
*
*    INPUT
*       R              : array of 5 random numbers between 0 and 1
*       P_in           : 4-momentum of mother particle
*       M_in           : mass**2 of mother particle
*       weight         : phasespace weight
*
*    OUTPUT
*       weight         : new phasespace weight (multiplied by Jacobi factor)
*       P_1, P_2, P_3  : 4-momenta of daughter particles
****************************************************************************
      IMPLICIT NONE
      DOUBLE PRECISION R(*)
      double precision P_in(0:3), M2_in
      double precision P_1(0:3), P_2(0:3), P_3(0:3), weight
      double precision twopi
      parameter (twopi = 2d0*3.141592653589793d0)
      integer*4 mu
      double precision alpha, sin_alpha, cos_alpha
      double precision gamma, sin_gamma, cos_gamma
      double precision cos_beta, sin_beta, s1, c1
      double precision E1, E2, M
      
      alpha = twopi*R(1)
      cos_beta = 1-2*R(2)
      gamma = twopi*R(3)
      
      if((1-cos_beta**2) .le. 1D-10) then
	weight = 0D0
	return 
      else	
	sin_beta = sqrt(1-cos_beta**2)
      end if
      sin_alpha = sin(alpha)
      cos_alpha = cos(alpha)
      sin_gamma = sin(gamma)
      cos_gamma = cos(gamma)
      
      M = sqrt(M2_in)
      E1 = 0.5d0*R(4)*M
      E2 = 0.5d0*M*(1-R(4)+R(5)*R(4))
      
      P_1(0) = E1
      P_1(1) = E1*(  cos_alpha*cos_gamma - cos_beta*sin_alpha*sin_gamma )
      P_1(2) = E1*( -cos_alpha*sin_gamma - cos_beta*cos_gamma*sin_alpha )
      P_1(3) = E1*sin_alpha*sin_beta

      c1 = -(2d0*E1*E2 - 2d0*(E1 + E2)*M + M2_in)/(2d0*E1*E2)
      if ((1D0-c1**2) .le. 1D-10) then
         weight = 0D0
         return
      end if
      s1 = sqrt(1-c1**2)

      P_2(0) = E2
      P_2(1) = E2*( s1*(cos_gamma*sin_alpha + cos_alpha*cos_beta*sin_gamma) - 
     &              c1*(cos_gamma*cos_alpha - sin_alpha*cos_beta*sin_gamma) )
      P_2(2) = E2*( s1*( cos_alpha*cos_beta*cos_gamma - sin_alpha*sin_gamma) - 
     &              c1*(-sin_alpha*cos_beta*cos_gamma - cos_alpha*sin_gamma) )
      P_2(3) = E2*sin_beta*(-s1*cos_alpha - c1*sin_alpha)


      CALL BOOSTN(P_1, P_in, P_1)  
      CALL BOOSTN(P_2, P_in, P_2)

      do mu = 0,3
         P_3(mu) = P_in(mu)-P_1(mu)-P_2(mu)
      enddo

      weight = weight*(1d0/twopi**3)*(M2_in*R(4))/(16d0)
      
      RETURN
      END
****************************************************************************
      SUBROUTINE TwoToOne(R,ECM,Qsqr,K1,K2,X1,X2,Q,weight)
****************************************************************************
*     Calculates the phasespace for 2 -> 1 massive
*
*    INPUT
*       R            : random number
*       ECM          : center-of-mass energy
*       Qsqr         : Q**2 of the outgoing particle
*       weight       : phasespace weight
*
*    OUTPUT
*       K1,K2  : 4-vectors for 2 incoming partons
*       X1,X2  : corresponding values of Feynman x
*       Q      : 4-vector of outgoing particle
*       weight : new phasespace weight (multiplied by Jacobi factor)
****************************************************************************

      IMPLICIT NONE

      DOUBLE PRECISION R, ECM, Qsqr
      DOUBLE PRECISION K1(0:3), K2(0:3), Q(0:3)
      DOUBLE PRECISION X1, X2, weight
      DOUBLE PRECISION TPI, PI
      PARAMETER (PI=3.141592653589793d0)
      PARAMETER (TPI=2D0*PI )

c local variables
      INTEGER*4 mu
      DOUBLE PRECISION s, tau, y, y_border

      S = ecm**2

c Generate y, tau -> x1, x2, shat from R first **************
      tau = Qsqr/s  
      weight = weight / s 
      y_border = log(1/tau)
      y = y_border*( R - 0.5d0 )
      weight = weight * y_border !(prop density)^-1 of y

      x1 = sqrt(tau) * exp(y)
      x2 = sqrt(tau) * exp(-y)

c generate incoming momenta in lab frame out of x1 x2
      k1(0) = 0.5d0 * x1 * ecm
      k1(1) = 0d0
      k1(2) = 0d0
      k1(3) = 0.5d0 * x1 * ecm

      k2(0) = 0.5d0 * x2 * ecm
      k2(1) = 0d0
      k2(2) = 0d0
      k2(3) = - 0.5d0 * x2 * ecm

      Q(1) = 0d0
      Q(2) = 0d0
      do mu = 0, 3, 3
         Q(mu) = k1(mu)+k2(mu)
      enddo

      weight = weight * 0.5d0/Qsqr        ! Flux
      weight = weight * TPI                  ! (2*pi)^(4-3) from delta- and ps-int
      weight = weight * 3.89379304d11        ! conversion into fb

      return
      end


****************************************************************************
      SUBROUTINE TwoToTwo(R,RN,ECM,Qsqr1,Qsqr2,K1,K2,X1,X2,Q1,Q2,weight)
****************************************************************************
*     Calculates the phasespace for 2 -> 2 massive
*
*    INPUT
*       R            : array with 3 random numbers between 0 and 1
*       RN           : one additional random number 
*       ECM          : center-of-mass energy square
*       Qsqr1, Qsqr2 : Q**2 of the two outgoing particles
*       weight       : phasespace weight
*
*    OUTPUT
*       K1,K2  : 4-vectors for 2 incoming partons
*       X1,X2  : corresponding values of Feynman x
*       Q1,Q2  : 4-vectors of outgoing particles
*       weight : new phasespace weight (multiplied by Jacobi factor)
****************************************************************************

      IMPLICIT NONE

      DOUBLE PRECISION R(*), RN, ECM, Qsqr1, Qsqr2
      DOUBLE PRECISION K1(0:3), K2(0:3), Q1(0:3), Q2(0:3)
      DOUBLE PRECISION X1, X2, weight
      DOUBLE PRECISION TPI, PI
      PARAMETER (PI=3.141592653589793d0)
      PARAMETER (TPI=2D0*PI )

c local variables
      INTEGER*4 mu
      DOUBLE PRECISION k(0:3)
      DOUBLE PRECISION s, ptmax, y2_border, temp, pt, ET, y_border, y, y2
      DOUBLE PRECISION phi, sinhy2, coshy2, abs_p, pL, tau
      DOUBLE PRECISION A, ptc, ptmin

c test variables
      COMMON /ptdistr/ pt
      
      !print*,"R,RN:",R(1),R(2),R(3),RN
      !print*,"ECM=",ECM
      !print*,"Qsqr1,Qsqr2=",Qsqr1,Qsqr2
      !print*,"weight =",weight

      S = ecm**2

c generate y, p_T, y*, phi from plain distributed R(1),R(2),R(3), rn

      phi = rn * tpi - tpi/2d0

c ptc is the 'changeover point' of the phasespace generation - 
c below we have a constant value, and above it's more complicated
      ptc = 10.0d0
      pTmax = ecm/2.0d0
c      ptmax = 0.5d0*sqrt((s-(Qsqr1+Qsqr2))**2 - 4d0*Qsqr1*Qsqr2)/ecm
      A = pTmax/ptc*1.0d0/(2*pTmax-ptc)

      if (r(1) .le. A * ptc) then
         pT = r(1)*1.0d0/A
         Weight = Weight/A
      else
         pT = pTmax*ptc*(1-A*ptc)/(pTmax-A*ptc**2-r(1)
     &               *(pTmax-ptc))
         Weight = Weight/A/ptc**2*pT**2
      endif
c pT-distribution not continuous------------------------------
c      if (r(1).le.0.2) then
c         pTmax = ptc
c         pTmin = 0.0d0       
c         pT = r(1)*(pTmax-ptmin)/0.2d0 + ptmin
c         Weight= Weight/0.2d0
c      else
c         pTmin = 2.0d0/ecm
c         pTmax = 1/ptc
c         pT = pTmin + (r(1)-0.2)/0.8d0*(pTmax-pTmin)
c         pT = 1/pT
c         Weight = Weight/0.8d0*pT**2
c      endif
c      Weight = Weight*(ptmax-ptmin)


      ET = sqrt(pT**2+Qsqr1)
      temp = (s+Qsqr1-Qsqr2)/(2*ecm*ET)

c calculate y*
      if ((temp**2 - 1D0) .gt. 0D0) then
         y2_border = log( temp+sqrt(temp**2-1) )
         y2 = y2_border*2d0*( R(2) - 0.5d0 ) ! y*
         weight = weight * y2_border * 2d0 !(prop density)^-1 of y2

         sinhy2 = sinh(y2)
         coshy2 = sqrt(sinhy2**2+1)
         abs_p = sqrt(pT**2*coshy2**2+Qsqr1*sinhy2**2)
         pL = sinhy2*ET
         tau=((coshy2*ET+sqrt(coshy2**2*ET**2-Qsqr1+Qsqr2))**2)/s
         weight = weight * pT/(sqrt(tau*s))
         weight = weight/(sqrt(coshy2**2*ET**2-Qsqr1+Qsqr2))

         y_border = log(1/tau)
         y = y_border*( R(3) - 0.5d0 ) !y
         weight = weight * y_border !(prop density)^-1 of y


c------------------------------------------------------------------------
c generate incoming momenta in lab frame + x1 x2

         x1 = sqrt(tau) * exp(y)
         x2 = sqrt(tau) * exp(-y)

         k1(0) = 0.5d0 * x1 * ecm
         k1(1) = 0d0
         k1(2) = 0d0
         k1(3) = 0.5d0 * x1 * ecm

         k2(0) = 0.5d0 * x2 * ecm
         k2(1) = 0d0
         k2(2) = 0d0
         k2(3) = - 0.5d0 * x2 * ecm

         do mu = 0, 3
            K(mu) = k1(mu)+k2(mu)
         enddo
c Momenta of V1

         q1(0) = sqrt( pT**2 + pL**2 + Qsqr1 )
         q1(1) = pT * cos(phi)
         q1(2) = pT * sin(phi)
         q1(3) = pL


c Momenta of V2

         q2(0) = sqrt( pT**2 + pL**2 + Qsqr2 )
         q2(1) = - pT * cos(phi)
         q2(2) = - pT * sin(phi)
         q2(3) = - pL

c boost in lab frame

         call boostn( q1(0), K, q1(0) )
         call boostn( q2(0), K, q2(0) )

         weight = weight * 1d0/tpi/4d0/s
         weight = weight * 3.89379304d11 !conversion into fb

      else
         weight = 0D0
      end if


      RETURN
      END
****************************************************************************
      LOGICAL FUNCTION Resonance(xm2, xmg, rm2min, rm2max, R, weight, qsq)
****************************************************************************
*     Calculates the Q-square for a resonance
*
*    INPUT
*       xm2       : mass**2 of the resonance
*       xmg       : mass*width of the resonance
*       rm2min    : lowerBound**2 of the resonance
*       rm2max    : upperBound**2 of the resonance
*       R         : a random number between 0 and 1
*       weight    : phasespace weight
*
*    OUTPUT
*       weight    : new phasespace weight (multiplied by Jacobi factor)
*       qsq       : the Q**2 
*       Resonance : .true.,  if everything is ok 
*                   .false., if kinematics not ok -> weight = 0
****************************************************************************
      IMPLICIT NONE

      DOUBLE PRECISION xm2, xmg, rm2min, rm2max, R, weight, qsq, tanx
      DOUBLE PRECISION twopi
      PARAMETER (twopi = 2d0*3.141592653589793d0)
      DOUBLE PRECISION xmin, xmax, deltax, x 

#ifdef DEBUG
         print*, "Begin Resonance"
         print*, "R =", R
#endif

      IF (rm2max.LT.rm2min) THEN
         weight = 0d0 
         Resonance = .false.
         RETURN
      ENDIF
      
      xmin = ATAN( (rm2min - xm2)/xmg )
      xmax = ATAN( (rm2max - xm2)/xmg )
      deltax = xmax - xmin
      x = xmin + deltax * R
      tanx = tan(x)
      qsq = xm2 + tanx * xmg

      
      IF (qsq.LE.0 .or. deltax.le.0) THEN
         weight = 0d0 
         Resonance = .false.
         RETURN 
      ENDIF

      weight = weight * deltax / twopi
      weight = weight * xmg * (1 + tanx**2 )

      Resonance = .true.

      RETURN 
      END
      

****************************************************************************
      LOGICAL FUNCTION NonResonance(rm2min, rm2max, R, weight, qsq)
****************************************************************************
*     Calculates the Q-square for a non-resonance (a flat distribution)
*
*    INPUT
*       rm2min       : lowerBound**2 of the resonance
*       rm2max       : upperBound**2 of the resonance
*       R            : a random number between 0 and 1
*       weight       : phasespace weight
*
*    OUTPUT
*       weight       : new phasespace weight (multiplied by Jacobi factor)
*       qsq          : the Q**2 
*       NonResonance : .true.,  if everything is ok 
*                      .false., if kinematics not ok -> weight = 0
****************************************************************************
      IMPLICIT NONE

      DOUBLE PRECISION rm2min, rm2max, R, weight, qsq
      DOUBLE PRECISION pi
      PARAMETER (pi = 3.141592653589793d0)
      DOUBLE PRECISION minMass, maxMass, deltaMass, mass

      IF (rm2max.LT.rm2min) THEN
         weight = 0d0 
         NonResonance = .false.
         RETURN 
      ENDIF
      
      minMass = sqrt(rm2min)
      maxMass = sqrt(rm2max)
      deltaMass = maxMass - minMass
      mass = minMass + deltaMass * R
      qsq = mass**2

      weight = weight * deltaMass * mass / pi

      NonResonance = .true.
      RETURN 
      END

****************************************************************************
      LOGICAL FUNCTION NonResonance_Log(rm2min, rm2max, R, weight, qsq)
****************************************************************************
*     Calculates the Q-square for a non-resonance (a flat distribution)
*
*    INPUT
*       rm2min       : lowerBound**2 of the resonance
*       rm2max       : upperBound**2 of the resonance
*       R            : a random number between 0 and 1
*       weight       : phasespace weight
*
*    OUTPUT
*       weight       : new phasespace weight (multiplied by Jacobi factor)
*       qsq          : the Q**2 
*       NonResonance : .true.,  if everything is ok 
*                      .false., if kinematics not ok -> weight = 0
****************************************************************************
      IMPLICIT NONE

      DOUBLE PRECISION rm2min, rm2max, R, weight, qsq
      DOUBLE PRECISION pi
      PARAMETER (pi = 3.141592653589793d0)
      DOUBLE PRECISION minMass, maxMass, delta, mass
      DOUBLE PRECISION minLog, maxLog, t

      IF (rm2max.LT.rm2min) THEN
         weight = 0d0 
         NonResonance_Log = .false.
         RETURN 
      ENDIF
      
      minMass = sqrt(rm2min)
      maxMass = sqrt(rm2max)
      minLog = Log(minMass)
      maxLog = Log(maxMass)

      delta = maxLog - minLog
      t = minLog + delta*R
      mass = exp(t) 
      qsq = mass**2

      weight = weight*delta*qsq/pi

      NonResonance_Log = .true.
      RETURN 
      END
****************************************************************************
      LOGICAL FUNCTION NonResonance_Pow(a, rm2min, rm2max, R, weight, qsq)
****************************************************************************
*     Calculates the Q-square for a non-resonance (a flat distribution)
*
*    INPUT
*       a: transformation parameter (t = q^a) [suggestion: a = 0.5]
*       rm2min       : lowerBound**2 of the resonance
*       rm2max       : upperBound**2 of the resonance
*       R            : a random number between 0 and 1
*       weight       : phasespace weight
*
*    OUTPUT
*       weight       : new phasespace weight (multiplied by Jacobi factor)
*       qsq          : the Q**2 
*       NonResonance : .true.,  if everything is ok 
*                      .false., if kinematics not ok -> weight = 0
****************************************************************************
      IMPLICIT NONE

      DOUBLE PRECISION a, rm2min, rm2max, R, weight, qsq
      DOUBLE PRECISION pi
      PARAMETER (pi = 3.141592653589793d0)
      DOUBLE PRECISION minMass, maxMass, delta, mass
      DOUBLE PRECISION minT, maxT, t, b

      IF (rm2max.LT.rm2min) THEN
         weight = 0d0 
         NonResonance_Pow = .false.
         RETURN 
      ENDIF

      IF (a.eq.0) THEN
         print*,"NonResonance_Pow(): wrong input: a = 0" 
         weight = 0d0 
         NonResonance_Pow = .false.
         RETURN 
      ENDIF
      
      minMass = sqrt(rm2min)
      maxMass = sqrt(rm2max)
      minT = minMass**a
      maxT = maxMass**a

      delta = maxT - minT
      t = minT + delta*R
      b = 1d0/a
      mass = t**b
      qsq = mass**2

      weight = weight*delta*b*qsq/(pi*t)

      NonResonance_Pow = .true.
      RETURN 
      END

****************************************************************************
      LOGICAL FUNCTION NonResonance_kk(xm2, xmg, expo, rm2min, rm2max,
     &                                                R, weight, qsq)
****************************************************************************
*    Calculates q² for atan((q²-xm2)/xmg)**expo / ((q²-xm2)²+xmg²) distribution
*    !!!!!!!!!!!   xm2 needs to be lower than rm2min  !!!!!!!!!!!!!!!!!!!!!
*    INPUT
*       rm2min       : lowerBound**2 of the distribution
*       rm2max       : upperBound**2 of the distribution
*       R            : a random number between 0 and 1
*       weight       : phasespace weight
*       xm2, xmg, expo: parameters of distribution
*    OUTPUT
*       weight       : new phasespace weight (multiplied by Jacobi factor)
*       qsq          : the Q**2
*       NonResonance : .true.,  if everything is ok     
*                      .false., if kinematics not ok -> weight = 0
****************************************************************************

      IMPLICIT NONE

      double precision rm2min, rm2max, R, weight, qsq, xmg, xm2, expo
      double precision twopi
      parameter (twopi = 2d0*3.141592653589793d0)
      double precision xmax, x, tanx, deltax, xmin, n

      IF (rm2max.LT.rm2min) THEN
         weight = 0d0
         NonResonance_kk = .false.
         RETURN
      ENDIF

      IF (xm2 .ge. rm2min) THEN
         weight = 0d0
         NonResonance_kk = .false.
         RETURN
      ENDIF

      n=1d0/(expo+1d0)

      xmin = ATAN( (rm2min-xm2)/xmg )**(1/n)
      xmax = ATAN( (rm2max-xm2)/xmg )**(1/n)
      deltax= xmax -xmin
      x =xmin + deltax * R

      tanx = tan(x**n)

      qsq = xm2 + tanx * xmg

      weight = weight * deltax / twopi
      weight = weight * xmg * (1 + tanx**2 )
      weight = weight * n * x**(n-1)

      NonResonance_kk = .true.

      RETURN
      END

****************************************************************************
      LOGICAL FUNCTION TwoBodyDecay(R1, R2, P_in, M2_in, M2_out1, M2_out2, 
     &                              P_out1, P_out2, weight)
****************************************************************************
*     Calculates the two-body decay P_in -> P_out1 + P_out2
*
*    INPUT
*       R1, R2         : two random numbers between 0 and 1
*       P_in           : 4-momentum of mother particle
*       M_in           : mass**2 of mother particle
*       M_out1, M_out2 : mass**2 of daughter particles
*       weight         : phasespace weight
*
*    OUTPUT
*       weight         : new phasespace weight (multiplied by Jacobi factor)
*       P_out1, P_out2 : 4-momenta of daughter particles
*       TwoBodyDecay   : .true.,  if everything is ok 
*                        .false., if kinematics not ok -> weight = 0
****************************************************************************
      IMPLICIT NONE

      DOUBLE PRECISION R1, R2, P_in(0:3), M2_in, M2_out1, M2_out2
      DOUBLE PRECISION P_out1(0:3), P_out2(0:3), weight
      DOUBLE PRECISION twopi
      PARAMETER (twopi = 2d0*3.141592653589793d0)
      INTEGER*4 mu
      double precision Q(0:3), CTH, STH, PHID, M_in, M_out1, M_out2, E1, P1

#ifdef DEBUG
         print*, "Begin TwoBodyDecay"
         print*, "R1, R2 =", R1, R2
#endif

      M_in = sqrt(M2_in)
      M_out1 = sqrt(M2_out1)
      M_out2 = sqrt(M2_out2)

      IF (M_in .LT. (M_out1 + M_out2) ) THEN
         weight = 0d0 
         TwoBodyDecay = .false.
         RETURN 
      ENDIF

c generate outgoing momenta in CMS and boost back
      CTH = -1d0 + 2d0 * R1
      if ((1D0 - CTH**2) .le. 1D-10) then
         weight = 0D0
         return
      end if
      STH = sqrt(1d0 - CTH**2 )
      PHID = twopi * R2 - twopi/2.0d0
      
      E1 = (M2_in - M2_out2 + M2_out1) / (2d0*M_in)
      if ((E1**2 - M2_out1) .le. 1D-10) then
         weight = 0d0
         TwoBodyDecay = .false.
         return
      end if
      P1 = sqrt(E1**2 - M2_out1)
      Q(0) = E1
      Q(1) = P1 * STH*COS(PHID)
      Q(2) = P1 * STH*SIN(PHID)
      Q(3) = P1 * CTH

      call boostn(Q,P_in,P_out1)

      do mu = 0,3
         P_out2(mu) = P_in(mu) - P_out1(mu)
      enddo

      weight = weight * P1/(M_in * 2d0 * twopi) 
      
      TwoBodyDecay = .true.
      RETURN
      END

****************************************************************************
      SUBROUTINE TwoBodyDecay0(R1, R2, P_in, M2_in, P_out1, P_out2, weight)
****************************************************************************
*     Calculates the two-body decay P_in -> P_out1 + P_out2 
*     for MASSLESS decay particles 
*
*    INPUT
*       R1, R2         : two random numbers between 0 and 1
*       P_in           : 4-momentum of mother particle
*       M_in           : mass**2 of mother particle
*       weight         : phasespace weight
*
*    OUTPUT
*       weight         : new phasespace weight (multiplied by Jacobi factor)
*       P_out1, P_out2 : 4-momenta of daughter particles
****************************************************************************
      IMPLICIT NONE
      double precision R1, R2, P_in(0:3), M2_in
      double precision P_out1(0:3), P_out2(0:3), weight
      double precision twopi
      parameter (twopi = 2d0*3.141592653589793d0)
      integer*4 mu
      double precision Q(0:3), CTH, STH, PHID, M_in, E1

#ifdef DEBUG
         print*, "Begin TwoBodyDecay0"
         print*, "R1, R2 =", R1, R2
#endif

      M_in = sqrt(M2_in)

c generate outgoing momenta in CMS and boost back
      CTH = -1d0 + 2d0 * R1
      if ((1D0 - CTH**2) .le. 1D-10) then
         weight = 0D0
         return
      end if
      STH = sqrt(1d0 - CTH**2 )
      PHID = twopi * R2 - twopi/2d0
      
      E1 = 0.5d0 * M_in
      Q(0) = E1
      Q(1) = E1 * STH*COS(PHID)
      Q(2) = E1 * STH*SIN(PHID)
      Q(3) = E1 * CTH

      call boostn(Q,P_in,P_out1)
      do mu = 0,3
         P_out2(mu) = P_in(mu) - P_out1(mu)
      enddo

      weight = weight / (4d0 * twopi) 
      
      RETURN
      END

****************************************************************************
      LOGICAL FUNCTION TwoToJetsPlusX(N,R,RN,ECM,PTmin,XQ2,K1,K2,X1,X2,Q,P,weight)
****************************************************************************
*     Calculates the phasespace for the 2 -> N massless jets + X
*     X is a pseudo-particle with arbitrary mass>0
*
*    INPUT
*       R      : array with 3*N random numbers between 0 and 1
*       RN     : one additional random number 
*       ECM    : center-of-mass energy
*       PTmin  : minimum PT of jets
*       XQ2    : Q**2 of X-particle
*       weight : phasespace weight
*
*    OUTPUT
*       K1,K2  : 4-vectors for 2 incoming partons
*       X1,X2  : corresponding values of Feynman x
*       Q      : 4-vector of X-particle
*       P      : 4-vectors for the N massless jets
*       weight : new phasespace weight (multiplied by Jacobi factor)
****************************************************************************
      IMPLICIT NONE
      INTEGER*4 N, N1
      PARAMETER (N1 = 5)
      DOUBLE PRECISION TPI,PI
      PARAMETER (PI=3.141592653589793d0)
      PARAMETER (TPI=2D0*PI )
      DOUBLE PRECISION R(*),RN,ECM,PTmin,XQ2,
     1       K1(0:3),K2(0:3),Q(0:3),P(0:3,N1),
     2       X1,X2,weight
      
      INTEGER*4 I,J
      DOUBLE PRECISION S,FAC,DELY,XMIN,XMAX,X,DELX,YCM,YMIN,YMAX,YBAR
      DOUBLE PRECISION PCM(0:3), PNJ(0:4), Y(N1), PT(N1), PHI(N1)
      DOUBLE PRECISION SINHY, COSHY, M2, M2min, QV, XMJ, PJT2, SHA, ACM
      DOUBLE PRECISION SUMPST, Q0ST, RSHAT
      PARAMETER (m2min = 0.1d0**2)

      double precision dotrr
      external dotrr

#ifdef DEBUG
         print*, "Begin TwoToJetsPlusX"
         print*, "N =", N
	 do I = 1, 3*N
         print*, "R(i) =", I, R(I)
	 enddo
         print*, "RN =", RN
         print*, "ECM =", ECM
         print*, "PTmin,XQ2 =", PTmin, XQ2
#endif

      S = ECM**2
      FAC = 1d0/(16d0*PI**3)
      weight = weight*TPI

C  Reset local variables and incoming momenta      
      DO I=0,3
         PCM(I) = 0d0
         PNJ(I) = 0d0
         K1(I) = 0d0
         K2(I) = 0d0
      ENDDO
      PNJ(4) = 0d0

C  Get 4-momenta of the N massless jets
      DO I = 1,N
         weight = weight*FAC
         J = 3*(I-1)
         
         if (R(J+2).gt.0.2) then
            XMIN = 2D0/ECM
            XMAX = 1/PTmin
            DELX = XMAX-XMIN
            X = XMIN + (R(J+2)-0.2)/0.8d0*DELX
            PT(I) = 1/X
            weight = weight/0.8d0*DELX*PT(I)**3
         else
            xmin = 0
            xmax = PTmin
            DELX = XMAX-XMIN
            PT(i) = XMIN + R(J+2)/0.2d0*DELX
            weight = weight/0.2d0*DELX*PT(I)
         endif
         YMAX = ECM/(2d0*PT(I))
         YMAX = LOG(YMAX+SQRT(YMAX**2-1d0))

         DELY = 2*YMAX

         IF (I.LE.3) THEN
            Y(I) = -YMAX+R(J+1)*DELY
         ELSE
C  fix random number to 0 or 1 for delta-y = 0 with previous jet: R=0 or 1 
C  corresponds to lowest dijet invariant mass
            IF (Y(I-1).LT.0D0) THEN
               IF (Y(I-1).LT.-YMAX) THEN
                  Y(I) = -YMAX + R(J+1)*DELY
               ELSEIF (Y(I-1).LT.0d0) THEN
                  Y(I) = Y(I-1) + R(J+1)*DELY
                  IF (Y(I).GT.YMAX) Y(I) = -2*YMAX + Y(I)
               ELSE
                  Y(I) =  - R(J+1)*DELY
                  IF (Y(I).LT.-YMAX) Y(I) = 2*YMAX + Y(I)
               ENDIF
            ELSE
               IF (Y(I-1).GT.YMAX) THEN
                  Y(I) = YMAX - R(J+1)*DELY
               ELSEIF (Y(I-1).GT.0D0) THEN
                  Y(I) = Y(I-1) - R(J+1)*DELY
                  IF (Y(I).LT.-YMAX) Y(I) = 2*YMAX + Y(I)
               ELSE
                  Y(I) = R(J+1)*DELY
                  IF (Y(I).GT.YMAX) Y(I) = -2*YMAX + Y(I)
               ENDIF
            ENDIF
         ENDIF

         SINHY = SINH(Y(I))
         COSHY = SQRT(1+SINHY**2)
         weight = weight*DELY

         P(0,I) = PT(I)*COSHY

         IF (I.EQ.1) THEN
            PHI(I) = PI*(2*RN-1)
         ELSE
            PHI(I) = PHI(I-1) + TPI*R(J+3)
            if (phi(i).gt.pi) phi(i) = phi(i)-tpi
         ENDIF

         weight = weight*TPI

         P(1,I) = PT(I)*COS(PHI(I))
         P(2,I) = PT(I)*SIN(PHI(I))
         P(3,I) = PT(I)*SINHY

         DO J = 0,3
            PNJ(J) = PNJ(J)+P(J,I)
         ENDDO
      ENDDO

C  determine N-jet invariant mass
      IF (N.GT.1) THEN
         PNJ(4) = PNJ(0)**2-PNJ(1)**2-PNJ(2)**2-PNJ(3)**2
      ELSE
         PNJ(4) = 1D-18
      ENDIF
      IF (PNJ(4).GE.S) THEN
         weight=0d0
         TwoToJetsPlusX = .false.
         RETURN
      ENDIF

c  check invariant masses of all parton pairs: eliminate events which 
c  may cause numerical problems 
      do I = 1,N-1
         do J = I+1,N
            m2 = 2*dotrr(p(0,I),p(0,J))
            if (m2.lt.m2min) then
               weight = 0d0
               TwoToJetsPlusX = .false.
               RETURN
            endif
         enddo
      enddo
      
C  Next get the center of mass rapidity YCM
      QV = SQRT(XQ2)
      XMJ = SQRT(PNJ(4))
      YBAR = 0.5D0*LOG( (PNJ(0)+PNJ(3))/(PNJ(0)-PNJ(3)) )
      PJT2 = PNJ(1)**2+PNJ(2)**2
      SHA = (S-(QV+XMJ)**2-4*PJT2) / (4*(PNJ(4)+PJT2)) 
      IF (SHA.LE.0D0) THEN
        weight=0
        TwoToJetsPlusX = .false.
        RETURN
      ENDIF
      SHA = SQRT(SHA)
      ACM = LOG( SHA+SQRT(SHA**2+1) ) 
      YMIN = YBAR - ACM
      YMAX = YBAR + ACM

      DELY = YMAX-YMIN
      YCM = YMIN + R(3)*DELY
      SINHY = SINH(YCM)
      COSHY = SQRT(1+SINHY**2)
      weight = weight*DELY

C  now get Q and SHAT etc.
      SUMPST = PNJ(1)**2+PNJ(2)**2+(PNJ(3)*COSHY-PNJ(0)*SINHY)**2
      Q0ST = SQRT(XQ2+SUMPST)
      RSHAT = Q0ST + SQRT(PNJ(4)+SUMPST)
      PCM(0) = RSHAT*COSHY
      PCM(3) = RSHAT*SINHY

      DO J = 0,3
        Q(J) = PCM(J)-PNJ(J)
      ENDDO

      X1 = (PCM(0)+PCM(3))/ECM
      X2 = (PCM(0)-PCM(3))/ECM

      IF (X1.GE.1.0D0 .OR. X2.GE.1.0D0) THEN
         weight=0d0
         TwoToJetsPlusX = .false.
         RETURN
      ENDIF
      K1(0) =  X1*ECM/2
      K1(3) =  K1(0)
      K2(0) =  X2*ECM/2
      K2(3) = -K2(0)

c  check initial-final momentum transfer for all partons: eliminate events 
c  which may cause numerical problems 
!       do I = 1,N         
!          m2 = 2*min( p(0,I)*k1(0)-p(3,I)*k1(3),
!      1               p(0,I)*k2(0)-p(3,I)*k2(3) )
!          if (m2.lt.m2min) then
!             weight = 0d0
!             TwoToJetsPlusX = .false.
!             RETURN
!          endif
!       enddo

C  insert Jacobian for dtau * delta(PCM(0)-Ql(0,0)-PNJ(0))/(2 Ql(0,0))
      weight = weight * RSHAT/(S*Q0ST)
      
C  finally put in flux factor for the parton collisions and conversion to fb
      weight = weight / (2*RSHAT**2) * 3.89379304D+11

      TwoToJetsPlusX = .true.
      RETURN
      END

***********************************************************************
C********** BOOSTN ****************************************************
C
      SUBROUTINE BOOSTN(P,R,Q)
C
C
C     The four vector P is assumed to be given in the rest frame of R,
C     which must be a timelike vector.
C     output Q is the vector P boosted to the frame in which R is given.
C                                              Compare Jackson, p.517
C                                              D. Zeppenfeld (28.6.1985)
C     New version with energy stored in zeroth component of four vector
C     arrays. Checked on May 27, 1988.
C
      REAL*8 P(0:3),R(0:3),Q(0:3)
      REAL*8 BETA(3), X, Y, GAMMA
      INTEGER I
      X = 0D0
      Y = 0D0
      DO I = 1,3
         BETA(I) = R(I)/R(0)
         X = X + BETA(I)**2
         Y = Y + BETA(I)*P(I)
      ENDDO
      IF (X.GT.1D-16.AND.X.LT.(1D0-1D-12)) THEN
         GAMMA = 1D0/DSQRT(1D0-X)
         DO I = 1,3
            Q(I) = P(I)+BETA(I)*(Y*(GAMMA-1D0)/X + GAMMA*P(0))
         ENDDO
         Q(0) = GAMMA*(P(0) + Y)
      ELSE
         DO I = 0,3
            Q(I) = P(I)
         ENDDO
         IF(X.GE.(1D0-1D-12)) 
     *      WRITE(6,1000) R,R(0)**2-R(1)**2-R(2)**2-R(3)**2
      ENDIF
 1000 FORMAT (" The reference vector ",4G12.3," is not timelike."/
     1        " R**2 = ",G12.3)
      END

c**********************************************************************
      function dotrr(p1,p2)
c***************************************************************************
c
c     dotrr(p1,p2) = p1.p2
c
c***************************************************************************
      implicit none

      double precision dotrr,p1(0:3),p2(0:3)

      dotrr = p1(0)*p2(0) - p1(1)*p2(1) - p1(2)*p2(2) - p1(3)*p2(3)

      end
