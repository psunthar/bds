!!! Time-stamp: <gsipceq.f90 14:12, 22 Sep 2005 by P Sunthar>

!________________________________________________
!   Bead Spring Configuration Space Utilities
!________________________________________________

!!! $Log: gsipceq.f90,v $
!!! Revision 1.4  2005/05/24 01:06:20  pxs565
!!! With variable NTOPB
!!!
!!! Revision 1.3  2005/05/02 11:24:40  pxs565
!!! Tuned to work with intel/lapack
!!!
!!! Revision 1.2  2005/03/16 00:49:09  pxs565
!!! Preaveraged RPY working
!!!
!!! Revision 1.1  2005/02/25 10:21:22  pxs565
!!! Initial revision
!!!


Module csputls 
!!! This module contains subroutines used to manipulate the 
!!! configuration variable, R and its derivatives
  Use bspglocons

Contains

  Subroutine b2bvector_sym_up(N,R,b2b)
    Implicit None
    ! Calling arguments
    Integer, Intent(in)  :: N
    Real (FPREC) , Intent(in)  :: R(:,:)
    Real (FPREC) , Intent(out) :: b2b(:,:,:)

!!!_________________________________________________________|
!!! This routine returns the vector one bead to another bead
!!! as a symmetric matrix but only the upper elements filled
!!!_________________________________________________________|

    Integer mu,nu

    ! according to /usr/bin/prof, this takes max time. so commenting out
    ! since only upper diagonal is anyway reqd
    ! b2b = 0

    ! note only upper half (nu > mu) of the matrix is filled
    ! and the convention is R_12 = R_2 - R_1
    ! where R_12 is the vector from 1 to 2
    Do nu = 2,N
       Do mu = 1,nu-1
          b2b(:,mu,nu) = R(:,nu) - R(:,mu)
       End Do
    End Do
  End Subroutine b2bvector_sym_up

  Subroutine modr_sym_up(N,b2bvec,deltaR) 
    Implicit None
    ! calling arguments
    Integer, Intent (in)  :: N
    Real (FPREC) , Intent (in)  :: b2bvec(:,:,:)
    Real (FPREC) , Intent (out) :: deltaR(:,:) 
!!!_________________________________________________________|
!!! This subroutine returns the magnitude of each of the bead 
!!! to bead vector
!!!_________________________________________________________|

    Integer mu,nu,i
    Real (FPREC)  r12(Ndim),modr


    deltaR = 0

    ! note that we cant use a forall since dot_product is a
    ! transformational function
    Do nu = 2,N
       Do mu = 1,nu-1
          r12 = b2bvec(:,mu,nu)
          modr = 0.0
          Do i=1,Ndim
             modr = modr + r12(i) * r12(i)
          End Do
          If (modr < 1.0e-12) modr = 1.0e-12   !so that initial dr is != 0
          deltaR(mu,nu) = Sqrt(modr)
       End Do
    End Do
  End Subroutine modr_sym_up


  Subroutine tensor_prod_of_vec(alpha,vec,tensor)
    Implicit None
    Real (FPREC) , Intent (in) :: alpha
    Real (FPREC) , Intent (in) :: vec(:)
    Real (FPREC) , Intent (out) ::  tensor(:,:)
    Integer i,j
    Real (FPREC)  temp

    Do j = 1,Ndim
       Do i = 1,j
          temp = alpha * vec(i) * vec(j)
          tensor(i,j) = temp
          tensor(j,i) = temp
       End Do
    End Do
  End Subroutine tensor_prod_of_vec
End Module csputls



Module BSpModel
!!! This module contains the constants and functions which are specific
!!! to the model of spring and EV.
  Use bspglocons
  Implicit None

  ! These are paramters that will be used  by the following procedures
  Real (FPREC)  sqrtb

  ! paramters for EV
  Real (FPREC)  zsbyds5,pt5bydssq

  ! preaveraged HI
  !Real (FPREC), save, dimension(Ndim,NBeads,Ndim,NBeads) ::  Bpmat,Dpmat
  Real (FPREC), save, Allocatable, Dimension(:,:,:,:) ::  Bpmat,Dpmat


Contains 

  Subroutine force_sans_hookean(sptype,sr,ff)
    Integer, Intent (in) :: sptype 
    Real (FPREC) , Intent (in) :: sr
    Real (DOBL) , Intent (out) :: ff

    Real (DOBL) r

    r = sr 
    ! quick way to do double prec calculation.  otherwise, the complete 
    ! code needs to be Dprec to incorporate accuracy at every level

    ! compute the force factor, ff
    Select Case (sptype)
    Case (HOOK) ! Hookean
       ff = 1.0

    Case (FENE) ! Warner Spring
       ff = 1./(1 - r*r)

    Case (ILC)  ! Inverse Langevin, Pade approx 
       ff = (3-r*r)/(1 - r*r)/3.

    Case (WLC) ! Worm-like Chain, Marko-Siggia interpolation
       ff = 1./(6*r) * ( 4*r+ 1./(1-r)/(1-r) - 1)

    End Select
  End Subroutine force_sans_hookean

  Subroutine Spring_force(sptype,N,b2bvec,dr,Fs)
    Integer, Intent (in) :: sptype
    Integer, Intent (in)  :: N
    Real (FPREC) , Intent (in)  :: b2bvec(:,:,:)! Ndim x Nbeads x Nbeads
    Real (FPREC) , Intent (in)  :: dr(:,:)      ! Nbeads x Nbeads
    Real (FPREC) , Intent (out) :: Fs(:,:)      ! Ndim x Nbeads 
!!! Purpose:
    ! computes the net spring force on a bead due to 
    ! connector springs of neighbouring beads

    Integer nu
    Real (FPREC)  r ! = |Q| / sqrtb
    Real (DOBL)  Fcnu(Ndim), Fcnu1(Ndim), nonhook

    Fs  = 0.0

    Fcnu1 = 0.0 ! Fc of nu-1

    Do nu = 1,N-1 ! over each of the connector

       r = dr(nu,nu+1)/sqrtb ! only upper diagonal
       If (Abs(r - 1.0) < MYEPS) r = 1 - MYEPS

       Call force_sans_hookean(sptype,r,nonhook)

       ! connector force from 1 to 2 (nu to nu+1)
       Fcnu = b2bvec(:,nu,nu+1) * nonhook

       ! total spring force on a bead
       Fs(:,nu) = Fcnu - Fcnu1
       Fcnu1 = Fcnu ! save for the next bead
    End Do

    ! and the last bead
    Fs(:,N) = - Fcnu1

  End Subroutine Spring_force


  Subroutine solve_implicit_r(sptype,dtby4,gama,r,ff)
    Integer, Intent (in) :: sptype
    Real (FPREC) , Intent (in) :: dtby4, gama
    Real (FPREC) , Intent (out) :: r
    Real (DOBL) , Intent (out) :: ff

    !! Purpose
    ! solves the equation r ( 1 + dt/4 ff ) == gama
    ! where ff depends on spring force law F = H Q ff,
    ! for r and returns, r and ff(r)

    Real (FPREC)  coeff(4),denom

    ! set up the polynomial equation (cubic), and the guess for
    ! gama >> 1, obtained by the asymptotic behaviour

    coeff(4) = 1.0

    Select Case (sptype)
    Case (HOOK) ! Hookean
       r = gama/(1+dtby4)
       ff = 1
       Return

    Case (FENE) ! FENE
       coeff(1) = gama
       coeff(2) = -(1. + dtby4)
       coeff(3) = -gama
       r = 1 - dtby4/2/gama 

    Case (ILC) ! ILC
       denom = 3. + dtby4
       coeff(1) = 3. * gama / denom
       coeff(2) = -(3. + 3.* dtby4) / denom
       coeff(3) = -3. * gama / denom
       r = 1 - dtby4/3/gama

    Case (WLC) ! WLC
       denom = 2*(1.5 + dtby4)
       coeff(1) = -3 *gama/denom
       coeff(2) =  3 * (1. + dtby4 + 2.*gama) / denom
       coeff(3) = -1.5 * (4. + 3.*dtby4 + 2.*gama) / denom
       r = 1 - Sqrt(dtby4/6/gama)

    End Select

    ! the Hookean guess is same for all the force laws
    If (gama < 1.0) Then 
       r = gama/(1.+dtby4)  
    End If

    ! all the common forces laws yeild a cubic in the implicit form
    ! polish the guessed root by a Newt-Raph for polynomials

    ! for WLC, the newt raph, does not converge to the correct root,
    ! for gama > 100 , therefore simply ignore polishing
    If (sptype == WLC .And. gama > 100) Then
       r = r
    Else 
       Call polish_poly_root(coeff,r,1e-6)
    End If


    Call force_sans_hookean(sptype,r,ff)

  End Subroutine solve_implicit_r


  Subroutine Excluded_volume_force(N,b2bvec,dr,Fev)
    Use bspglocons

!!! This routine returns excluded volume force
    Integer, Intent (in) :: N
    Real (FPREC) , Intent (in)  :: b2bvec(:,:,:)! Ndim x Nbeads x Nbeads
    Real (FPREC) , Intent (in)  :: dr(:,:)      ! Nbeads x Nbeads
    Real (FPREC) , Intent (out) :: Fev(:,:)     ! Ndim x Nbeads 


    Integer nu,mu
    Real (FPREC)  Fpair12(Ndim)

    Fev  = 0.0

    If (zsbyds5 > 0.0) Then
       ! caution donot use forall here, i means F12 will b evaluated for 
       ! all the values of mu,nu and then assigned, which is not what we 
       ! want.  The rule to use forall statements is that lhs must be 
       ! unique for-all the values of loop index, and RHS must not depend 
       ! on any other index of lhs, other than the lhs itself
       ! ie, we can have  foo(mu) = foo(mu) + bar(nu)
       ! but not  foo(mu) = foo(mu+3) + bar(nu)
       Do nu = 2,N
          Do mu  = 1,nu-1 
             ! force between the pair mu,nu, since the EV force is
             ! repulsive, we take it to be in the opposite direction
             ! of the connector vector mu -> nu
             ! convention followed: force is positive for attraction
             Fpair12 = - b2bvec(:,mu,nu) & 
                  * zsbyds5*Exp(-dr(mu,nu)*dr(mu,nu) * pt5bydssq) 
             Fev(:,mu) = Fev(:,mu) + Fpair12
             Fev(:,nu) = Fev(:,nu) - Fpair12
          End Do
       End Do
    End If
  End Subroutine Excluded_volume_force

  Subroutine preavDnB(N,hs)

    Integer , Intent (in) :: N
    Real (FPREC), Intent (in) :: hs

    Integer mu,nu,i,j,k
    Integer lwork, info
    Real (FPREC)  x, erfx, Hbar(N,N), Hp(N*(N+1)/2)  &
         , Evec(N,N), Eval(N) &
         , Bpreav(N,N)

    Real (FPREC), Allocatable :: work(:)


    Real (FPREC), save :: hsave=0,nsave=0

    if (hsave ==  hs .and. nsave == N) Then
       Return  ! to avoid recalculation of Bpmat and Dpmat
    end if

    If (.Not. Allocated(Bpmat) ) Then
       Allocate(Bpmat(Ndim,N,Ndim,N),  &
            Dpmat(Ndim,N,Ndim,n))
    Else if (nsave /= N) Then
       Deallocate(Bpmat,Dpmat)
       Allocate(Bpmat(Ndim,N,Ndim,N),  &
            Dpmat(Ndim,N,Ndim,n))
    End If

    hsave = hs
    nsave = N

    !-- Calculate the preaveraged diffusion tensor (only upper diagonal)
    If (hs == 0)  Then
       Write (*,*) 'Hstar = 0, called in preavDnB, check call sequence'
       Return
    end If

    Hbar = 0

    Do mu = 1,N
       Do nu = mu+1,N
          x = sqrt(2*Pi/(nu-mu)) * hs

          Select Case (FPREC)
          Case(SNGL)
	     erfx = erf(x) ! plain fortran (calling libc ?)
             !call vserf(1,x,erfx) ! Intel fortran
             !erfx = erf(x) ! IBM XLF
          Case(DOBL)
	     erfx = erf(x) ! plain fortran (calling libc ?)
             !call vderf(1,x,erfx) ! Intel fortran
             !erfx = derf(x) ! IBM XLF
             !--
             !  Note that with IBM, this line needs to be commented out
             !  for single precision, the compiler does not ignore the
             !  derf(real x) call
             !--
          end Select

          !-- For RPY tensor (Fixman, JCP 78, 1594)
          Hbar(mu,nu) =  erfx - (1-exp(-x*x))/sqrt(Pi)/x
       end Do

       Hbar(mu,mu) = 1
    end Do


    !-- eigensystem of Hbar

    !!*** IBM (ESSL) only ***!!
!!$    !-- 
!!$    !  ESSL: convert to upper packed storage mode, see example in
!!$    !   http://publib.boulder.ibm.com/clresctr/docs/essl/essl_42/
!!$    !   200410/am501402/am501402186.html
!!$    !--
!!$    k=0
!!$    Do j=1,N
!!$       Do i=1,j
!!$          k = k+1
!!$          Hp(k) = Hbar(i,j)
!!$       end Do
!!$    end Do
!!$
!!$    lwork = 0 ! automatic allocation
!!$    Select Case (FPREC)
!!$    Case(SNGL)
!!$       call sspev(21,Hp,Eval,Evec,N,N,work,lwork)
!!$    Case(DOBL)
!!$       call dspev(21,Hp,Eval,Evec,N,N,work,lwork)
!!$    end Select

    !!*** LAPACK  ***!!

    Evec = Hbar
     lwork = 3*N
     Allocate(work(lwork)) ! not required for ESSL

     Select Case (FPREC)
     Case(SNGL)
        call ssyev('V','U', N, Evec, N, Eval, work, lwork, info)
     Case(DOBL)
        call dsyev('V','U', N, Evec, N, Eval, work, lwork, info)
     end Select

     if (info < 0) Then
        Write (*,*) 'Error in syev, illegal value for parameter no',&
              -info
        Stop
     end if
      if (info > 0) Then
        Write (*,*) 'Error in syev, failed to converge:', info
        Stop
     end if


    Bpreav = 0
    Forall (mu=1:N) 
       Bpreav(mu,mu) = sqrt(Eval(mu))
    end Forall

    Bpreav = Matmul(Evec, Bpreav)
    Bpreav = Matmul(Bpreav, Transpose(Evec))


    Dpmat = 0
    Bpmat = 0

    Do mu = 1,N
       Do nu = mu,N
          Do i = 1,Ndim
             Dpmat(i,mu,i,nu) = Hbar(mu,nu)
             !-- for Dpmat only upper diagonal reqd, but full for Bpmat
             Bpmat(i,mu,i,nu) = Bpreav(mu,nu)
             Bpmat(i,nu,i,mu) = Bpreav(mu,nu)
          end Do
       end Do
    end Do


!!$  write (*,*) 'Hbar'
!!$  write (*,81) Hbar
!!$  write (*,81) 'Eval'
!!$  write (*,81) Eval
!!$  write (*,81) 'Evec'
!!$  write (*,81) Evec
!!$  write (*,*) 'Bpreav'
!!$  write (*,81) Bpreav
!!$  stop

!81  Format (<N>(G10.4,2x))
  end Subroutine preavDnB


End Module BSpModel


Subroutine Thermalise_Chain(NBeads, R_Bead, spring_type, &
     Ntcur, Ntsteps, Delts, &
     Hstar, Zstar, Dstar, L0s, FHI, &
     TrajType, nseed, &
     Ntcsamp, Ntsamp, propsums, Nsdone)  ! last line are optional
  Use bspglocons
  Use bspglovars, Only : Imploop_tol
  Use bspintfacs
  Use bspmodel
  Use csputls

  Implicit None

  Integer, Intent(in) :: NBeads
  Real (FPREC) , Intent (inout), Dimension(:,:) :: R_Bead 
  ! Initial positions of Beads

  Integer, Intent (in) :: spring_type ! Spring force law type
  ! Hookean = 1, FENE = 2, ILC = 3, WLC = 4

  Integer(k4b) , Intent (in) :: Ntcur ! Current time steps completed
  Integer , Intent (in) :: Ntsteps ! time steps to be advanced to

  Real (FPREC) , Intent (in) :: Delts  ! Integration time interval specs

  Real (FPREC) ,  Intent (in) :: Hstar         ! HI parameter (non-dim)
  Real (FPREC) ,  Intent (in) :: Zstar,Dstar   ! EV parameters (non-dim)
  Real (FPREC) ,  Intent (in) :: L0s           ! finite ext., param sqrt(b)
  Logical , Intent (in) :: FHI  ! fluctuating/preaveraged HI 

  Integer, Intent (in) :: Trajtype ! EQuilibriation, PRoduction

  Integer (k4b), Intent (inout) :: nseed ! seed for rnd number

  Integer, Intent (in), Optional :: Ntcsamp, Ntsamp ! Time steps per sample
  Real (FPREC), Intent (inout), Optional :: propsums(:) ! prop sums for avgs
  Integer, Intent (inout), Optional :: Nsdone ! samples taken


!!!_________________________________________________________|
  !     This driver uses:                                              
  !     Predictor-corrector        scheme for excluded volume          
  !     Fully implicit scheme for the spring force                     
  !           based on Ottinger's and Somasi et al.'s suggestions      
  !     Rotne-Prager-Yamakawa HI tensor with Chebyshev polynomial     
  !         approx. using Fixman's approximation to get conditioner    
  !     Narrow Gaussian EV potential                                   
!!!_________________________________________________________|

  ! external BLAS functions
  Real (SNGL) snrm2, sdot
  Real (DOBL) dnrm2, ddot

  Integer Ntime
  Real (FPREC)  fd_err, fd_err_max, om1, om2, rs,&
       rsbyhs, hsbyrs, hsbyrs2, L_max, L_min, d_a, d_b,&
       gama_mag, lt_err, dtsby4, pcerr

  Real (FPREC),  dimension(Ndim) :: cofm !, ddisp 

  ! various scratch storages 
  Real (FPREC),  dimension (Ndim,NBeads) ::  &
       R_pred, DR_pred, Ups_pred, R_pred_old, R_corr, Rtemp, delta_R &
       , F_ev, F_spring, F_tot & ! forces
       , DelS &
       , X_0, X_l_1 &
       , X_l, X_lp1 

  Real (FPREC) &
       b2b_sup(Ndim,Nbeads,Nbeads) & ! bead to bead vector (super diagonal)
       , deltaR_sup(Nbeads,Nbeads) & ! bead to bead dist (super diagonal)
       ! symmetric diffusion tensor, only upper diagonal elements stored
       , Diffusion_sup(Ndim,NBeads,Ndim,NBeads) &
       ! another scratch storage for diffusion tensor
       , Dp_sym_up(Ndim,NBeads,Ndim,NBeads) &
       , kappa(Ndim,Ndim)&  !flow field divergence
       ! for Chebyshev polynomials
       , cheba(0:MAXCHEB), gama_mu(Ndim) &
       , diffD(Ndim,NBeads-1,Ndim,NBeads) &
       , diffUps(Ndim,NBeads-1)

  ! some connector forces and related requiring higher precision
  Real (DOBL)  ff,  &
       F_con_mu(Ndim), &    ! connector force between beads
       F_con_mu1(Ndim)    ! connector force, save for mu-1 

  Logical fd_check, ncheb_flag

  !     Other definitions...
  Integer i, j, l, mu, nu, ncheb, ncheb_old, lt_count
  Real (FPREC)  temp1, sqrt2inv, TPI, RPI, TRPI, &
       C1, C2, C3, C4, C5, C6, C7, RPI3by4 , r , rinv, trdiff


  !     Definitions for the BLAS-2 routine "dsymv"
  Integer Ndof, lda, ldadiff, incx, incy
  Real (FPREC)  alpha, beta

  Real (FPREC)  pinst(Ndim,Ncorrel)


 
       


  sqrtb = L0s ! sqrtb will be needed by other modules

  ! for EV
  zsbyds5 = Zstar/(Dstar**5.0)
  pt5bydssq = 0.5/(Dstar*Dstar)


  sqrt2inv = 1.0/(2.0**0.5)
  TPI = 2.0*PI
  RPI = Sqrt(PI)
  TRPI = 2.0*RPI
  C1 = TPI/3.0
  C2 = 1.0
  C3 = 0.75/RPI*0.375 !C3 = 9./32/RPI
  C4 = 8.0*PI
  C5 = 14.14855378
  C6 = 1.21569221
  C7 = 0.09375/RPI
  RPI3by4 = RPI*0.75


  incx = 1
  incy =1
  Ndof = Ndim*NBeads  ! degrees of freedom
  lda = Ndof

  fd_err_max = 0.0025 
  pcerr = 0.0


  delta_R = 0.0


  Diffusion_sup = 0.0
  ! diagonal elements
  Forall (mu = 1:Nbeads, i = 1:Ndim )
     Diffusion_sup(i,mu,i,mu) = 1.0
  End Forall

  !-- Calculate the preaveraged diffusion tensor (only upper diagonal)
  If (.not. FHI .and. Hstar > 0 .and. TrajType /= EQ) Then
     call preavDnB(NBeads, Hstar) ! calculate Bpmat, Dpmat (Bspmodel)
     !--debug print *, 'calculated DnB'
  end If
  
    
  
!!$  ddisp = 0




  !     Get initial estimates of L_max and L_min using
  !     Kroeger et al's approx. of Zimm theory's preds.;
  !     the number of Chebyshev terms, and the Chebyshev coefficients
  L_max = 2*(1 + PI*Sqrt(Dble(NBeads))*Hstar)
  L_min = (1 - 1.71*Hstar)/2     ! Note Hstar < 0.58
  ncheb = Int(Sqrt(L_max/L_min)+0.5)+1
  ncheb_flag = .False.
  d_a = 2/(L_max-L_min)
  d_b = -(L_max+L_min)/(L_max-L_min)
  Call chbyshv(ncheb, d_a, d_b, cheba)
  DelS = 0.0      

  !___________________________________________________________________|
  !                The Time integration loop begins here              |
  !___________________________________________________________________|


  Overtime: Do Ntime = 0,Ntsteps-1

     ! find the center of mass
     cofm = 0.0
     Do mu = 1, NBeads
        cofm = cofm + R_bead(:,mu)
     End Do

     cofm = cofm/NBeads

     ! shift the origin to the current center of mass
     Forall (mu = 1: NBeads)
        R_Bead(:,mu) = R_Bead(:,mu) - cofm
     End Forall

     ! calculate the bead to bead vector and its magnitude
     ! for use in the forces
     Call b2bvector_sym_up(NBeads,R_Bead,b2b_sup)
     Call modr_sym_up(NBeads,b2b_sup,deltaR_sup)


     ! change the nature of the spring, in the specific subroutine
     ! in Spring_force(), in module Bead_spring_model
     Call Spring_force(spring_type,NBeads,b2b_sup,deltaR_sup,F_spring)
     Call Excluded_volume_force(NBeads,b2b_sup,deltaR_sup,F_ev)


     ! Therefore, total force on each bead...    
     F_tot = F_spring + F_ev


     !        Calculate the RPY diffusion tensor
     If (Hstar > 0) Then
        If(FHI) Then ! Fluctuating HI
        Do nu = 2,NBeads
           Do mu = 1,nu-1
              rs = deltaR_sup(mu,nu)
              hsbyrs = Hstar/rs
              rsbyhs = rs/Hstar
              hsbyrs2 = hsbyrs*hsbyrs
              If (rsbyhs >= TRPI) Then
                 om1 = RPI3by4*hsbyrs*(1.0+C1*hsbyrs2)
                 om2 = RPI3by4*(hsbyrs/(rs*rs))*(1-TPI*hsbyrs2)
              Else
                 om1 = 1-C3*rsbyhs
                 om2 = C7*rsbyhs/(rs*rs)
              End If

              ! store only the upper triangular in the mu,nu index.
              Call tensor_prod_of_vec(om2,b2b_sup(:,mu,nu), &
                   Diffusion_sup(:,mu,:,nu) )

              ! additional factor for diagonal elements, in each mu,nu
              Forall (i = 1:Ndim )
                 Diffusion_sup(i,mu,i,nu) = Diffusion_sup(i,mu,i,nu) + om1
              End Forall
           End Do
        End Do
     Else
        Diffusion_sup = Dpmat ! Preaveraged HI
     End If

     
     End If ! Hstar


     !--- The following will be used for the predictor step as well as

     !    for calculating 1/4 D.F for diffusivity
     If (Hstar > 0) Then
        alpha = 0.25
        beta = 0.0
        ! Assigns DR_pred <- 0.25*D.F      
        Select Case (FPREC)
        Case(SNGL)
           Call ssymv('U', Ndof,alpha,Diffusion_sup,lda, F_tot,incx, &
                beta,DR_pred,incy)   
        Case(DOBL)
           Call dsymv('U', Ndof,alpha,Diffusion_sup,lda, F_tot,incx, &
                beta,DR_pred,incy)   
        End Select

     Else
        DR_pred = 0.25 * F_tot
     End If



     !!***  Sampling  ***!!                                  

     If (TrajType == PR) Then

        !*** Time correlation functions
        If ( Mod(Ntcur+Ntime, Ntcsamp) == 0) Then
           !-- End-to-End Distance
           pinst(:,1) = R_Bead(:,NBeads) - R_Bead(:,1)

           !-- 1/4 D . F for calculating long time diffusivity
           pinst(:,2) = Sum(DR_pred,DIM=2)

           !-- Stress (note that origin is at the cofm already)
           Do i = 1,Ndim
              j = Mod(i,Ndim)+1
              pinst(i,3) = Sum(R_Bead(i,:)  * F_tot(j,:))
              pinst(i,4) = Sum(R_Bead(j,:)  * F_tot(i,:))
           End Do

           Call sample_tcorr(SAMPLE, pval=pinst) 
           !order of args changed, so include keyword=
        end if
        

        !*** Averages
        If ( Mod(Ntcur+Ntime, Ntsamp) == 0) Then
           Nsdone = Nsdone +1


        
           Call eqchain_props (NBeads, R_Bead, F_tot, propsums) 

           ! Kirkwood-Diffusivity from trace of diffusion tensor
           ! obtain the inverse r(m,n)

           rinv = 0
           Do nu = 2,NBeads
              Do mu = 1,nu-1
                 rinv = rinv + 1/deltaR_sup(mu,nu)
              End Do
           End Do

           propsums(4) = propsums(4) + 2*rinv /NBeads/NBeads


           ! Kirkwood-Diffusivity from trace of diffusion tensor


           trdiff = 0
           Do mu = 1,Nbeads-1
              !-- nu > mu
              Do nu = mu+1,NBeads
                 Do i=1,Ndim
                    trdiff = trdiff + 2*Diffusion_sup(i,mu,i,nu)
                 End do
              end Do
           end Do
           
           Do mu = 1,Nbeads
              !-- nu = mu
              Do i=1,Ndim
              trdiff = trdiff + Diffusion_sup(i,mu,i,mu)
	      End Do
           end Do
        

           propsums(5) = propsums(5) + trdiff/4/Ndim/NBeads/NBeads

        End If  !whether to sample
     End If !for production runs

     !____________________________________________________________________|
     !    Chebyshev polynomial approximation of DelS begins here          |
     !____________________________________________________________________|

     !Generate the random vector X_0
     X_0 = 0.0
     Call ran_1(0,Ndof, X_0, nseed)
     X_0 = X_0 - 0.5
     X_0 = (X_0*X_0*C5 + C6)*X_0! Element-wise multiplications

     ! Has check for deviation from fluctuation-dissipation theorem
     ! been performed?
     fd_check = .False.

     If (Hstar > 0) Then          

        If (FHI) Then ! Fluctuating HI
        FDloop: Do
           ! Update DelS vector
           DelS = cheba(0) * X_0

           ! Shift the D matrix
           Dp_sym_up = d_a * Diffusion_sup

           ! diagonal elements
           Forall (mu = 1:Nbeads, i = 1:Ndim )
              Dp_sym_up(i,mu,i,mu) = Dp_sym_up(i,mu,i,mu) + d_b
           End Forall

           ! Calculate the second Chebyshev vector            
           X_l_1 = X_0
           alpha = 1.0
           beta = 0.0

           ! BLAS2 symmetric matrix-vector multiplication
           Select Case (FPREC)
           Case(SNGL)
              Call ssymv('U', Ndof, alpha,Dp_sym_up,lda,  X_l_1,incx, &
                   beta,X_l,incy)  
           Case(DOBL)
              Call dsymv('U', Ndof, alpha,Dp_sym_up,lda,  X_l_1,incx, &
                   beta,X_l,incy)  
           End Select


           ! Update DelS vector            
           DelS = DelS + cheba(1)*X_l

           Do l = 2,ncheb
              alpha = 2.0
              beta = 0.0
              Select Case (FPREC)
              Case(SNGL)
                 Call ssymv('U', Ndof, alpha,Dp_sym_up,lda, X_l,incx,  &
                      beta,X_lp1,incy)
              Case(DOBL)
                 Call dsymv('U', Ndof, alpha,Dp_sym_up,lda, X_l,incx,  &
                      beta,X_lp1,incy)
              End Select

              X_lp1 = X_lp1-X_l_1
              X_l_1 = X_l

              ! The l-th Chebyshev vector
              X_l = X_lp1

              ! Update DelS vector               
              DelS = DelS + cheba(l)*X_l

           End Do

           ! Calculate the deviation from 
           ! the fluctuation-dissipation theorem

           If (.Not.fd_check) Then

              !Use D:X_0X_0 = X_0.D.X_0 
              !Get D.X_0 first, and then get the dot product of X_0 with the
              ! resulting vector. X_l is reused 

              alpha = 1.0
              beta = 0.0

              Select Case (FPREC)
              Case(SNGL)
                 fd_err = sdot(Ndof, DelS, 1, DelS, 1) ! BLAS-1 function
                 Call ssymv('U', Ndof,alpha,Diffusion_sup,lda, &
                      X_0,incx, beta,X_l,incy)       
                 temp1 = sdot(Ndof, X_0, 1, X_l, 1)
              Case(DOBL)
                 fd_err = ddot(Ndof, DelS, 1, DelS, 1) ! BLAS-1 function
                 Call dsymv('U', Ndof,alpha,Diffusion_sup,lda, &
                      X_0,incx, beta,X_l,incy)       
                 temp1 = ddot(Ndof, X_0, 1, X_l, 1)
              End Select


              fd_err = Abs((fd_err - temp1)/temp1)
           End If

           If ((fd_err <= fd_err_max) .Or. fd_check) Then
              Exit FDloop
           Else 
              !  If the fd_check has been performed and deviation is large
              !  recaclulate the number of Chebyshev terms required.
              !  First, get the maximum and minimum eigen values of D
              !  using Fixman's suggestion.
              fd_check = .True.
              Call maxminev_fi(Ndof, Diffusion_sup, L_max, L_min)
              ncheb_old = ncheb
              ncheb = Int(Sqrt(L_max/L_min)+0.5)+1
              If ((ncheb/ncheb_old) > 2) ncheb_flag = .True.
              If (ncheb > 500) ncheb = 500
              d_a = 2/(L_max-L_min)
              d_b = -(L_max+L_min)/(L_max-L_min)
              Call chbyshv(ncheb, d_a, d_b, cheba)
           End If
        End Do FDloop
     
     Else
     !-- Preaveraged HI
     DelS = Reshape(Matmul( &
           Reshape(Bpmat,(/ Ndim*NBeads, Ndim*Nbeads /) ),  &
           Reshape(X_0,               (/ Ndim*NBeads /) )   &
           ), (/ Ndim,NBeads/) )
  end If 


!!$     write(*,*) 'Dpmat'
!!$     write(*,82) Dpmat
!!$
!!$     write(*,*) 'D_sup'
!!$     write(*,82) Diffusion_sup
!!$
!!$     write(*,*) 'Bpmat'
!!$     write(*,82) Bpmat
!!$
!!$     write (*,*) 'X_0'
!!$     write (*,82) X_0
!!$     write (*,*) 'dels'
!!$     write (*,82) dels


  Else ! Hstar == 0, the free-draining case              
        DelS = X_0                                                
     End If

     DelS = DelS * Sqrt(delts)

     !____________________________________________________________________|
     !        The predictor step                                          |
     !____________________________________________________________________|


     DR_pred = DR_pred * Delts

     ! Eq (14)
     R_pred = R_Bead + DR_pred + sqrt2inv*DelS

     R_pred_old = R_pred


     !____________________________________________________________________|
     !        The first semi-implicit corrector step                      |
     !____________________________________________________________________|

     ! initialise with part of Eq (18)
     Ups_pred = R_Bead + 0.5*DR_pred + sqrt2inv*DelS

     If (Zstar > 0) Then
        ! Calculate distances between beads using R_pred
        ! reuse variables
        Call b2bvector_sym_up(NBeads,R_pred,b2b_sup)
        Call modr_sym_up(NBeads,b2b_sup,deltaR_sup)
        Call Excluded_volume_force(NBeads,b2b_sup,deltaR_sup,F_ev)

        ! Calculate D.FEV using R_pred and update
        If (Hstar > 0) Then
           alpha = 0.125*Delts      ! The prefactor is 1/8 and not 1/4
           beta = 1.0               ! Add to existing
           Select Case (FPREC)
           Case(SNGL)
              Call ssymv('U', Ndof, alpha,Diffusion_sup,lda, F_ev,incx,  &
                   beta,Ups_pred,incy) 
           Case(DOBL)
              Call dsymv('U', Ndof, alpha,Diffusion_sup,lda, F_ev,incx,  &
                   beta,Ups_pred,incy) 
           End Select

        Else
           Ups_pred = Ups_pred + 0.125 * Delts * F_ev 
        End If
     End If

     !        Calculate the 0.5*Delts*K.R_pred vector 

     ! Eq (18) is completely assembled now

     ! Generating the matrix obtained by using the D_nu
     ! operator on the D super-matrix

     ! (in the following comments indices i,j are omitted for clarity)
     !    diffD(mu,nu) = D(mu+1,nu) - D(mu,nu)
     ! this is true only for nu > mu, since D's elements are computed only
     ! for nu >= mu, and the RHS depends on mu+1.  for nu <= mu+1, the RHS
     ! is rewritten in terms of the symmetric matrix D
     !    diffD(mu,nu) = D(nu,mu+1) - D(nu,mu)
     ! so that all the elements of the RHS are the computed ones

     Forall (mu=1:NBeads-1)
        diffUps(:,mu) = Ups_pred(:,mu+1) - Ups_pred(:,mu) 
     End Forall

     Do mu = 1,Nbeads-1
        !-- nu > mu
       Do nu = mu+1,NBeads
          diffD(:,mu,:,nu) = Diffusion_sup(:,mu+1,:,nu)-Diffusion_sup(:,mu,:,nu)
       end Do
       !-- nu <= mu
       Do nu = 1,mu
          diffD(:,mu,:,nu) = Diffusion_sup(:,nu,:,mu+1)-Diffusion_sup(:,nu,:,mu)
       end Do
    end Do
     

     ! Start the loop for solving for connector vectors
     R_corr = R_pred

     ldadiff = Ndim*(Nbeads-1)

     lt_count = 0

     dtsby4 = Delts/4.0


     Keepdoing: Do

        F_con_mu1 = 0.0

        oversprings: Do mu = 1,NBeads-1


           ! the connector force for this spring is obtained from
           ! Fs(mu) = Fc(mu) - Fc(mu-1)

           F_con_mu = F_con_mu1 + F_spring(:,mu) 
           ! note:  F_spring contains forces evaluated with the 
           !          predictor Q for all beads < mu
           !          corrector Q for all beads > mu

           If ( Hstar > 0 ) Then
              gama_mu = 0.125*Delts* &
                   Matmul( &
                   Reshape(diffD(:,mu,:,:),(/ Ndim, Ndim*Nbeads /) ),  &
                   Reshape(F_spring,             (/ Ndim*Nbeads /)) &
                   )
           Else
              gama_mu  = 0.125*Delts*(F_spring(:,mu+1) - F_spring(:,mu))
           End If

           ! the remaining terms on the RHS of Eq.(20)
           gama_mu = gama_mu + diffUps(:,mu) + 0.25 * F_con_mu * Delts 

           Select Case (FPREC)
           Case(SNGL)
              gama_mag = snrm2(Ndim,gama_mu,1) ! BLAS single normal two
           Case(DOBL)
              gama_mag = dnrm2(Ndim,gama_mu,1) ! BLAS single normal two
           End Select


           ! r = Q_nu_mag/sqrt(b) varies from (0,1)
           Call solve_implicit_r(spring_type,dtsby4,gama_mag/sqrtb,r,ff)
           ! implicit solution of Eq.(21)
           !   r ( 1 + deltat/4 ff) -  |gama_mu|/sqrtb == 0
           ! and returns r and ff, the factor in the spring force other 
           ! than hookean F = H Q ff

           If (Abs(r-1.)  <  1e-6) r = 1 - 1e-6 ! TODO, y this

           ! unit vector for Q_mu same as for Gama_mu
           gama_mu = gama_mu/gama_mag 
           gama_mu = gama_mu*r*sqrtb; ! connector vector update

           ! corrector positions vector update
           R_corr(:,mu+1) = R_corr(:,mu) + gama_mu

           ! Connector force 
           F_con_mu1 = gama_mu * ff ! to b used for next spring

           ! update the spring forces, and connector force for this spring
           ! which will b used at the start of the loop
           F_spring(:,mu)   = F_spring(:,mu)   + (F_con_mu1 - F_con_mu)
           F_spring(:,mu+1) = F_spring(:,mu+1) - (F_con_mu1 - F_con_mu)


        End Do oversprings

        Rtemp = R_corr - R_pred

        Select Case (FPREC)
        Case(SNGL)
           lt_err = snrm2(Ndof,Rtemp,1)/snrm2(Ndof,R_Bead,1)
        Case(DOBL)
           lt_err = dnrm2(Ndof,Rtemp,1)/dnrm2(Ndof,R_Bead,1)
        End Select


        lt_count = lt_count + 1

        If ((lt_err < Imploop_tol) .Or. (lt_count > 10*Nbeads)) Exit

        R_pred = R_corr
     End Do Keepdoing

     If (lt_count > 10*NBeads) Write (*,1024) lt_count, Ntcur, Ntime
1024 Format ('Loop exceeded ', I3, ' for start = ', G11.4, ' at step ', G11.4)

     ! Final Major updates, new position vector and time
     R_Bead = R_corr

     Rtemp = R_Bead - R_pred_old

     Select Case (FPREC)
     Case(SNGL)
        pcerr = snrm2(Ndof,Rtemp,1)/NBeads
        pcerr = snrm2(Ndof,Rtemp,1)/snrm2(Ndof,R_Bead,1)
     Case(DOBL)
        pcerr = dnrm2(Ndof,Rtemp,1)/NBeads
        pcerr = dnrm2(Ndof,Rtemp,1)/dnrm2(Ndof,R_Bead,1)
     End Select


     If (ncheb_flag) Then
        ncheb = ncheb_old
        ncheb_flag = .False.
     End If
  End Do    Overtime

!!$  !-- The current cumulative displacement (which includes the jump 
!!$  !-- at the end of the last Brownian step = cofm) is added to the 
!!$  !-- position vector so that it may be used in the next continuation
!!$  !-- when the then cofm would be = ddisp + cofm
!!$  Forall (mu = 1: NBeads)
!!$     R_Bead(:,mu) = R_Bead(:,mu) + ddisp
!!$  End Forall

!82 Format (<NBeads*Ndim>(G10.4,2x))

End Subroutine Thermalise_Chain

