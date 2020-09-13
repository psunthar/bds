!!! Time-stamp: <dmodules.f90 12:35, 22 May 2005 by P Sunthar>

!________________________________________________________________
!   Global modules and interface headers
!________________________________________________________________

!!! $Log: dmodules.f90,v $
!!! Revision 1.3  2005/05/24 01:07:24  pxs565
!!! variable NTOPB
!!!
!!! Revision 1.2  2005/03/16 00:55:07  pxs565
!!! Preaveraged RPY working
!!!
!!! Revision 1.1  2005/02/25 21:25:26  pxs565
!!! Initial revision
!!!


Module Bspglocons ! Bead Spring simulation, Global constant (parameters)
  Save  ! always use to ensure contents remein unchanged between calls

  !-- Physical
  Integer, Parameter :: Ndim = 3 ! Dimension of simulation
  Integer, Parameter ::  MaxXdist = 101

  !-- Computational dimensions
  Integer, Parameter :: k4b = Selected_int_kind(9)
  Integer, Parameter :: SNGL = Selected_real_kind(4)
  Integer, Parameter :: DOBL = Selected_real_kind(8)

  ! Choose the precision of computation here
  !Integer, Parameter :: FPREC = DOBL
  Integer, Parameter :: FPREC = SNGL

  !-- Math Constants
  Real (FPREC) , Parameter :: PI = 3.14159265358979323846
  Real (FPREC) , Parameter :: TINI = 1e-25
  Real (FPREC) , Parameter :: MYEPS = 1e-6
  Integer, Parameter :: MAXCHEB = 500


  Integer, Parameter :: Hook = 1, FENE = 2, ILC = 3, WLC = 4

  !-- Averages
  Integer, Parameter :: NAvgs = 5
  ! for Eqbm and Prod runs of TrajType
  Integer, Parameter :: EQ = 1, PR = 2

  !-- Time correlation
  Integer, Parameter :: Ncorrel = 4
  Integer, Parameter :: ALLOCAT=0, INIT=1, SAMPLE=2, PRNT=3  &
       , GETDAT=4, SAVDAT=5

  !-- for ran_1()
  Integer, Parameter :: RESET=1, SHOW=2

End Module Bspglocons


Module Bspglovars ! Bead Spring simulation, Global variables
  Use Bspglocons
  Save  ! always use to ensure contents remein unchanged between calls




  Real (FPREC)  Imploop_tol
End Module Bspglovars


Module BlockTcorr
  Use Bspglocons
  Save  ! always use to ensure contents remain unchanged between calls

  Integer, Parameter :: MAXBLKCOR = 2 ! Maximum number of blocks
  Integer NTOPB 

  Real (FPREC) , Parameter :: TDIFFMAX = 3e6

  Integer   Nto(MAXBLKCOR)     &     ! Number of elements filled
       , Iblmc                       ! counter to track MaxBlock
  

  Integer, Allocatable :: Npp(:,:)   ! number of TC sampled
  Real (FPREC), Allocatable ::   PPcorr(:,:,:) ! correl funcs 

End Module BlockTcorr



Module Bspintfacs
  Use Bspglocons

  !---------------------------------------------------------------------c	
  !     Driver subroutines                                              c
  !---------------------------------------------------------------------c

  Interface
     Subroutine Initial_position(Stype,N,L0s,R,seed)
       Use bspglocons
       Integer, Intent (in) :: Stype
       Integer(k4b), Intent (in) :: N
       Real (FPREC) , Intent (in) :: L0s
       Integer(k4b), Intent (inout) :: seed
       Real (FPREC) , Intent (out) :: R(:,:)
       !Real, intent (out) :: R(Ndim,N)
     End Subroutine Initial_position
  End Interface

  Interface
     Subroutine GaussRand(N,GX,seed)
       Use bspglocons
       Integer(k4b), Intent (in) :: N
       Integer(k4b), Intent (inout) :: seed
       Real (FPREC) , Intent (out), Dimension(N) :: GX
     End Subroutine GaussRand
  End Interface


  Interface
     Subroutine Thermalise_Chain(NBeads, R_Bead, spring_type, &
          Ntcur, Ntsteps, Delts, &
          Hstar, Zstar, Dstar, L0s, FHI, &
          TrajType, nseed, Ntcsamp, Ntsamp, propsums, Nsdone)
       Use bspglocons
       Implicit None
       Integer, Intent(in) :: NBeads
       Real (FPREC) , Intent (inout), Dimension(:,:) :: R_Bead 

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

       Integer (k4b), Intent (inout) :: nseed 

       Integer, Intent (in), Optional :: Ntcsamp, Ntsamp 
       Real (FPREC), Intent (inout), Optional :: propsums(:) 
       Integer, Intent (inout), Optional :: Nsdone 

     End Subroutine Thermalise_Chain


  End Interface


  Interface
     Subroutine polish_poly_root(c,xin,atol)
       Use Bspglocons
       Implicit None
       ! use newton raphson to polish the root of a polynomial
       Integer, Parameter :: n=4
       Real (FPREC) , Intent (in) :: c(n)
       Real (FPREC) , Intent (inout) :: xin
       Real (FPREC) , Intent (in) :: atol
     End Subroutine polish_poly_root
  End Interface

  !---------------------------------------------------------------------c	
  !     Utility subroutines                                             c
  !---------------------------------------------------------------------c

  Interface
     Subroutine ran_1(action, n, R, idum)
       Use bspglocons
       Integer, Intent(in) :: action
       Integer(k4b), Intent(in) :: n
       Real (FPREC) , Intent(inout) :: R(n)
       Integer(k4b), Intent(inout) ::idum
     End Subroutine ran_1
  End Interface

  Interface
     Subroutine maxminev_fi(n, A, maxev, minev)
       Use Bspglocons
       Integer n
       Real (FPREC)  A(:,:,:,:), maxev, minev
     End Subroutine maxminev_fi
  End Interface

  Interface
     Subroutine chbyshv (L, d_a, d_b, a)
       Use bspglocons
       Integer L
       Real (FPREC)  d_a, d_b, a(0:MAXCHEB)
     End Subroutine chbyshv
  End Interface

  Interface
     Subroutine cublu
       Use bspglocons
     End Subroutine cublu
  End Interface
  !---------------------------------------------------------------------c
  !     Properties estimators                                           c
  !---------------------------------------------------------------------c	

  Interface
     Subroutine eqchain_props(N, Rbead, F, propsums)
       Use bspglocons
       Implicit None
       Integer, Intent (in) :: N
       Real (FPREC) , Intent (in), Dimension(:,:)  :: Rbead, F
       Real (FPREC) , Intent (inout),Dimension(:) ::propsums
     End Subroutine eqchain_props
  End Interface

  Interface
     Subroutine numint(f,dx,nmin,nmax,nord,sumf)
       Use Bspglocons
       Implicit None

       Real (FPREC) , Intent(in), Dimension(:) :: f
       Real (FPREC) , Intent(in) :: dx
       Real (FPREC) , Intent(out) :: sumf
       Integer, Intent (in) :: nmin, nmax,nord
     End Subroutine numint
  End Interface

  Interface
     Subroutine meanerr(vec,mean,err)
       Use Bspglocons
       Real (FPREC) , Intent(in) :: vec(:)
       Real (FPREC) , Intent(out) :: mean
       Real (FPREC) , Intent(out) :: err
     End Subroutine meanerr
  End Interface

  Interface
     Function normeqpr(Stype, r, b)
       Use Bspglocons
       Implicit None
       Integer, Intent (in) :: Stype
       Real (FPREC) , Intent (in) :: b, r
       Real (FPREC)   normeqpr
     End Function normeqpr
  End Interface

  Interface
     Function  lam1_th(hs, N)
       Use Bspglocons
       Real (FPREC)  lam1_th
       Real (FPREC) ,  Intent (in) :: hs
       Integer,  Intent (in) :: N
     End Function lam1_th
  End Interface



  Interface
     Subroutine Sample_Tcorr(action, sfile, pval)
       Use Bspglocons
       Use BlockTcorr
       Implicit None

       Integer, Intent(in) :: action
       Real (FPREC) , Intent(in), Optional :: pval(Ndim,Ncorrel) 
       Character (*), Intent (in), Optional :: sfile

     End Subroutine Sample_Tcorr

  End Interface

End Module Bspintfacs

