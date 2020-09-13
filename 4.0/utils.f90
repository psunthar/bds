!!! Time-stamp: <utils.f90 10:09, 18 Oct 2005 by P Sunthar>

!_______________________________________
!  Contains general purpose utilities
!_______________________________________

!!! $Log: utils.f90,v $
!!! Revision 1.3  2005/10/18 00:08:23  pxs565
!!! With dotlockfile external file locking program
!!!
!!! Revision 1.2  2005/05/24 01:07:10  pxs565
!!! variable NTOPB
!!!
!!! Revision 1.1  2005/03/16 00:54:06  pxs565
!!! Initial revision
!!!



Subroutine Initial_position(SType,N,L0s,R,myseed)
  Use bspglocons
  Use bspglovars
  Use bspintfacs
  !Use iflport ! for portable intel intrinsic functions, here ran()
  Implicit None

  Integer, Intent (in) :: Stype
  Integer(k4b), Intent (in) :: N
  Integer(k4b), Intent (inout) :: myseed
  Real (FPREC) , Intent (in) :: L0s
  Real (FPREC) , Intent (out), Dimension(:,:) :: R

  Integer, Save :: varseed = 201235

  Integer mu, done
  Real (FPREC) , Parameter :: upbound = 2
  Real (FPREC)  modr,b2b(Ndim),rv,eqpr, rvchk, b

  Call GaussRand(Ndim*N,R,myseed)



  R(:,1) = 0

  ! intrinsic ran() modifies the seed for every call
  ! so use a separate varible seed, so that the ran_1 sequence is
  ! not altered by this variation.
  !varseed = seed 

  b = L0s*L0s


!!$  do rv = 0,0.99,.01
!!$     eqpr = normeqpr(Stype,rv,b);
!!$     write (*,*) rv,eqpr
!!$  end do



  Do mu = 2,N
     b2b = R(:,mu) ! * sqrt(sigma=1), gaussian dist betwn two adj beads


     If (SType .Ne. HOOK) Then

        done = 0;
        Do While (done .Ne. 1)

           !rv = ran(varseed)
	   call Random_Number(rv) ! standards conforming IBM
           eqpr = normeqpr(Stype,rv,b);


           ! use another uniform random number to decide acceptance
           ! ie, accept only eqpr fraction at this r
           !rvchk = ran(varseed)
	   call Random_Number(rvchk) ! standards conforming IBM
           !write (*,*) rv, rvchk , eqpr
           If (rvchk .Le. eqpr/upbound) done = 1
        End Do

        modr = Sqrt(Sum(b2b*b2b))
        b2b = b2b/modr * (rv * L0s)
     End If

     R(:,mu) = R(:,mu-1) + b2b 

     ! in the case of a finitely extensible spring, a gaussian distribution
     ! can lead to incorrect forces some exceptional b2b vectors, we limit
     ! this by altering the distance when preseving the random direction.
     ! and an additional random factor between 0 and 1
     ! A good guess of the factor cud b obtained more rigorously.
  End Do
End Subroutine Initial_position

Subroutine GaussRand(N,GX,seed)
!!! Returns a normalised gaussian random variable with mean zero 
!!! and variance 1: f(r) = exp(-r^2/2) / sqrt(2*Pi)
!!! To obtain a dist with a given sigma use:  sqrt(sigma) * Gx

  Use bspglocons
  Use bspintfacs
  Implicit None

  Integer(k4b), Intent (in) :: N
  Integer(k4b), Intent (inout) :: seed
  Real (FPREC) , Intent (out), Dimension(N) :: GX

  Real (FPREC)  rnd(2), x,y, rsq, fac
  Integer i,n2
  Real (FPREC)  :: rtemp(N+1) 


  Do i = 1,N,2
     rsq = 0.0
     Do While (rsq .Ge. 1.0  .Or.  rsq .Eq. 0.0)
        Call ran_1(0,2,rnd,seed)
        x = 2*rnd(1) -1
        y = 2*rnd(2) -1
        rsq = x*x + y*y
     End Do
     fac = Sqrt(-2*Log(rsq)/rsq);
     rtemp(i) = x*fac
     rtemp(i+1) = y*fac
  End Do

  GX = rtemp(1:N)

End Subroutine GaussRand




Subroutine ran_1(action, n, X, idum)
  Use bspglocons
  Implicit None
  Integer, Intent(in) :: action
  Integer(k4b), Intent(inout) ::idum
  Integer(k4b), Intent(in) :: n
  Real (FPREC) , Intent(inout) :: X(n)
  !---------------------------------------------------------------------c
  !     Random number generator based on procedure given in             c
  !     Numerical Recipes in Fortran 90, Chapter B7, pp. 1141-1143      c
  !     This generator is supposed to have a period of about 3.1E18.    c
  !                                                                     c
  !     The calling sequence is the same as that of RANL,               c 
  !     the BLAS random number generator. The arguments are:            c
  !     n - gives the dimension of vector R                             c
  !     X - 1D array of dimension n; on exit                            c
  !         contains n unifromly distributed random numbers in (0,1)    c
  !     idum |- seeds which are updated on exit                         c
  !---------------------------------------------------------------------c

  Integer(k4b),Parameter :: IA=16807,IM=2147483647,IQ=127773,IR=2836
  Real (FPREC) , Save :: am
  Integer(k4b), Save :: ix=1, iy = -1, k, callcntr=0
  Integer(k4b) count


  Select Case (action)
  Case (RESET)
     ix = 1
     iy = -1
     k = 0
     callcntr = 0

  Case (SHOW)
     !-- pass the value of ix, which can be used to display
     X(1) = ix
     !-- Write (*,*) 'ran_1() status: ', callcntr,ix,iy,idum


  Case DEFAULT
     callcntr=callcntr+1
     If ((idum.Le.0) .Or. (iy.Lt.0)) Then
        am = Nearest(1.0,-1.0)/IM
        iy = Ior(Ieor(888889999,Abs(idum)),1)
        ix = Ieor(777755555,Abs(idum))
        idum = Abs(idum) + 1
     End If

     Do count = 1,n
        ix = Ieor(ix,Ishft(ix,13))
        ix = Ieor(ix, Ishft(ix,-17))
        ix = Ieor(ix,Ishft(ix, 5))
        k = iy/IQ
        iy = IA*(iy-k*IQ)-IR*k
        If (iy.Lt.0) iy = iy + IM
        X(count) = am*Ior(Iand(IM,Ieor(ix,iy)),1)
     End Do

  End Select

End Subroutine ran_1


Subroutine chbyshv (L, d_a, d_b, a)
  Use bspglocons
  Implicit None
  Integer L
  Real (FPREC)  d_a, d_b, a(0:MAXCHEB)
  !---------------------------------------------------------------------c
  !     This routine calculates the Chebyshev coefficients for the      c
  !     Chebyshev polynomial approximation of a square root function.   c
  !     The routine computes L+1 coefficients, which are returned in    c
  !     in the one-dimensional array a, for the approximation of the    c
  !     square root of any variable that lies in the interval           c 
  !     (d_a, d_b).                                                     c
  !---------------------------------------------------------------------c 


  Integer j, k
  Real (FPREC)  xks(0:L)


  !     Calculate the shift factors
  !	d_a = 2/(lmax-lmin)
  !      -d_b/d_a = (lmax+lmin)/2 

  !     Calculate the collocation points
  Do k = 0,L
     xks(k) = Cos(PI*(k+0.5)/(L+1))/d_a -d_b/d_a
  End Do

  !     Calculate the Chebyshev coefficients
  a = 0.0
  Do j = 0,L
     Do k = 0,L
        a(j) = a(j)+(xks(k)**0.5)*Cos(j*(k+0.5)*PI/(L+1))
     End Do
     a(j) = 2.0/(L+1)*a(j)
  End Do
  a(0) = a(0)/2

End Subroutine chbyshv



Subroutine maxminev_fi(n, A, maxev, minev)
  Use bspglocons
  Implicit None
  Integer n
  Real (FPREC)  A(:,:,:,:), maxev, minev
  !---------------------------------------------------------------------c
  !     This routine approximates the maximum and minimum eigen values  c
  !     of a symmetric matrix nxn matrix A using Fixman's suggestion.   c
  !---------------------------------------------------------------------c
  Integer lda, incx, incy, i
  Real (FPREC)  F(n), G(n), alpha, beta
  Real (SNGL) sdot
  Real (DOBL) ddot !! Blas

  F = 1.0
  G = 0.0

  lda = n
  alpha = 1.0
  beta = 0.0
  incx = 1
  incy = 1

  Select Case (FPREC)
  Case(SNGL)
     Call ssymv('U',n,alpha,A,lda,F,incx,beta,G,incy)
     maxev = sdot(n, F, 1, G, 1)
  Case(DOBL)
     Call dsymv('U',n,alpha,A,lda,F,incx,beta,G,incy)
     maxev = ddot(n, F, 1, G, 1)
  End Select


  maxev = 2.0*maxev/n

  !Forall (i = 1:n) F(i) = (-1.0)**i
  Do i = 1,n,2
     F(i)   = -1.
     F(i+1) =  1.
  End Do


  Select Case (FPREC)
  Case(SNGL)
     Call ssymv('U',n,alpha,A,lda,F,incx,beta,G,incy)
     minev = sdot(n, F, 1, G, 1)
  Case(DOBL)
     Call dsymv('U',n,alpha,A,lda,F,incx,beta,G,incy)
     minev = ddot(n, F, 1, G, 1)
  End Select



  minev = minev/2.0/n
  !        write (*,*) maxev, minev

End Subroutine maxminev_fi




Subroutine polish_poly_root(c,xin,atol)
  Use Bspglocons
  Implicit None
  ! use newton raphson to polish the root of a polynomial
  Integer, Parameter :: n=4 ! presently only for cubic
  Real (FPREC) , Intent (in) :: c(n)
  Real (FPREC) , Intent (inout) :: xin
  Real (FPREC) , Intent (in) :: atol

  Real (FPREC)  ::  p, p1,x0,x

  Integer i,iter
  Integer, Parameter :: IMAX = 30

  x = xin

  ! algo from NR: techniques for polishing, roots of polynomials
  Do iter=1,IMAX
     x0 = x
     p =  c(n)*x + c(n-1)
     p1 = c(n)
     Do i=n-2,1,-1
        p1 = p + p1*x
        p  = c(i) + p*x
     End Do

     !if (abs(p1) < atol) then ! f' = 0 is not solvable in newt-raph
     If (Abs(p1) .Eq. 0.0) Then ! f' = 0 is not solvable in newt-raph
        !! for WLC, in this case f'' also = 0, therefore
        p1 = 6 * c(4) ! f''' 
        p1 = 6*p/p1
        ! note sign(a,b) returns sgn(b) * mod(a)
        x = x - Sign(1._FPREC,p1) * Abs(p1)**(1./3.)
        !! we omit considering cases for ILC and FENE, as they
        !! donot seem to have this singularity
     Else
        x = x - p/p1
     End If

     If (Abs(p) .Le. atol) Exit
  End Do
  If (iter>IMAX) Then
     !write (5,*) 'poly-root: loop exceeded'
  End If


  xin = x
End Subroutine polish_poly_root


Subroutine numint(f,dx,nmin,nmax,nord,sumf)
  Use Bspglocons
  Implicit None

  Real (FPREC) , Intent(in), Dimension(:) :: f
  Real (FPREC) , Intent(in) :: dx
  Real (FPREC) , Intent(out) :: sumf
  Integer, Intent (in) :: nmin, nmax,nord



  Integer, Parameter :: stdout=5, stderr=6
  Integer nbound,nint,j

  nbound = Ubound(f,1)
  nint = nmax-nmin

  !write (*,*) 'nbound = ', nbound


  If (nint > nbound .Or. nmin < 1 .Or. nmax > nbound ) Then
     Write(stderr,*) 'Array out of bounds: ', nmin,nmax,nbound
     Return
  End If

  sumf = 0.

  Select Case (nord) 
  Case (1) ! trapezoidal rule
     sumf = 0.5 * (f(nmax) + f(nmin))

  Case(2)
     sumf = 5./12 * (f(nmin) + f(nmax)) + &
          13./12 *(f(nmin+1) + f(nmax-1)) 
  Case (3)
     sumf = 3./8 * (f(nmin) + f(nmax)) + &
          7./6 *(f(nmin+1) + f(nmax-1)) + &
          23./24 *(f(nmin+2) + f(nmax-2))
  End Select

  Do j = nmin+nord, nmax-nord
     sumf = sumf + f(j)
  End Do
  sumf =  dx*sumf

End Subroutine numint

Subroutine meanerr(vec,mean,err)
  Use Bspglocons
  Real (FPREC), Intent(in) :: vec(:)
  Real (FPREC), Intent(out) :: mean
  Real (FPREC), Intent(out) :: err ! standard error of mean

  Integer n

  n = Ubound(vec,1)
  mean = Sum(vec)/n

  ! though the following is correct, 
  ! err = sqrt((sum(vec*vec) - mean*mean*n)/(n-1))
  ! it does not always numerically evaluate to a quantity > 0 
  ! for very small errors, such as in the case of diffusivity.  
  ! so use the sure shot formula.
  err = 2*Sqrt(Sum((vec-mean)*(vec-mean))/(n-1)/n)
  ! the standard error of mean is approximately in the 95% confidence
  ! interval, giving the factor 2 (from the t-distribution)


End Subroutine meanerr


Function normeqpr(Stype, r, b)
  Use bspglocons
  Implicit None
  Integer, Intent (in) :: Stype
  Real (FPREC) , Intent (in) :: b, r
  Real (FPREC)   normeqpr

  Real (FPREC)  expphi
  External expphi

  Integer, Save :: norm_computed = 0
  Real (FPREC)    , Save :: norm_b, bsav

  Integer , Parameter :: Ln = 21 ! 21 always see below for X,W data
  Real (DOBL) X(Ln), W(Ln), scratch(Ln)
  Real (DOBL)  endpts(2)

  Integer n
  Real (FPREC)  xfact,rg, fb , M

  If ( norm_computed .Eq. 0  .And. bsav .Ne. b) Then

     !! generate the gauss abscissa and weights for legendre
     !call gaussq(1, Ln, 0, 0, 0, endpts, scratch, X, W) ;

     !! generated from c code
     X(01) = -0.993752170620389;   W(01) = 0.016017228257774
     X(02) = -0.967226838566306;   W(02) = 0.03695378977085309
     X(03) = -0.920099334150401;   W(03) = 0.05713442542685689
     X(04) = -0.853363364583318;   W(04) = 0.0761001136283793
     X(05) = -0.768439963475678;   W(05) = 0.09344442345603385
     X(06) = -0.667138804197413;   W(06) = 0.1087972991671478
     X(07) = -0.55161883588722;   W(07) = 0.1218314160537286
     X(08) = -0.424342120207439;   W(08) = 0.1322689386333376
     X(09) = -0.288021316802401;   W(09) = 0.1398873947910733
     X(10) = -0.145561854160895;   W(10) = 0.14452440398997
     X(11) = -2.4782829604619e-16;   W(11) = 0.1460811336496907
     X(12) = 0.145561854160895;   W(12) = 0.1445244039899697
     X(13) = 0.288021316802401;   W(13) = 0.1398873947910732
     X(14) = 0.424342120207439;   W(14) = 0.1322689386333371
     X(15) = 0.55161883588722;   W(15) = 0.1218314160537288
     X(16) = 0.667138804197413;   W(16) = 0.1087972991671492
     X(17) = 0.768439963475678;   W(17) = 0.09344442345603389
     X(18) = 0.853363364583317;   W(18) = 0.07610011362837897
     X(19) = 0.920099334150401;   W(19) = 0.05713442542685797
     X(20) = 0.967226838566306;   W(20) = 0.03695378977085323
     X(21) = 0.993752170620389;   W(21) = 0.01601722825777395




     M = 18  ! some large number, obtained by trial and error

     ! limit the upper bound of M
     If (M > b/2) M = b/2

     xfact = Sqrt(2*M/b)

     norm_b = 0
     fb = 0
     Do n=1,Ln
        rg = (X(n) + 1)/2 * xfact
        norm_b = norm_b +  W(n) * expphi(Stype,rg,b) * rg * rg
        fb = fb +  W(n) * expphi(Stype,rg,b) * rg**4
     End Do

     !Write (*,*) '# b = ', b
     !write (*,*) '# computed normb = ', norm_b
     !Write (*,*) '# f(b) = ', b*fb/norm_b

     bsav = b
     norm_computed = 1
  End If

  normeqpr = r*r*expphi(Stype,r,b)/norm_b

End Function normeqpr

Function expphi(Stype, r, b)
!!! return the distribution function (without normalisation)
  Use bspglocons
  Implicit None
  Integer, Intent (in) :: Stype
  Real (FPREC) , Intent (in) :: r, b
  Real (FPREC)   expphi

  Select Case (Stype)
  Case (HOOK)
     Write (*,*) 'Hookean case called in expb: Check the code'
  Case (FENE)
     expphi = (1-r*r)**(b/2)

  Case (WLC)
     expphi = Exp( -b/6.0 * ( 2*r*r + 1/(1-r) -r -1 ));
  End Select
End Function expphi


Function  lam1_th(hs, N)
  Use bspglocons
  Implicit None
  Real (FPREC)  lam1_th
  Real (FPREC) ,  Intent (in) :: hs
  Integer,  Intent (in) :: N


  Real (FPREC) s,b,aj


  b   = 1 - 1.66*hs**0.78;
  s   =   - 1.40*hs**0.78;

  aj = 4*Sin(Pi/2./N)*Sin(Pi/2./N);
  lam1_th = 2./( aj * b / N**s);

End Function lam1_th


!!!===================================================================!!!
!!! This module contains utilities to lock files to prevent other programs
!!! to open and write data into it.  This is done by creating another file in
!!! the directory called <file>.lock, where <file> is the file to be locked.
!!! The module contains routines for locking, unlocking (removing the lockfile
!!! The locking and unlocking is achieved by an external executable 
!!!   dotlockfile
!!! Since an external program is called, make sure the PATH variable
!!! contains the location of dotlockfile
!!! Attempts to lock files from within this program was not entirely 
!!! successfull, and difficult to find where the error occurred.
!!!===================================================================!!!

Module flock_utils
Implicit none

!--a constant number to be added to the unit of the file to be locked
Integer, parameter :: lunit = 3141
Character(40) :: lfile,syscmd

Contains
  
  Function lockfilename(file)
    Character (*) file
    Character (40) lockfilename
    !--string suffix for the lock file
    Character(5), parameter :: sufix = '.lock'

    lockfilename = trim(file) // sufix
  end Function lockfilename
  

  Subroutine islocked(file,lstat)

    Character(*), intent(in) :: file
    Logical, intent(out) :: lstat


    lfile=lockfilename(file)
    Inquire (file=lfile, exist=lstat)

  end Subroutine islocked
  

  Subroutine getlock(file)
  !! Obtain a lock
  !  Either a new one, or wait till the present one is released.
  !  If a stale file (not used for 5 mins), it is deleted
  !  and a fresh lock is obtained (see details in manual page
  !  of dotlockfile
  
    Character(*), intent(in) :: file
    
    Character(8) date
    Character(10) time

    
    lfile = lockfilename(file)

    !Call Date_and_Time(date,time)
    !write(*,*) time,': Waiting to obtain lock'
    !-- 12 retries (390 secs) or stale file (300 secs) which ever is earlier
    syscmd= 'dotlockfile -r 12 ' // lfile
    call system(syscmd)
    !Call Date_and_Time(date,time)
    !write(*,*) time,': Obtained lock'

  end Subroutine getlock


  Subroutine unlockfile(file)

    Character(*), intent(in) :: file

    lfile = lockfilename(file)
    syscmd= 'dotlockfile -u ' // lfile
    call system(syscmd)

  end Subroutine unlockfile
end Module flock_utils
