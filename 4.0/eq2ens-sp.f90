! Time-stamp: <eq2ens-sp.f90 07:09, 15 Mar 2008 by P Sunthar>


!!! $Log: eq2ens-sp.f90,v $
!!! Revision 1.0  2005/10/17 23:17:49  pxs565
!!! Initial revision
!!!
!!! Revision 1.4  2005/05/24 01:50:18  pxs565
!!! With variable NTOPB
!!!
!!! Revision 1.3  2005/05/02 01:20:25  sunthar
!!! FHI only (for time scaling)
!!!
!!! Revision 1.2  2005/04/30 00:36:37  sunthar
!!! Based on xcscal RCS:1.4, PPcorr save included
!!!
!!! Revision 1.1  2005/04/29 01:44:38  sunthar
!!! Initial revision
!!!



Program Eqbm2Ensemble

!!!*** ====================   Declarations ======================= ***!!!

  !! *** Load Modules containing global variables and interfaces *** !!

  Use bspglocons
  Use bspglovars
  !Use BlockSums
  Use BlockTcorr
  Use bspintfacs
  Use Flock_utils
  Implicit None

  !!*** Variables  ***!!

!!$  !*** Additional debugging information
!!$  !Logical, parameter :: DEBUG = .TRUE. 
!!$  Logical, Parameter :: DEBUG = .False. 
!!$  Integer, Parameter :: debunit=16
!!$  Real (8) wtime, wtime0, wtimeinit


  Integer mpierr, Nprocs, MyId


  !** Model parameters and variables
  Integer  :: Stype, NBeads
  Real (FPREC) :: hstar, zstar, dstar, sqrtb
  Real (FPREC), Allocatable :: PosVecR(:,:), pv_main(:,:), pv_cv(:,:)

  !** Time-Integration  parameters

  !--
  !  Set the multiple of Rouse relaxation times used for the 
  !  EQuilibriation trajectory, one with equilibrium time step
  !  and another with production run time step
  !--

  Integer, parameter :: EqEqTime = 8, EqPrTime = 1
  Real (FPREC), parameter :: deltseq = 1

  Integer :: Ntraj, Ntrtot, Ntsteps
  Integer (k4b) :: Nitot, Ntdone 
  Real :: fNitot, fNtrtot
  Real (FPREC) :: deltspr


  !** Sampling parameters
  Integer  :: Maxlam, TrajType, Ntsamp, Ntcsamp, Nsdone
  Real (FPREC) :: Epslam, Dctime

  !** File i/o
  Integer, parameter :: SavMins = 15 ! time interval between saves
  Character (len=20) infile, outfile, gavgfile, corfile, gppfile &
       , resfile, contfile, xcfile
  Integer, Parameter ::   stdout=5, inunit=10, resunit=11, corunit = 12 &
       , gavunit=13, gppunit=14, cunit=17, xcunit=19, avunit=20
  Character (len=10), Parameter :: ErString="Error"

  !-- for gfortran
  Character (len=80) :: Format500, Format501, Format600, Format601, &
       Format650, Format651

  !** Trajectories related
  Logical Filexists, ToContinue, FHI, Flocked
  Integer (k4b) nseed, nseed_cv
  Real (FPREC) ::  Delts, Conf_t, dumX(1)

  !** Misc
  Character(10) :: Prop_names(NAvgs), Correl_names(Ncorrel)
  Integer i, j, nprops

  Real (FPREC) tlongest,  t1rouse, t1zimm

  Real (FPREC), Allocatable :: avgs(:)

  Integer ni,  iblock, nblock, gntrajdone, gntrthis, ntrajproc
  Real (FPREC), Allocatable, Dimension(:,:,:) :: xcorravgs,  gxcavgs
  Real (FPREC), Allocatable, Dimension(:,:) :: Xvar, Yvar


  Integer ppcount, ib
  Real (FPREC) ctime  
  Real (FPREC), Allocatable ::  ppbuf(:,:,:), Gppavg(:,:,:,:) 

  Real time_begin, time_sa, time_se, time_end


!!!*** =============== Execution begin ============================ ***!!!

  !  !!*** MPI Intialisation ***!!
  !
  !  Call MPI_Init(mpierr)                    
  !
  !  !-- Obtain the total number of processors
  !  Call MPI_Comm_Size (MPI_Comm_World, Nprocs, mpierr) 
  !
  !  !-- Obtain the rank or ID of the current processor
  !  Call MPI_Comm_Rank (MPI_Comm_World, MyId, mpierr)                

  Nprocs = 1
  MyId = 0


  !!*** Input DATA ***!!
  infile="par-eq.in"
  Open (unit=inunit,file=infile,status="old")
  Read (inunit,*)
  Read (inunit,*) SType, NBeads, hstar, zstar, dstar, sqrtb
  Read (inunit,*)
  Read (inunit,*) Ntraj, fNtrtot, fNitot, deltspr, Imploop_tol
  !^^ Note Imploop_tol is a global variable

  !-- 
  !  It is easier to enter a large number in e format, but
  !  gfotran does not automatically read it into an integer
  !--
  Nitot = fNitot
  Ntrtot = fNtrtot

  Read (inunit,*)
  Read (inunit,*) Epslam, Maxlam
  Close (unit=inunit)



  !!*** Memory Allocation  ***!!

  !-- longest relaxation times obtained from Rouse and Zimm Models
  t1rouse =   0.5/Sin(PI/2/Nbeads)**2
  t1zimm = lam1_th(hstar,NBeads) ! Thurstons approximation

  !-- Number of time steps per sample (for time correl)
  Ntcsamp = Max(t1zimm * epslam/deltspr, 1._FPREC) 

  !--Number of time origins per block (number of points per block)
  NTOPB = Maxlam * t1zimm / (Ntcsamp * deltspr)

  If (NTOPB <= 1) Then
     print *, 'Error: NTOPB=', NTOPB, ' is <=1. Check the input'
     !call MPI_Abort(mpierr)
     stop
  End If




  Allocate(PosVecR(Ndim,NBeads), pv_main(Ndim,NBeads), pv_cv(Ndim,NBeads) )
  Allocate ( &
       avgs(NAvgs)  &
       )

  Allocate ( &
       xcorravgs(NAvgs,Nitot,3) &
       , Xvar(Navgs,Nitot) &
       , Yvar(Navgs,Nitot) &
       )



  call Sample_Tcorr(ALLOCAT)



  !** Variables specific to root processor only

  If (MyId.Eq.0) Then

     Call CPU_Time(time_begin)

     Allocate ( &
          gxcavgs(NAvgs,Nitot,3)  &
                                !-- receive buffer
          , ppbuf(Ncorrel,MAXBLKCOR,0:NTOPB) &  
                                !-- global avg, for main and cv
          , Gppavg(Ncorrel,MAXBLKCOR,0:NTOPB,2) & 
          )

     ! the standard error of mean obtained from t-distribution
     ! for sample mean and sample standard deviation
     Conf_t = 2 ! 95% confidence for degrees of freedom > 20

     Prop_names(1) = "R^2"        ! square of end-to-end vector          
     Prop_names(2) = "Rg^2"       ! sq. Radius of gyration
     Prop_names(3) = "X"          ! Avg Stretch 
     Prop_names(4) = "Rinv"       ! Sum <1/rmn> / N^2
     Prop_names(5) = "TrDiff"  ! Center of mass diffusivity

     Correl_names(1) = "Re"
     Correl_names(2) = "A=D.F/4"
     Correl_names(3) = "Sxy"

  End If


  !!*** Initialisation ***!!

  !** Generate a pseudo random seed
  call GetaSeed(nseed)

  !--debug suntemp nseed = 201271
  print *, 'Initial seed of ', MyId, ' = ', nseed 



  !** Other initialisations

  !-- the current number of time steps completed (not used in this 
  !-- algorithm)
  Ntdone = 0

  !-- Filenames
  gavgfile = "gavgseq.dat"
  gppfile  = "gppcorr.dat"

  resfile = "result.dat" ! only for the Master processor 
  contfile = "continue"
  xcfile = 'xcorr.dat' ! cross correlation function
  corfile = "acorr.dat"

  !--
  !        Initialization variable format expressions
  !         <> language extension not available in gfortran
  !--


  !500 Format ('#', 2(A13, 2x), <2*(Ncorrel-1)>(A13, 2x))
  Write(Format500,"(a,I3,a)") "('#', 2(A13, 2x), ", 2*(Ncorrel-1), "(A13, 2x))"

  !501 Format (2(G13.6, 2x), <2*(Ncorrel-1)>(G13.6, 2x) )
  Write(Format501,"(a,I3,a)") "(2(G13.6, 2x), ", 2*(Ncorrel-1), "(G13.6, 2x) )"



  !600 Format (A, 2x, <NAvgs-1>(A13,2x))
  Write(Format600,"(a,I3,a)") "(A, 2x,", NAvgs-1, "(A13,2x))"

  !601 Format (G13.6, 2x, <NAvgs-1>(G13.6, 2x) )
  Write(Format601,"(a,I3,a)") "(G13.6, 2x, ", NAvgs-1, "(G13.6, 2x) )"

  !650 Format (A, 2x, <2*(NAvgs-1)>(A13,2x))
  Write(Format650,"(a,I3,a)") "(A, 2x, ", 2*(NAvgs-1), "(A13,2x))"
  
  !651 Format (G13.6, 2x, <2*(NAvgs-1)>(G13.6, 2x) )
  Write(Format651,"(a,I3,a)") "(G13.6, 2x, ", 2*(NAvgs-1), "(G13.6, 2x) )"


  !** Continuation run (if possible)

  !-- obtain the # completed trajectories from the disk
  Inquire (file=gavgfile, exist=Filexists)
  If (Filexists) Then
     Open (unit=gavunit,file=gavgfile,status='old')
     Read (gavunit,*) gntrajdone
     Close(gavunit)
  Else
     gntrajdone = 0
  End If

  !-- Ntraj is the number of trajectories for this run
  !-- Ntrtot is the total trajs desired
  Ntraj = Min(Ntraj,Ntrtot-gntrajdone)

  !-- trajectories per processor
  nblock = Max(Ntraj/Nprocs,1) ! Note integer division



  !-- Initialisation configuration of the main variate's PosVec

  Call Initial_position(SType,Nbeads,sqrtb,pv_main,nseed)


  !-- 
  !  Number of time steps: for one sampling production run,
  !  and for sampling time correlation function
  !--
  Ntsteps = t1zimm/deltspr ! approx decorrelation time steps
  Ntsamp = Ntsteps ! sample once after decorrelation

  !-- sampling interval for time-correl
  Dctime = deltspr*Ntcsamp




  !-- averages for the two variates and their cross correlation function
  xcorravgs = 0
  !-- number of trajectories completed in the current run (in each processor)
  ntrajproc = 0


  !-- number of items to be transferred in mpi_reduce
  ppcount = Ncorrel*MAXBLKCOR*(NTOPB+1)

  !-- Global time correlation fn data
  If (MyID ==0) Gppavg = 0

  !-- Start time to save data at regular intervals
  Call CPU_Time(time_sa)

!!!*** =================== Trajectories Begin ======================== ***!!!

  Trajectories: Do iblock = 1,nblock

     !-debug print *,  'Begin trajectory ' , iblock, ' in proc ', MyID

     !!***  Main Variate (with fluctuating HI) ***!!

     !-- 
     !  A fresh trajectory is started from a configuration which
     !  has been thermalised for Nitot steps from the current,
     !  see below
     !--
     PosVecR = pv_main   

     !-- save for future use for control variate below
     pv_cv = PosVecR  

     !-- 
     !  for a variation across trajectories,'cos ran_1 is reset below.
     !  It is possible that nseed remains unaltered across one trajectory,
     !  if sufficient number of calls are not made.
     !--
     nseed = nseed+1  
     nseed_cv = nseed
     Call ran_1(RESET,1,dumX,nseed)


     !** Equilibriate chain configuration 

     !--
     !  An exact synchronisation of the control variate with the main is
     !  not valid for cases other than Hookean in theta solvent.  In all 
     !  other cases, it is essential to equilibriate the control variate
     !  to its equilibrium distribution, starting from the synchronised value.
     !  The sampling is done only after this.  In order to do this both 
     !  the main and the CV are thermalised for a predetermined duration
     !  of time (typically about 10 T_Zimm).  Since we have seen
     !  earlier that the variates remain correlated for more than 300 
     !  T_Zimm, we can be sure that after 10 Tzimm, both will remain 
     !  correlated, and obey the respective equilibrium distribution 
     !  functions.  Note the equilibriation is required only for the 
     !  control variate, however, in order to maintain the same RN 
     !  sequence, the main also needs to go through the same sequence.
     !--


     If (zstar .Ne. 0 .or. SType /= Hook) Then
        TrajType = EQ

        !--debug print *, 'eqbting in main'

        Ntsteps = EqEqTime*t1rouse/deltseq

        !-- Free draining at a larger time step
        Call Thermalise_Chain(NBeads, PosVecR, SType,  &
             Ntdone, Ntsteps, deltseq, &
             0._FPREC, zstar, dstar, sqrtb, FHI, TrajType, nseed)


        !-- Free draining at the production time step
        Ntsteps = EqPrTime*t1rouse/deltspr

        Call Thermalise_Chain(NBeads, PosVecR, SType,  &
             Ntdone, Ntsteps, deltspr, &
             0._FPREC, zstar, dstar, sqrtb, FHI, TrajType, nseed)

     End If

     !--debug Call ran_1(SHOW,1,dumX,nseed)
     !--debug print *, MyID, 'RNG status (main):', dumX

     !**  Production run

     FHI=.True.

     !-- Production trajectory for sampling
     TrajType = PR

     !-- Number of time steps to thermalise
     Ntsteps = t1zimm/deltspr ! approx decorrelation time
     Ntsamp = Ntsteps ! one sample after decorrelation

     !-- Initialise static averages and time corr variables
     Xvar = 0
     avgs = 0
     Call Sample_Tcorr(INIT) ! PPcorr is initialised to 0

     Do ni=1,Nitot

        !--
        !  avgs is a cumulative sum in the following call; to get the 
        !  instantaeneous value it is added here and will be 
        !  subtracted out below. Note that Ntsamp = Ntsteps implies
        !  that only one sample is taken for one call.(at the begining)
        !--
        Xvar(:,ni) = avgs 

        !-- Actual Brownian dynamics done only here!
        Call Thermalise_Chain(NBeads, PosVecR, SType,  &
             Ntdone, Ntsteps, deltspr, &
             hstar, zstar, dstar, sqrtb, FHI, TrajType, nseed, &
             Ntcsamp, Ntsamp, avgs, Nsdone)

        Xvar(:,ni) = avgs - Xvar(:,ni)
     End Do

     !-- Time correlations from sample_tcorr
     Do i=1,Ncorrel
        Where (Npp > 0) &
             PPcorr(i,:,:) = PPcorr(i,:,:) / (Npp / Ncorrel)
        ! since Npp is incremented by Ncorrel for each sample
     End Do

     !Call MPI_Barrier (MPI_Comm_World, mpierr)
     !Call GMPI_Redsum(PPcorr,ppbuf,ppcount,mpierr)
     ppbuf = PPCorr


     !-- 
     !  transfer to the global average of main as ppbuf will be
     !  overwritten 
     !--
     If (MyID == 0) Gppavg(:,:,:,1) = Gppavg(:,:,:,1) + ppbuf


     !-- 
     !  save for the next trajectory. Note: PosVecR will b 
     !  overwritten in the control variate below
     !--
     pv_main = PosVecR 

     !-- 
     !  display the RNG status, use for development stage, to compare
     !  with the control variate
     !--
     !--debug Call ran_1(SHOW,1,dumX,nseed)
     !--debug print *, MyID, 'RNG status (main):', dumX




     !!***  Repeat for the control variate (preaveraged HI) ***!!

     !-debug print *, 'control variate'
     FHI = .False.

     !-- Restore saved values for the control variate
     PosVecR = pv_cv

     nseed = nseed_cv
     Call ran_1(RESET,1,dumX,nseed)

     !** Equilibriate chain configuration 

     !--
     !  See notes in the FHI equilibriation section above
     !--


     If (zstar .Ne. 0 .or. SType /= Hook) Then
        TrajType = EQ

        !--debug print *, 'eqbting in cv'

        Ntsteps = EqEqTime*t1rouse/deltseq

        !-- Free draining at a larger time step
        Call Thermalise_Chain(NBeads, PosVecR, Hook,  &
             Ntdone, Ntsteps, deltseq, &
             0._FPREC, 0._FPREC, 1._FPREC, sqrtb, FHI, TrajType, nseed)


        !-- Free draining at the production time step
        Ntsteps = EqPrTime*t1rouse/deltspr
        Call Thermalise_Chain(NBeads, PosVecR, Hook,  &
             Ntdone, Ntsteps, deltspr, &
             0._FPREC, 0._FPREC, 1._FPREC, sqrtb, FHI, TrajType, nseed)

     End If

     !--debug Call ran_1(SHOW,1,dumX,nseed)
     !--debug print *, MyID, ': RNG status (cv):', dumX


     !--debug print *, 'Production in cv'

     !** Production run

     !-- Production trajectory for sampling
     TrajType = PR

     !-- Number of time steps to thermalise
     Ntsteps = t1zimm/deltspr ! approx decorrelation time
     Ntsamp = Ntsteps ! one sample after decorrelation

     !-- Initialise static averages and time corr variables
     Yvar = 0
     avgs = 0
     Call Sample_Tcorr(INIT)

     Do ni=1,Nitot
        Yvar(:,ni) = avgs 

        !--debug print *, 'Here', ni
        !-- 
        !  Note that for the control variate, we use Hookean springs
        !  at theta conditions with preaveraged HI.  
        !  'sqrtb' is therefore ignored in the  procedure
        !--
        Call Thermalise_Chain(NBeads, PosVecR, Hook,  &
             Ntdone, Ntsteps, deltspr, &
             hstar, 0._FPREC, 1._FPREC, sqrtb, FHI, TrajType, nseed, &
             Ntcsamp, Ntsamp, avgs, Nsdone)

        Yvar(:,ni) = avgs - Yvar(:,ni)

     End Do

     !-- debug print *, 'End pr cv'   

     !-- Time correlations from sample_tcorr
     Do i=1,Ncorrel
        Where (Npp > 0) &
             PPcorr(i,:,:) = PPcorr(i,:,:) / (Npp / Ncorrel)
        ! since Npp is incremented by Ncorrel for each sample
     End Do

     !Call MPI_Barrier (MPI_Comm_World, mpierr)
     !Call GMPI_Redsum(PPcorr,ppbuf,ppcount,mpierr)
     ppbuf = PPcorr

     !-- 
     !  transfer to the global average of CV as ppbuf will be
     !  overwritten 
     !--
     If (MyID == 0) Gppavg(:,:,:,2) = Gppavg(:,:,:,2) + ppbuf



     !-- 
     !  display the RNG status, use for development stage, to compare
     !  with the main
     !--
     !--debug Call ran_1(SHOW,1,dumX,nseed)
     !--debug print *, MyID, ': RNG status (cv):', dumX


     !-- store for the mean and cross-correlation 
     xcorravgs(:,:,1) =  xcorravgs(:,:,1) +  Xvar
     xcorravgs(:,:,2) =  xcorravgs(:,:,2) +  Yvar
     xcorravgs(:,:,3) =  xcorravgs(:,:,3) +  Xvar*Yvar

     !-- number of trajectories completed in this processor
     ntrajproc = ntrajproc + 1

     !-- 
     !  While it is expected that for equilibrium runs, atleast,
     !  each of the processors will take the same amount of time 
     !  irrespective of the starting configuration,  we however impose a
     !  barrier here to ensure synchronised computation.  This is because
     !  it is possible (in principle) that one of the processors 
     !  completes all the 'nblock' trajectories within the time 
     !  saving interval and exits the loop.  This would mean that the
     !  consolidation of data across processors done in 'xcorrsav()' 
     !  will be from different calls, one following and one outside the loop.
     !  The "fast" processors will end up calling xcorrsav only once!
     !--
     !Call MPI_Barrier (MPI_Comm_World, mpierr)
     !--debug print *, 'END trajectory : ', iblock, Myid

     !--End time for writing to disk at regular intervals
     Call CPU_Time(time_se)

     !-- Save data at regular intervals
     If ( (time_se-time_sa)/60  >= SavMins) Then
        !If ( (time_se-time_sa)  >= 5) Then
        !--
        !  the internal procedure xcorrsav, adds the current data
        !  to that in the disk (if present). Note: in the process
        !  the variables xcorravgs, ntrajproc is overwritten and 
        !  is reset to zero in the end.  Since xcorrsav involves
        !  consolidation of data across processors, no specific
        !  processor condition is used 
        !--

        !-- 
        !  no of trajdone cud vary across procs for the same time interval,
        !  a simple multiple of Nprocs will not hold good in general
        !  Note that since ntrajproc is transferred, it will be used
        !  as a buffer in the following calls; it nolonger retains
        !  the meaning of its name!
        !--
        !Call MPI_REDUCE (ntrajproc, gntrthis, 1,      &
        !     MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
        gntrthis = ntrajproc

        Call xcorrsav

        Call PPcorrsav


        time_sa  = time_se
     End If !-- to save or not to save

  End Do Trajectories
!!!*** ====================== Trajectories End ======================== ***!!!

  !** Final consolidation
  !Call MPI_REDUCE (ntrajproc, gntrthis, 1,      &
  !     MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
  gntrthis = ntrajproc

  Call xcorrsav
  Call PPcorrsav


  If (MyID == 0) Then
     Call CPU_Time(time_end)
     Print *, 'CPU_time usage:',  (time_end-time_begin)/3600.0, 'hrs'
  End If


  !--debug print *, Myid, iblock
  !Call MPI_Finalize(mpierr)      


!!!*** ======================= Format Statements ====================== ***!!!


Contains 

!!!*** ============ Internal Procedures Begin ======================= ***!!

  !!*** Subroutine GetaSeed ***!!
  Subroutine GetaSeed(aseed)
    Integer (k4b) aseed
    Character (len=12) clk(3)
    Integer clok(8)

    !-- clock(8) gives the milliseconds
    Call Date_and_time(clk(1), clk(2), clk(3), clok) 
    aseed = clok(8)*100000+clok(7)*1000+clok(6)*10+clok(5)

    !-- an unique seed incase clok is the same
    aseed = (aseed + 201271)*(MyId +1)  

    !-- shuffling, lifted from ran_1
    aseed = Ieor(aseed, Ishft(aseed,13))
    aseed = Ieor(aseed, Ishft(aseed,-17))
    aseed = Ieor(aseed, Ishft(aseed,5))

  End Subroutine GetaSeed



  !  Subroutine GMPI_Redsum(var,buf,N,mpierr)
  !    !--
  !    !  This procedure is a generic wrapper for MPI_Reduce(),
  !    !  for various floating point precisions
  !    !--
  !    Use Bspglocons
  !    Implicit None
  !
  !    Integer N,mpierr 
  !    Real (FPREC):: var(N), buf(N)
  !
  !    Include 'mpif.h'
  !
  !    !--debug print *, 'count requested = ', N, '. size(var) = ', size(var)
  !    if (N /= size(var)) Then
  !       print *, 'Mismatch of requested size and buffer size', N, size(var)
  !       return
  !    end if
  !    
  !
  !    Select Case(FPREC)
  !    Case(SNGL)
  !       Call MPI_Reduce (var, buf, N,      &
  !            MPI_Real, MPI_Sum, 0, MPI_Comm_World, mpierr)
  !    Case(DOBL)
  !       Call MPI_Reduce (var, buf, N,      &
  !            MPI_Double_Precision, MPI_Sum, 0, MPI_Comm_World, mpierr)
  !    Case DEFAULT
  !       Write (*,*) 'Wrong precision in GMPI_redsum()'
  !       Call MPI_Abort(MPI_Comm_World,mpierr)
  !    End Select
  ! End Subroutine GMPI_Redsum



  !!*** Subroutine xcorrsav ***!!
  Subroutine xcorrsav
    Integer count

    !-- -----------
    !  This internal procedure, adds the current static data to that 
    !  in the disk (if present), thru the following steps
    !  * consolidate current results from processors 
    !       xcorravgs -> gxcavgs
    !  * read data from disk into xcorravgs
    !  * add to gxcavgs
    !  * write gxcavgs to disk
    !  * reset xcorravgs to zero
    !
    !  Note on usage:  
    !    Since a consolidation step is involved, this procedure must
    !    be called by all processors, and not exclusive to a particular MyID.
    !-- -----------

    !! *** consolidate results from the other processors *** !!

    count = Size(xcorravgs)

    !Call GMPI_Redsum(xcorravgs,gxcavgs,count,mpierr)
    gxcavgs = xcorravgs


    If (MyID == 0) Then

       !-- obtain global values from the disk
       Inquire (file=gavgfile, exist=Filexists)
       If (Filexists) Then
          !-- obtain a lock (new or wait for one to be cleared) before 
          !-- performing i/o
          call getlock(gavgfile)
          !^^ to be unlocked only after new data is written

          Open (unit=gavunit,file=gavgfile,status='old')
          !-- 
          !  read in ntrajproc,xcorravgs as the current run's 
          !  values have already been transferred to global values
          !--
          Read (gavunit,*) ntrajproc
          Read (gavunit,*) xcorravgs
          Close (gavunit)
       Else ! if Globalavgs File_does_not_exist
          !-- obtain a lock 
          Call getlock(gavgfile) 
          !^^ tobe unlocked after new data is written

          ntrajproc = 0
          xcorravgs = 0
       End If

       !-- Note that the 2nd term in RHS is obtained from the disk
       gxcavgs = gxcavgs + xcorravgs
       gntrajdone = gntrthis + ntrajproc
       !-- gntrthis was consolidated in the main program, b4 this call


       !--
       !  Save the updated global averages data on to the disk
       !  Note that the gavunit file is still locked till this
       !  is completed
       !--

       Open(unit=gavunit,file=gavgfile,status='replace')

       Write (gavunit,*) gntrajdone
       Write (gavunit,*) gxcavgs

       Close(gavunit)

       !-- Release global avgs file for use to other instances of this program
       Call unlockfile(gavgfile)


       !-- Compute averages and cross correlations
       gxcavgs = gxcavgs/gntrajdone

       Open (unit = xcunit, file = xcfile, status="unknown")
       Write (xcunit, 590) '# Cross correlation functions separated &
            by one t1zimm.'
       Write (xcunit, 592) '# Trajectories: ', gntrajdone
       Write (xcunit, Format600) "# Time     ", (Prop_names(i), i = 1,NAvgs-1) 

       Do ni = 1,Nitot
          Write (xcunit, Format601) ni, &
               ( (gxcavgs(i,ni,3) -  &
               gxcavgs(i,ni,1) * gxcavgs(i,ni,2)), &
               i = 1,Navgs-1)
       End Do
       Close(xcunit)


       Open (unit = avunit, file = "avgs.dat", status="unknown")
       Write (avunit, 640) '#  Ensemble averages of functions &
            separated by one t1zimm, for the FHI and PREAV.' 
       Write (avunit,642) '# Trajectories:', gntrajdone

       Write (avunit, Format650 ) "# Time       ", (Prop_names(i), Prop_names(i) &
            , i = 1,NAvgs-1) 

       Do ni = 1,Nitot
          Write (avunit, Format651) ni, &
               (  gxcavgs(i,ni,1) , gxcavgs(i,ni,2) , &
               i = 1,Navgs-1)
       End Do
       Close(avunit)

    End If ! Root processor


    !-- 
    !  Since the data has been written to disk the following
    !  variables need to be reset before continuing with the simulation
    !--
    xcorravgs = 0
    ntrajproc = 0



590 Format (A)
592 Format (A,G12.5)
!600 Format (A, 2x, <NAvgs-1>(A13,2x))
!601 Format (G13.6, 2x, <NAvgs-1>(G13.6, 2x) )

640 Format (A)
642 Format (A,G12.5)
!650 Format (A, 2x, <2*(NAvgs-1)>(A13,2x))
!651 Format (G13.6, 2x, <2*(NAvgs-1)>(G13.6, 2x) )

  End Subroutine xcorrsav


  !!*** Subroutine PPcorrsav ***!!
  Subroutine PPcorrsav

    !-- -----------
    !  This internal procedure, adds the current Time correl data to that 
    !  in the disk (if present), thru the following steps
    !  * read data from disk into ppbuf
    !  * add to gppavg
    !  * write gppavg to disk
    !  * write time correl to disk
    !  * reset gppavg to zero
    !
    !-- -----------

    !-- 
    !  Note that in the case of PPcorr, the results have been consolidated
    !  from each processor after each trajectory
    !-- 



    If (MyID ==0) Then

       !-- obtain global values from the disk
       Inquire (file=gppfile, exist=Filexists)
       If (Filexists) Then
          !-- obtain a lock 
          Call getlock(gppfile) 

          Open (unit=gppunit,file=gppfile,status='old')

          !-- 
          !  read in ntrajproc, ppbuf as the current run's 
          !  values have already been transferred to global avgs
          !--

          Read (gppunit,*) ntrajproc
          gntrajdone = gntrthis + ntrajproc
          !-- gntrthis was consolidated in the main program, b4 this call

          !-- 
          !  Given that the total length and interval of
          !  sampling remains same for main and cv, and
          !  across several runs, the maximum sampled
          !  points in a block Nto(ib) (see sample_tcorr),
          !  remains same.  So it does not matter when it is
          !  read into the same variable.  However, not reading
          !  this will not print out the time correlation function
          !  in corfile for the case: when no trajectories
          !  are sampled in this run and there is only data in
          !  the disk to be processed.
          !--
          Read (gppunit,*) Iblmc, Nto

          Read (gppunit,*) ppbuf
          Gppavg(:,:,:,1) = Gppavg(:,:,:,1) + ppbuf

          !-- for the CV
          Read (gppunit,*) ppbuf
          Gppavg(:,:,:,2) = Gppavg(:,:,:,2) + ppbuf

          Close (gppunit)
       Else ! if Globalavgs File_does_not_exist
          !-- lock the file b4 opening a new one for writing
          Call getlock(gppfile) 

          !-- These statements are redundant, but just there to
          !-- show the similarities with xcorrsav()
          !-- ntrajproc = 0
          !-- ppbuf = 0 
       End If ! gppfile exists

       !--
       !  Save the new global averages data on to the disk
       !  Note that the gppunit file is still locked till this
       !  is completed
       !--

       Open(unit=gppunit,file=gppfile,status='replace')

       Write (gppunit,*) gntrajdone
       Write (gppunit,*) Iblmc, Nto
       Write (gppunit,*) Gppavg(:,:,:,1)
       Write (gppunit,*) Gppavg(:,:,:,2)

       Close(gppunit)

       !-- Release global avgs file for use to other programs (jobs)
       Call unlockfile(gppfile)


       !** Calculate time correlation functions

       Gppavg = Gppavg/gntrajdone

       !--
       !  The stress correlation function comes from six components
       !  of the stress tensor, and their dot product is stored
       !  in Ncorrel=3 and Ncorrel=4.  Average them out here,
       !  and store in Ncorrel=3
       !--
       Gppavg(3,:,:,:) = (Gppavg(3,:,:,:) + Gppavg(3,:,:,:))/6

       Open (unit = corunit, file = corfile, status="unknown")
       Write (corunit,*) '# Time correlation functions for the main and CV'
       Write (corunit,*) '# Trajectories completed:', gntrajdone
       Write (corunit,Format500) "Time       ", "Time/t1zimm"  &
            , (correl_names(i), correl_names(i), i = 1,Ncorrel-1)

       ib=1
       Do j = 0, 1
          ctime = j*Dctime*(NTOPB**(ib-1))
          Write (corunit, Format501 ) ctime, ctime/t1zimm &
               , (Gppavg(i,ib,j,1), Gppavg(i,ib,j,2), i = 1,Ncorrel-1)
       End Do

       Do ib = 1, Min(MAXBLKCOR, Iblmc)
          Do j = 2, Min(Nto(ib), NTOPB)
             ctime = j*Dctime*(NTOPB**(ib-1))
             Write (corunit, Format501 ) ctime, ctime/t1zimm &
                  , (Gppavg(i,ib,j,1), Gppavg(i,ib,j,2), i = 1,Ncorrel-1)
          End Do
       End Do

       Close(corunit)

       !** Write continuation information to disk 

       contfile="continue"
       Open(unit=cunit,file=contfile,status='replace')

       If (Ntrtot - gntrajdone > 0) Then
          Write (cunit,*) (Ntrtot - gntrajdone) &
               , '  more trajectories to be completed of ', Ntrtot
          Close(cunit)
       Else 
          Close(cunit,status='delete')
       End If

       !-- Note that gppavg is defined only in the root proc
       Gppavg = 0

    End If ! Root processor


    !-- 
    !  Since the data has been written to disk the following
    !  variables need to be reset before continuing with the simulation
    !  Gppavg, howver needs to b reset only in the root processor.
    !  Setting it to zero in other procs (without allocation)
    !  leads to the process hanging
    !--
    ntrajproc = 0



    !500 Format ('#', 2(A13, 2x), <2*(Ncorrel-1)>(A13, 2x))
    !501 Format (2(G13.6, 2x), <2*(Ncorrel-1)>(G13.6, 2x) )

  End Subroutine PPcorrsav


!!!*** ============ Internal Procedures End ======================= ***!!


End Program Eqbm2Ensemble



! some conventions followed for comments

!-- To obtain an outline, use grep "\!\*" <filename>

!!!*** ====== Main divisions (Parts) ======== ***!!!

!!***   Section  ***!!
!** Subsection

!-- Comment for the line(s) below
!^^ Comment for the line(s) above

!--
!   Block comment paragraph 
!   spanning
!   several 
!   lines 
!--

