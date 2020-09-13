!!!*** ====================   Declarations ======================= ***!!!
  !!*** Variables  ***!!
!!$  !*** Additional debugging information
  !** Model parameters and variables
  !** Time-Integration  parameters
  !** Sampling parameters
  !** File i/o
  !** Trajectories related
  !** Misc
!!!*** =============== Execution begin ============================ ***!!!
  !  !!*** MPI Intialisation ***!!
  !!*** Input DATA ***!!
  !!*** Memory Allocation  ***!!
  !** Variables specific to root processor only
  !!*** Initialisation ***!!
  !** Generate a pseudo random seed
  !** Other initialisations
  !** Continuation run (if possible)
!!!*** =================== Trajectories Begin ======================== ***!!!
     !!***  Main Variate (with fluctuating HI) ***!!
     !** Equilibriate chain configuration 
     !**  Production run
     !!***  Repeat for the control variate (preaveraged HI) ***!!
     !** Equilibriate chain configuration 
     !** Production run
!!!*** ====================== Trajectories End ======================== ***!!!
  !** Final consolidation
!!!*** ======================= Format Statements ====================== ***!!!
!!!*** ============ Internal Procedures Begin ======================= ***!!
  !!*** Subroutine GetaSeed ***!!
  !!*** Subroutine xcorrsav ***!!
  !!*** Subroutine PPcorrsav ***!!
       !** Calculate time correlation functions
       !** Write continuation information to disk 
!!!*** ============ Internal Procedures End ======================= ***!!
!!!*** ====== Main divisions (Parts) ======== ***!!!
!!***   Section  ***!!
!** Subsection
