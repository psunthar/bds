!!! Time-stamp: <properties.f90 06:58, 15 Mar 2008 by P Sunthar>

!________________________________________________________________
!   utilities to extract properties from a given configuration
!________________________________________________________________

!!! $Log: properties.f90,v $
!!! Revision 1.2  2005/05/24 01:06:44  pxs565
!!! variable NTOPB
!!!
!!! Revision 1.1  2005/02/25 21:25:47  pxs565
!!! Initial revision
!!!



Subroutine eqchain_props(NBeads, Rbead, F, propsums)
  Use bspglocons
  Use bspglovars
  Implicit None
  Integer, Intent (in) :: NBeads
  Real (FPREC) , Intent (in), Dimension(:,:)  :: Rbead, F
  Real (FPREC) , Intent (inout),Dimension(:) ::propsums

  !---------------------------------------------------------------------
  !    This routine estimates the contributions of a single bead-spring     
  !    chain in the configuration specified by the vector R, to the
  !    equilibrium and non-equilibrium properties of the solution. 
  !    (Time correlations are evaluated by a separte procedure.)
  !    
  !    The idea behind using this routine is to have a code wherein
  !    the main program can be run using different property estimation 
  !    routines (such as this) without any serious changes to the
  !    main program itself. In addition, the user should have the 
  !    freedom to choose properties of interest.
  !
  !    In the present routine, the properties are returned in the
  !    vector props whose elements are given below:
  !
  !    propsums(1) ---  square of end-to-end vector
  !    propsums(2) ---  radius of gyration
  !    propsums(3) ---  Avg stretch
  !--------------------------------------------------------------------
  Integer i, j, k, mu, nu
  Real (FPREC)  R(Ndim,NBeads), Rc(Ndim),dist(Ndim)
  Real (FPREC) lsum





  R = RBead

  !     Calculation of square of end-to-end vector

  dist = R(:,Nbeads) - R(:,1)

  propsums(1) = propsums(1) + Sum(dist*dist)


  !     Calculation of shape tensor

  !       First, calculate position of centre of mass
  Rc = 0.0
  Do mu = 1, NBeads
     Rc = Rc + R(:,mu)
  End Do
  Rc = Rc/NBeads

  ! shift the origin to the current center of mass
  Forall (mu = 1: NBeads)
     R(:,mu) = R(:,mu) - Rc
  End Forall



  ! prop (2)
  ! radius of gyration
  lsum = 0
  Do mu = 1, NBeads
     lsum = lsum + Sum(R(:,mu)*R(:,mu))
  End Do

  propsums(2) =  propsums(2) +  lsum/NBeads  


  !-- Average stretch 
  lsum = 0 
  Do i = 1,Ndim
     lsum = lsum + Maxval(R(i,:)) - Minval(R(i,:))
  End Do
  propsums(3) = propsums(3) + lsum/Ndim


End Subroutine eqchain_props





!!$ Subroutine Sample_blocksums(action, sfile, ddisp)
!!$   Use Bspglocons
!!$   Use Bspglovars
!!$   Use Blocksums
!!$   Implicit None
!!$ 
!!$   Integer, Intent(in) :: action
!!$   Character (*),  Intent (in), Optional :: sfile
!!$ 
!!$   !-- displacement since last sample
!!$   Real (FPREC) , Intent(in), Optional :: ddisp(Ndim) 
!!$ 
!!$ 
!!$ !!!  This soubroutine takes samples for diffusivity using the
!!$ !!!  block-sum algorithm suggested in Frenkel Smit book on simulations
!!$ 
!!$   Integer, Save ::   tvacf                    ! time counter
!!$   Real (FPREC) , Save :: dsum(Ndim, MAXBLKDIF, NCOARSE)    ! displacement sums 
!!$ 
!!$ 
!!$ 
!!$   Integer, Parameter :: sunit = 40
!!$   Integer ii, ib, iblock, inmax, in, inp
!!$ 
!!$   Real (FPREC)  ctime, deld(Ndim)
!!$ 
!!$ !!! BLAS-1 functions
!!$   Real (SNGL)  sdot
!!$   Real (DOBL)  ddot
!!$ 
!!$   Select Case (action)
!!$ 
!!$   Case (INIT) ! Initialise variables and counters
!!$      tvacf = 0
!!$      Nbl = 0
!!$      Dsum = 0 
!!$      Nsqd = 0
!!$      Sqdisp = 0
!!$ 
!!$ 
!!$   Case (SAMPLE) ! sample
!!$ 
!!$ 
!!$      tvacf = tvacf + 1
!!$ 
!!$      !---determine current maximum number of blocks: Iblm
!!$      Iblm = 1
!!$      ii = tvacf/NCOARSE
!!$      Do While (ii.Ne.0)
!!$         Iblm = Iblm + 1
!!$         ii = ii/NCOARSE
!!$         !---test maximum time not longer than tdifmax:
!!$         ctime = (NCOARSE**(Iblm))
!!$         If (ctime > TDIFFMAX) ii = 0
!!$      End Do
!!$      !---limit the maximum number of blocks
!!$      Iblm = Min(Iblm, MAXBLKDIF)
!!$      Do ib = 1, Iblm
!!$         iblock = NCOARSE**(ib-1)
!!$         !---test for blocking operation
!!$         If (Mod(tvacf,iblock) == 0) Then
!!$            Nbl(ib) = Nbl(ib) + 1
!!$            !---limit to length n (=NCOARSE)
!!$            inmax = Min(Nbl(ib), NCOARSE)
!!$            !---for error calculation
!!$            If (ib == 1) Then
!!$               !---zero block: ordinary displacement
!!$               deld = ddisp
!!$            Else
!!$               !---(i)th block: coarsed displacement, previous block
!!$               deld = dsum(:,ib-1,1)
!!$            End If
!!$            Do in = 1, inmax
!!$               inp = in
!!$               If (Nbl(ib) > NCOARSE) inp = in + 1
!!$               If (in < inmax) Then
!!$                  dsum(:,ib,in) = dsum(:,ib,inp) + deld
!!$               Else
!!$                  dsum(:,ib,in) = deld
!!$               End If
!!$            End Do
!!$            Do in = 1, inmax
!!$               Nsqd(ib, in) = Nsqd(ib, in) + 1
!!$ 
!!$               Select Case (FPREC)
!!$               Case (SNGL)
!!$                  Sqdisp(ib, in) = Sqdisp(ib, in)  &
!!$                       +  sdot(Ndim,dsum(:,ib,inmax-in+1),1  &
!!$                       , dsum(:,ib,inmax-in+1),1)
!!$               Case (DOBL)
!!$                  Sqdisp(ib, in) = Sqdisp(ib, in)  &
!!$                       +  ddot(Ndim,dsum(:,ib,inmax-in+1),1  &
!!$                       , dsum(:,ib,inmax-in+1),1)
!!$               End Select
!!$            End Do
!!$         End If
!!$      End Do
!!$ 
!!$ 
!!$   Case (SAVDAT)
!!$ 
!!$      Open(unit=sunit, file = sfile, status = "unknown")
!!$      Write (sunit,*) tvacf
!!$      Write (sunit,*) Nbl
!!$      Write (sunit,*) Dsum
!!$      Write (sunit,*) Nsqd
!!$      Write (sunit,*) Sqdisp
!!$      Close(sunit)
!!$ 
!!$   Case (GETDAT)
!!$      Open(unit=sunit, file = sfile, status = "old")
!!$      Read (sunit,*) tvacf
!!$      Read (sunit,*) Nbl
!!$      Read (sunit,*) Dsum
!!$      Read (sunit,*) Nsqd
!!$      Read (sunit,*) Sqdisp
!!$      Close(sunit)
!!$ 
!!$ 
!!$   Case DEFAULT 
!!$      Write (*,*) 'Wrong Case called in Sample_blocksums'
!!$   End Select
!!$ 
!!$ End Subroutine Sample_blocksums



Subroutine Sample_Tcorr(action, sfile, pval)
  Use Bspglocons
  Use Bspglovars
  Use BlockTcorr
  !doing TODO
  Implicit None

  Integer, Intent(in) :: action

  Character (*),  Intent (in), Optional :: sfile

  ! current values of properties
  Real (FPREC) , Intent(in), Optional :: pval(Ndim,Ncorrel) 


!!!  This soubroutine calculates the time correlation function of
!!!  various properties, using the 
!!!  blocking algorithm suggested in Frenkel Smit book on simulations
!!!  Adapted from sample2.f, obtainable from their website

  Integer, Save ::   tacf                    ! time counter

  Real (FPREC) , Save, Allocatable :: p0val(:,:,:,:)  
  !values at time origins
  Real (FPREC)  :: curp(Ndim)


  Integer, Parameter :: sunit = 40
  Integer pp, ii, ib, iblock, inmin, inmax, in, inp

  Real (FPREC)  ctime, deld(Ndim)

  !-- for gfortran
  Character (len=80) :: Format9876


  !--
  !        Initialization variable format expressions
  !         <> language extension not available in gfortran
  !--


  !9876 Format (2(I3,1x),4x,<NTOPB+1>(2x,F8.0))
  Write(Format9876,"(a,I3,a)") "(2(I3,1x),4x, ",  NTOPB+1, "(2x,F8.0))"



  Select Case (action)

  Case (ALLOCAT) ! Initialise variables and counters
     call allocate_variables

  Case (INIT) ! Initialise variables and counters
     tacf = -1
     Nto = -1
     p0val = 0
     Npp = 0
     PPcorr = 0


  Case (SAMPLE) ! sample


     tacf = tacf + 1


     !---determine current maximum number of blocks: Iblmc
     Iblmc = 1
     ii = tacf/NTOPB
     Do While (ii.Ne.0)
        Iblmc = Iblmc + 1
        ii = ii/NTOPB
        !---test maximum time not longer than tdifmax:
        ctime = (NTOPB**(Iblmc))
        If (ctime > TDIFFMAX) ii = 0
     End Do
     !---limit the maximum number of blocks
     Iblmc = Min(Iblmc, MAXBLKCOR)

     Blocks: Do ib = 1, Iblmc
        iblock = NTOPB**(ib-1)
        !---test for blocking operation
        If (Mod(tacf,iblock) == 0) Then
           Nto(ib) = Nto(ib) + 1

           !-- leave out repeated calculation of <P0 P0> and <P0P1>
           inmin = 0
           If (ib > 1) inmin = 2
           !---limit to length n (=NTOPB)
           inmax = Min(Nto(ib), NTOPB)


           Properties: Do pp=1,Ncorrel
              curp = pval(:,pp)

              If (Nto(ib) > NTOPB)  Then
                 !-- Shift the time origins by one toward left
                 Do in = 0, NTOPB-1
                    p0val(:,pp,ib,in) = p0val(:,pp,ib,in+1)
                 End Do
              End If

              !--- save time origin values
              p0val(:,pp,ib,inmax) = curp



              !---- compute the correlation function  <P(t) P(0)>

              Do in = inmin, inmax
                 Npp(ib, in) = Npp(ib, in) + 1 
                 ! note increments by Ncorrel for each sample
                 PPcorr(pp, ib, in) = PPcorr(pp,ib, in)  &
                      +  Sum(p0val(:,pp,ib,inmax-in) * curp(:))
              End Do
           End Do Properties
        End If
     End Do Blocks

     If (tacf == 0) Then ! fill other blocks with the foremost velocity
        Do ib = 2, MAXBLKCOR
           Nto(ib) = Nto(ib) + 1
           Do pp = 1,Ncorrel
              p0val(:,pp,ib,0) = pval(:,pp)
           End Do
        End Do
     End If


  Case (PRNT)

     Do ib = 1, Min(MAXBLKCOR, Iblmc)
        If (Mod(tacf,NTOPB**(ib-1)) == 0) Then
           Write (10+ib,Format9876) tacf, Nto(ib), (p0val(1,1,ib,in), in=0,NTOPB)
           Write (20+ib,Format9876) tacf, Npp(ib,0), (PPcorr(1,ib,in), in=0,NTOPB)
        End If
     End Do
     !9876 Format (2(I3,1x),4x,<NTOPB+1>(2x,F8.0))

  Case (SAVDAT)

     Open(unit=sunit, file = sfile, status = "unknown")
     Write (sunit,*) tacf
     Write (sunit,*) Nto
     Write (sunit,*) p0val
     Write (sunit,*) Npp
     Write (sunit,*) PPCorr
     Close(sunit)

  Case (GETDAT)
     Open(unit=sunit, file = sfile, status = "old")
     Read (sunit,*) tacf
     Read (sunit,*) Nto
     Read (sunit,*) p0val
     Read (sunit,*) Npp
     Read (sunit,*) PPCorr

     Close(sunit)

  Case DEFAULT 
     Write (*,*) 'Wrong Case called in Sample_Tcorr'
  End Select


Contains
!!! =======================================================================
!!!                    Internal Procedures Begin
!!! =======================================================================

  Subroutine allocate_variables

    !-- Allocate variables after NTOPB has been calculated from input
    
    If( .not. Allocated(p0val) ) &
         Allocate(p0val(Ndim, Ncorrel, MAXBLKCOR, 0:NTOPB))
    
    If( .not. Allocated(PPcorr) ) &
         Allocate(PPcorr(Ncorrel,MAXBLKCOR,0:NTOPB))

    If( .not. Allocated(Npp) ) &
         Allocate(Npp(MAXBLKCOR,0:NTOPB))

  End Subroutine allocate_variables
  

End Subroutine Sample_Tcorr



