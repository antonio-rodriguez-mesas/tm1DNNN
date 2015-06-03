!********************************************************************
!
! TM1DNNN - Transfer matrix method for the Anderson model
! in 1D with finite-range hopping
!
!********************************************************************

!********************************************************************
!
! $Header: /home/cvs/phsht/GrapheneTMM/src/Restart/ALL/mainGR.f90,v 1.1 2011/07/22 17:49:19 ccspam Exp $
!
!********************************************************************

!********************************************************************
!
!      Comments:
!
!      IBCFlag         0        hardwall or
!                      1        periodic or
!                      2        antiperiodic boundary conditions
!
!      IKeepFlag       0        overwrite old .raw files
!                      1        keep old .raw files and skip to next
!                               configuration
!
!      IWriteFlag      0        no protocolling
!                      1        write out final wave function for each
!                               strip width and flux (disorder/U)
!                      2        additionally, write out each nGamma
!                               value after reorthonormalisation
!
!      ISortFlag       0        NO sorting of eigenvalues/vectors after
!                               each reorthonormalization
!                      1        YES, sorting is done.
!
!      IFluxFlag       0        DiagDis loop
!                      1        Energy loop
!
!      Notation:
!
!      see B. Kramer and M. Schreiber, "Transfer-Matrix Methods and Finite-
!	Size Scaling for Disordered Systems" and K. Frahm et al., EPL 31,
!      169 (1995).
!
!      Also, Martin Hennecke, "Anderson Uebergang im Magnetfeld", PTB
!      Bericht, PTB-PG-6, Braunschweig, Dec 94.
!
!********************************************************************

PROGRAM TM1DNNN

  !--------------------------------------------------------------------
  ! parameter and global variable definitions
  !--------------------------------------------------------------------

  USE MyNumbers

  USE CConstants
  USE IConstants
  USE IChannels

  USE IPara
  USE DPara

  USE Randoms
  USE Extensions
  USE Gammas

  USE RNG

  !--------------------------------------------------------------------
  ! local variable definitions
  !--------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER IWidthSquared, IWidthRL, Width0_X, &
       IChannelMax, index,jndex,kndex,lndex, jjndex,kkndex, andex,bndex, &
       Iter2, NOfG,iG, Ilayer

  REAL(KIND=RKIND) flux, flux0,flux1,dflux
  REAL(KIND=RKIND) KappaA, KappaB, sum

  REAL(KIND=RKIND), DIMENSION(:), ALLOCATABLE ::             &
       !nGamma, gamma, gamma2, acc_variance,
       OnsitePotVec

  REAL(KIND=RKIND), DIMENSION(:,:), ALLOCATABLE ::             &
       HopMatiLR, HopMatiL, EMat, dummyMat

  !COMPLEX(KIND=CKIND), DIMENSION(:,:), ALLOCATABLE :: PsiA, PsiB

  CHARACTER*40 surname
  CHARACTER*3 EnStr

  ! timing variables
  INTEGER :: IHours,IMinutes,ISeconds,IMilliSeconds, &
       IStartTime, ICurrentTime ,IRate

  REAL(RKIND) :: Duration, time

  !REAL(KIND=RKIND) dummy, dummyB
  INTEGER IErr

  !General Parameters
  LOGICAL :: LExist
  INTEGER(KIND=IKIND) :: TMM_CONVERGED
  REAL(KIND=RKIND) :: start, finish, hours

  !-------------------------------------------------------------------
  ! constants
  !-------------------------------------------------------------------

  CALL Init_Numbers

  !--------------------------------------------------------------------
  ! protocal feature startup
  !--------------------------------------------------------------------

  PRINT*,"TM1DNNN ASYMP (double precision)"
  PRINT*,RStr,DStr,AStr
  PRINT*,"--------------------------------------------------------------"

  !--------------------------------------------------------------------
  ! DBG mode
  !--------------------------------------------------------------------

  !MS$IF DEFINED (DBG)
  !PRINT*,"DBG mode"
  !MS$ENDIF

  !--------------------------------------------------------------------
  ! timing startup
  !--------------------------------------------------------------------

  CALL SYSTEM_CLOCK(count_rate=IRate)
  CALL SYSTEM_CLOCK(IStarttime)

  !--------------------------------------------------------------------
  ! input handling
  !--------------------------------------------------------------------

  CALL Input( IErr )
  !PRINT*, "DBG: IErr=", IErr
  IF( IErr.NE.0 ) THEN
     PRINT*,"main: error in Input()"
     STOP
  ENDIF

  IF (ISortFlag.EQ.1) THEN
     PRINT*,"ReSort() is done after each renormalization !"
     PRINT*,"-------------------------------------------------------"
  ENDIF

  IRestart = IRestart - 1

  !-----------------------------------------------------------------
  ! choose between square or ZZ
  !-----------------------------------------------------------------

  IF(IKappaFlag.EQ.0) THEN
     KappaA     = Kappa
     KappaB     = Kappa
  ELSE
     KappaA     = Kappa
     KappaB     = 0.0D0
  ENDIF

  !--------------------------------------------------------------------
  ! select the quantity to be varied
  !--------------------------------------------------------------------

  SELECT CASE(IFluxFlag)
  CASE(0)
     PRINT*,"main: Varying Diagonal Disorder" 
     flux0     = DiagDis0
     flux1     = DiagDis1
     dflux     = dDiagDis
  CASE(1)
     PRINT*, "main: Varying Energy"
     flux0     = Energy0
     flux1     = Energy1
     dflux     = dEnergy
  END SELECT

  !--------------------------------------------------------------------
  ! main parameter sweep
  !--------------------------------------------------------------------

  range_loop: &
  DO IRange= IRange0,IRange1,dIRange

     !--------------------------------------------------------------
     ! the # Lyapunov exponents is maximally .EQ. to IRange
     !--------------------------------------------------------------
     
     NOfG= MIN( NOfGamma, IRange )
     
     !--------------------------------------------------------------------
     ! Set-Up AVG Files
     !--------------------------------------------------------------------
     
     SELECT CASE(IKeepFlag)
     CASE(0)
        CALL OpenOutputAvg(IErr)
        IF( IErr.EQ.1 ) THEN
           PRINT*,"main: Error in OpenOutputAvg()"
           STOP
        ENDIF
     CASE(1)
        CALL CheckOutputAvg(IErr,LExist)
        IF( IErr.EQ.1 ) THEN
           PRINT*,"main: Error in CheckOutputAvg()"
           STOP
        ENDIF
        SELECT CASE(LExist)
        CASE(.TRUE.)
           CALL ReOpenOutputAvg(IErr)
           IF( IErr.EQ.1 ) THEN
              PRINT*,"main: Error in ReOpenOutputAvg()"
              STOP
           ENDIF
        CASE(.FALSE.)
           CALL OpenOutputAvg(IErr)
           IF( IErr.EQ.1 ) THEN
              PRINT*,"main: error in OpenOutputAvg()"
              STOP
           ENDIF
        END SELECT
     END SELECT
     
     !--------------------------------------------------------------
     ! allocate the memory for the TMM, gamma and RND vectors
     !--------------------------------------------------------------
     
     ALLOCATE(PsiA(IRange,IRange), STAT = IErr)
     ALLOCATE(PsiB(IRange,IRange), STAT = IErr)
     
     ALLOCATE(nGamma(IRange), STAT = IErr)
     ALLOCATE(gamma(IRange), STAT = IErr)
     ALLOCATE(gamma2(IRange), STAT = IErr)
     ALLOCATE(acc_variance(IRange), STAT = IErr)

     ALLOCATE(HopMatiLR(IRange,IRange), STAT = IErr)
     ALLOCATE(HopMatiL(IRange,IRange), STAT = IErr)
     ALLOCATE(EMat(IRange,IRange), STAT = IErr)
     ALLOCATE(dummyMat(IRange,IRange), STAT = IErr)
     ALLOCATE(OnsitePotVec(IRange), STAT = IErr)
     
     !PRINT*, "DBG: IErr=", IErr
     IF( IErr.NE.0 ) THEN
        PRINT*,"main: error in ALLOCATE()"
        STOP
     ENDIF
     
     !--------------------------------------------------------------
     ! flux loop
     !--------------------------------------------------------------
     
flux_loop: &
     DO flux= flux0,flux1,dflux

        !--------------------------------------------------------------
        ! set values for the physical quantities
        !--------------------------------------------------------------
  
        SELECT CASE(IFluxFlag) 
        CASE(0)
           DiagDis = flux
           Energy  = Energy0
        CASE(1)
           Energy  = flux
           DiagDis = DiagDis0
        END SELECT

        !--------------------------------------------------------------
        ! save the current parameters
        !--------------------------------------------------------------

        CALL SaveCurrent(IErr)
        IF(IErr.NE.0) THEN
           PRINT*,"main: ERR in SaveCurrent()"
        ENDIF

        !--------------------------------------------------------------
        ! open the nGamma file
        !--------------------------------------------------------------

        IF(IWriteFlag.GE.2) THEN
           CALL OpenOutputGamma( IRange, DiagDis,Energy, IErr )
           IF( IErr.NE.0 ) THEN
              PRINT*,"main: error in OpenOutputGamma()"
              STOP
           ENDIF
        ENDIF

        !--------------------------------------------------------------
        ! protocoll feature
        !--------------------------------------------------------------

2600    WRITE(*,2610) IRange, DiagDis, Energy, Kappa
2610    FORMAT("START @ IW= ",I4.1,    &
             ", DD= ", G10.3,          &
             ", En= ", G10.3,          &
             ", KA= ", G10.3)

        !--------------------------------------------------------------
        ! Initialize the convergence flags
        !--------------------------------------------------------------

        TMM_CONVERGED = 0

        !--------------------------------------------------------------
        ! initialize the wave vectors, the gamma sums and RNG
        !--------------------------------------------------------------

        IF(IMidDump .EQ. 1 .AND. IKeepFlag .EQ. 1) THEN
           CALL LoadLoopParams(IErr)
           IKeepFlag = 0
           IF( IErr.NE.0 ) THEN
              PRINT*,"main: error in LoadLoopParams()"
              STOP
           ENDIF

        ELSE

           ! start the index
           index1= 1

           ! reset the wave vectors
           PsiA= ZERO
           PsiB= ZERO      

           !reset the gamma sum approximants for a SINGLE
           !configuration
           DO iG=1,IRange
              gamma(iG)          = ZERO
              gamma2(iG)         = ZERO
              acc_variance(iG)   = ZERO            
           ENDDO

           DO index=1,IRange
              PsiA(index,index)= ONE
           ENDDO

           !PRINT *, 'DBG: PsiA=', PsiA
           !PRINT *, 'DBG: PsiB=', PsiB

           !initialize the RNG
           CALL SETVALUES(ISeed)
        ENDIF

        !--------------------------------------------------------------
        ! initialize the connectivity and energy matrices
        !--------------------------------------------------------------

        HopMatiL=ZERO
        HopMatiLR=ZERO

        DO index=1,IRange
           DO jndex=1,IRange
              IF( index .LE. jndex) THEN
                 HopMatiL(index,jndex)= REAL(IRange - jndex + index,RKIND)
              ENDIF

              IF( index .GE. jndex) THEN
                 HopMatiLR(index,jndex)= -REAL(IRange - index + jndex,RKIND)
              ENDIF

              ! only off-diagonal elements of EMat
              IF( index .NE. jndex) THEN
                 EMat(index,jndex)= REAL(-ABS(index-jndex),RKIND)
              ENDIF
           ENDDO
        ENDDO

!!$        PRINT*,"HopMatL=", HopMatiL
!!$        PRINT*,"HopMatR=", HopMatiLR
!!$        PRINT*,"EMat=", EMat
!!$        PAUSE

        CALL INVERT(IRange,HopMatiL,dummyMat,IErr)
        IF( IErr.NE.0) THEN
           PRINT*,"main(): IErr=", IErr," in INVERT() --- aborting!"
           STOP
        ENDIF
        HopMatiL= dummyMat

        HopMatiLR= MATMUL(HopMatiL,HopMatiLR)

!!$        PRINT*,"HopMatiL=", HopMatiL
!!$        PRINT*,"HopMatiLR=", HopMatiLR

        !--------------------------------------------------------------
        ! iteration loop
        !--------------------------------------------------------------

tmm_loop:&
        DO Iter1= index1, MAX( (NOfIter)/(NOfOrtho), 1)

           DO Iter2= 1, NOfOrtho, 2
   
              ! Ilayer is the current horizontal position
              Ilayer= (Iter1-1)*NOfOrtho+Iter2
              
              ! do the TM multiplication             
              CALL TMMultNNN( PsiA, PsiB, Ilayer, &
                   HopMatiLR, HopMatiL, EMat, dummyMat, OnsitePotVec, &
                   Energy, DiagDis, IRange) 

              CALL TMMultNNN( PsiB, PsiA, Ilayer+1, &
                   HopMatiLR, HopMatiL, EMat, dummyMat, OnsitePotVec, &
                   Energy, DiagDis, IRange)  

              !CALL Swap( PsiA, PsiB, IRange)

           ENDDO !northo_loop

           !-------------------------------------------------------
           ! renormalize via Gram-Schmidt
           !-------------------------------------------------------

           CALL ReNorm(PsiA,PsiB,gamma,gamma2,IRange)

           IF(IWriteFlag.GE.MAXWriteFLAG+1) THEN   
              PRINT*,"DBG: Orth,iL,PsiA", iLayer, PsiA
              PRINT*,"DBG: Orth,iL,PsiB", iLayer, PsiB; !PAUSE
           ENDIF

           !-------------------------------------------------------
           ! sort the eigenvalues by LARGEST first AND also
           ! resort the eigenvectors accordingly
           !-------------------------------------------------------           

           IF (ISortFlag.EQ.1) THEN
              CALL ReSort(PsiA,PsiB,gamma,gamma2,IRange)          
           ENDIF

           !--------------------------------------------------------
           ! "Iter1" counts the number of actual renormalizations
           ! of the transfer matrix.
           !--------------------------------------------------------

           !-----------------------------------------------------------
           ! do the gamma computations
           !-----------------------------------------------------------         

           DO iG=1, IRange

              nGamma(IRange+1-iG)= gamma(iG)/REAL(NOfOrtho*Iter1)

              acc_variance(IRange+1-iG)=             &
                   SQRT( ABS(                        &
                   (gamma2(iG)/REAL(Iter1) -         &
                   (gamma(iG)/REAL(Iter1))**2 )      &
                   / REAL( MAX(Iter1-1,1) )          &
                   )) / ABS( gamma(iG)/REAL(Iter1) )
           ENDDO

           !-----------------------------------------------------------
           ! write the nGamma data to file
           !-----------------------------------------------------------

           IF(IWriteFlag.GE.2) THEN
              CALL WriteOutputGamma( Iter1, & 
                   nGamma, acc_variance, & 
                   NOfG, IErr )
              !PRINT*, "DBG: IErr=", IErr
              IF( IErr.NE.0 ) THEN
                 PRINT*,"main: error in WriteOutputGamma()"
                 STOP
              ENDIF
           ENDIF

           !--------------------------------------------------------
           ! write the nGamma data to stdout
           !--------------------------------------------------------

           IF(IWriteFlag.GE.1 .AND. MOD(Iter1,NOfPrint).EQ.0 ) THEN
              WRITE(*,4210) Iter1, nGamma(1), acc_variance(1)
4210          FORMAT(I9.1, G15.7, G15.7)
           ENDIF

           !-----------------------------------------------------------
           ! check accuracy and dump the result
           !-----------------------------------------------------------

           IF(Iter1.GE.IRange .AND. &
                Iter1.GE.MINIter) THEN

              IF(acc_variance(1).LE.Epsilon .AND. &
                   acc_variance(1).GE.TINY) THEN
                 TMM_CONVERGED=1
                 GOTO 4000
              ENDIF

           ENDIF

           CALL cpu_time(finish)
           hours=(finish-start)/3600.0

           IF(hours .GT. IWalltime-1) THEN

              !PRINT*, "flag: saveloopparams"

              IF(IMidDump .EQ. 1) THEN
                 CALL SaveLoopParams(IErr)
                 IF( IErr.NE.0 ) THEN
                    PRINT*,"main: error in SaveLoopParams()"
                    STOP
                 ENDIF
              ENDIF

              !PRINT*, "flag: saveloopparams complete"
              !PRINT*, "flag: restartroutine"   

              SELECT CASE (IFluxFlag)
              CASE (0) !Disorder flux    

                 SELECT CASE (IBCFlag)

                 CASE (0) !Hardwall

                    IF(IRestart .GE. 1) THEN
                       CALL RestartRoutineZZDISHW(IErr)
                       IF( IErr.NE.0 ) THEN
                          PRINT*,"main: error in RestartRoutineZZDISHW()"
                          STOP
                       ENDIF
                    ENDIF

                 CASE (1) !Periodic

                    IF(IRestart .GE. 1) THEN
                       CALL RestartRoutineZZDISPBC(IErr)
                       IF( IErr.NE.0 ) THEN
                          PRINT*,"main: error in RestartRoutineZZDISPBC()"
                          STOP
                       ENDIF
                    ENDIF

                 END SELECT

              CASE (1) !Energy flux

                 SELECT CASE (IBCFlag)

                 CASE (0) !Hardwall

                    IF(IRestart .GE. 1) THEN
                       CALL RestartRoutineZZENEHW(IErr)
                       IF( IErr.NE.0 ) THEN
                          PRINT*,"main: error in RestartRoutineZZENEHW()"
                          STOP
                       ENDIF
                    ENDIF

                 CASE (1) !Periodic

                    IF(IRestart .GE. 1) THEN
                       CALL RestartRoutineZZENEPBC(IErr)
                       IF( IErr.NE.0 ) THEN
                          PRINT*,"main: error in RestartRoutineZZENEPBC()"
                          STOP
                       ENDIF
                    ENDIF

                 END SELECT

              END SELECT

              EXIT
           ENDIF

        ENDDO tmm_loop

        IF(TMM_CONVERGED .NE. 1) THEN
           PRINT*,"main: No convergence in TMM_loop!"
        ENDIF

4000    CONTINUE  

        !--------------------------------------------------------------
        ! write the AVG data
        !--------------------------------------------------------------

        CALL WriteOutputAvg( IRange, &
             DiagDis,Energy, &
             nGamma, acc_variance, &
             NOfG, PsiA, IErr, TMM_CONVERGED)
        !PRINT*, "DBG: IErr=", IErr
        IF( IErr.NE.0 ) THEN
           PRINT*,"main: error in WriteOutputAvg()"
           STOP
        ENDIF

        !--------------------------------------------------------
        ! write the nGamma data to stdout
        !--------------------------------------------------------

        WRITE(*,5010) Iter1, &
             DiagDis,Energy, Kappa
        WRITE(*,5012) nGamma(1), acc_variance(1)

5010    FORMAT("END @ ", I9.1, &
             ",", G15.7, ",", G15.7,",", G15.7)
5012    FORMAT("     ", &
             G15.7, ",", G15.7)

        !--------------------------------------------------------------
        ! close the nGamma file
        !--------------------------------------------------------------

        IF(IWriteFlag.GE.2) THEN
           CALL CloseOutputGamma( IErr )
           !PRINT*, "DBG: IErr=", IErr
           IF( IErr.NE.0 ) THEN
              PRINT*,"main: error in CloseOutputGamma()"
              STOP
           ENDIF
        ENDIF

        IF(hours .GT. IWalltime-1) THEN
           EXIT
        ENDIF

        !--------------------------------------------------------------
        ! end of flux loop
        !--------------------------------------------------------------

     ENDDO flux_loop

     !-----------------------------------------------------------------
     ! close the AVG files
     !-----------------------------------------------------------------
     
     CALL CloseOutputAvg( IErr )
     !PRINT*, "DBG: IErr=", IErr
     IF( IErr.NE.0 ) THEN
        PRINT*,"main: error in CloseOutputAvg()"
        STOP
     ENDIF

     !--------------------------------------------------------------
     ! delete the current parameters
     !--------------------------------------------------------------

     CALL DeleteCurrentParameters(IRange, IErr)
     IF(IErr.NE.0) THEN
        PRINT*,"main: ERR in DeleteCurrentParameters()"
     ENDIF

     !--------------------------------------------------------------
     ! DEallocate the memory for the TMM, gamma and RND vectors
     !--------------------------------------------------------------

     DEALLOCATE(PsiA,PsiB,nGamma,gamma,gamma2,acc_variance)
     DEALLOCATE(HopMatiLR,HopMatiL,EMat,dummyMat,OnsitePotVec)

     ! --------------------------------------------------
     ! get time at the end of the process
     ! --------------------------------------------------

     CALL SYSTEM_CLOCK(ICurrentTime)
     Duration=REAL(ICurrentTime-IStartTime)/REAL(IRate)
     IHours = FLOOR(Duration/3600.0D0)
     IMinutes = FLOOR(MOD(Duration,3600.0D0)/60.0D0)
     ISeconds = MOD(Duration,3600.0D0)-IMinutes*60
     IMilliSeconds = INT((Duration-(IHours*3600+IMinutes*60+ISeconds))*1000,IKIND)
     
     PRINT*, "tm1dNNN: used time=", IHours, "hrs ", &
          IMinutes,"mins ",ISeconds,"secs ", IMilliSeconds,"millisecs"

     !-----------------------------------------------------------------
     ! end of range loop
     !-----------------------------------------------------------------

  ENDDO range_loop

  STOP "TM1DNNN $Revision:$"

END PROGRAM TM1DNNN
