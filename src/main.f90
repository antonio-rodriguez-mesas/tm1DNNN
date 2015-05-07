!********************************************************************
!
! TMSEGR - Transfer matrix method for the Anderson
! model with diagonal disorder in Graphene
!
!********************************************************************

!********************************************************************
!
! $Header: /home/cvs/phsht/GrapheneTMM/src/Restart/ALL/mainGR.f90,v 1.1 2011/07/22 17:49:19 ccspam Exp $
!
!********************************************************************

!**************************************************************************
!
! $Log: mainGR.f90,v $
! Revision 1.1  2011/07/22 17:49:19  ccspam
! Programs for ZZ and AC with the restart routine for francesca
!
! Revision 1.1  2011/05/31 13:55:07  ccspam
! *** empty log message ***
!
! Revision 1.1  2011/05/06 08:13:09  phsht
! 1st installement

!**************************************************************************

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
!      IEdgeFlag       0        Zigzag or
!                      1        Armchair edge
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

PROGRAM TMSEGR

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
  
  !REAL(KIND=RKIND), DIMENSION(:), ALLOCATABLE ::             &
       !nGamma, gamma, gamma2, acc_variance
  
  !COMPLEX(KIND=CKIND), DIMENSION(:,:), ALLOCATABLE :: PsiA, PsiB

  CHARACTER*40 surname
  CHARACTER*3 EnStr
  
  ! timing variables
  REAL(4), DIMENSION(1:3,1:1) :: TIME, TIMESUM
  REAL(4) STARTTIME(2), ENDTIME(2) 
  REAL(4) T, T2, ETIME
  
  EXTERNAL ETIME
  
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
  
  PRINT*,"TMSEGR ASYMP (double precision)"
  !      PRINT*,"$Revision: 1.1 $ $Date: 2011/07/22 17:49:19 $ $Author: ccspam $"
  PRINT*,RStr,DStr,AStr
  !      PRINT*,"$Revision: 1.1 $ $Date: 2011/07/22 17:49:19 $ $Author: ccspam $"
  PRINT*,"--------------------------------------------------------------"
  
  !--------------------------------------------------------------------
  ! DBG mode
  !--------------------------------------------------------------------
  
  !MS$IF DEFINED (DBG)
  !PRINT*,"DBG mode"
  !MS$ENDIF
  
  !--------------------------------------------------------------------
  ! input handling and timing
  !--------------------------------------------------------------------
  
  CALL cpu_time(start)

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
  
  width_loop: &
  DO IWidth= Width0,Width1,dWidth

     ! --------------------------------------------------     
     ! get time at start of the process
     ! --------------------------------------------------
     T = ETIME(STARTTIME)
  
     !--------------------------------------------------------------
     ! the # Lyapunov exponents is maximally .EQ. to IWidth
     !--------------------------------------------------------------
     
     NOfG= MIN( NOfGamma, IWidth )

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
     
        ALLOCATE(PsiA(IWidth,IWidth), STAT = IErr)
        ALLOCATE(PsiB(IWidth,IWidth), STAT = IErr)
        
        ALLOCATE(nGamma(IWidth), STAT = IErr)
        ALLOCATE(gamma(IWidth), STAT = IErr)
        ALLOCATE(gamma2(IWidth), STAT = IErr)
        ALLOCATE(acc_variance(IWidth), STAT = IErr)
     
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
           CALL OpenOutputGamma( IWidth, DiagDis,Energy, IErr )
           IF( IErr.NE.0 ) THEN
              PRINT*,"main: error in OpenOutputGamma()"
              STOP
           ENDIF
        ENDIF

        !--------------------------------------------------------------
        ! protocoll feature
        !--------------------------------------------------------------

        SELECT CASE (IEdgeFlag)

        CASE(0)

2500    WRITE(*,2510) IWidth, DiagDis, Energy, KappaA, KappaB, MagFlux
2510    FORMAT("START @ IW= ",I4.1,    &
             ", DD= ", G10.3,          &
             ", En= ", G10.3,          &
             ", KA= ", G10.3,          &
             ", KB= ", G10.3,          &
             ", Mf= ", G10.3)
        
        CASE(1)

2600    WRITE(*,2610) IWidth, DiagDis, Energy, Kappa, MagFlux
2610    FORMAT("START @ IW= ",I4.1,    &
             ", DD= ", G10.3,          &
             ", En= ", G10.3,          &
             ", KA= ", G10.3,          &
             ", Mf= ", G10.3)
 
       END SELECT


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
          DO iG=1,IWidth
             gamma(iG)          = 0.0D0
             gamma2(iG)         = 0.0D0
             acc_variance(iG)   = 0.0D0            
          ENDDO
          
          DO index=1,IWidth
             PsiA(index,index)= ONE
          ENDDO
        
        !PRINT *, 'DBG: PsiA=', PsiA
        !PRINT *, 'DBG: PsiB=', PsiB

            !initialize the RNG
            CALL SETVALUES(ISeed)
        ENDIF
        
        !--------------------------------------------------------------
        ! iteration loop
        !--------------------------------------------------------------
        
tmm_loop: &
        DO Iter1= index1, MAX( (NOfIter)/(NOfOrtho), 1)

           SELECT CASE (IEdgeFlag)

           CASE(0) !ZigZag

           DO Iter2= 1, NOfOrtho, 2

              ! Ilayer is the current horizontal position
              Ilayer= (Iter1-1)*NOfOrtho+Iter2

              ! do the TM multiplication             
              CALL TMMult2DZZ( PsiA, PsiB, Ilayer, &
                    Energy, DiagDis, KappaA, KappaB, IWidth) 

              CALL TMMult2DZZ( PsiB, PsiA, Ilayer+1, &
                    Energy, DiagDis, KappaB, KappaA, IWidth)  
              
              !CALL Swap( PsiA, PsiB, IWidth)
              
           ENDDO !northo_loop

           CASE(1) !ArmChair

           DO Iter2= 1, NOfOrtho, 4

              ! Ilayer is the current horizontal position
              Ilayer= (Iter1-1)*NOfOrtho+Iter2

              SELECT CASE (IBCFlag)

              CASE (0) !Hardwall

                 ! do the TM multiplication             
                  CALL TMMultACHWI( PsiA, PsiB, Ilayer, &
                       Energy, DiagDis, IWidth) 
                  CALL TMMultACHWII( PsiA, PsiB, Ilayer+1, &
                       Energy, DiagDis, IWidth)
                  CALL TMMultACHWIII( PsiA, PsiB, Ilayer+2, &
                       Energy, DiagDis, IWidth) 
                  CALL TMMultACHWIV( PsiA, PsiB, Ilayer+3, &
                       Energy, DiagDis, IWidth)                
              
                 ! do the TM multiplication as a function of Kappa             
                 !CALL TMMultKACI( PsiA, PsiB, Ilayer, &
                   !   Energy, DiagDis, Kappa, IWidth) 
                 !CALL TMMultKACII( PsiA, PsiB, Ilayer+1, &
               !       Energy, DiagDis, Kappa, IWidth)
                 !CALL TMMultKACIII( PsiA, PsiB, Ilayer+2, &
               !       Energy, DiagDis, Kappa, IWidth) 
                 !CALL TMMultKACIV( PsiA, PsiB, Ilayer+3, &
               !       Energy, DiagDis, Kappa, IWidth)                

              CASE (1) !Periodic

              ! do the TM multiplication             
              CALL TMMult2DACI( PsiA, PsiB, Ilayer, &
                    Energy, DiagDis, IWidth) 
              CALL TMMult2DACII( PsiA, PsiB, Ilayer+1, &
                   Energy, DiagDis, IWidth)
              CALL TMMult2DACIII( PsiA, PsiB, Ilayer+2, &
                    Energy, DiagDis, IWidth) 
              CALL TMMult2DACIV( PsiA, PsiB, Ilayer+3, &
                   Energy, DiagDis, IWidth)                
              
              ! do the TM multiplication as a function of Kappa             
              !CALL TMMultKACI( PsiA, PsiB, Ilayer, &
               !     Energy, DiagDis, Kappa, IWidth) 
              !CALL TMMultKACII( PsiA, PsiB, Ilayer+1, &
               !     Energy, DiagDis, Kappa, IWidth)
              !CALL TMMultKACIII( PsiA, PsiB, Ilayer+2, &
               !     Energy, DiagDis, Kappa, IWidth) 
              !CALL TMMultKACIV( PsiA, PsiB, Ilayer+3, &
               !     Energy, DiagDis, Kappa, IWidth)    
              
               END SELECT

           ENDDO !northo_loop

           END SELECT
            
           !-------------------------------------------------------
           ! renormalize via Gram-Schmidt
           !-------------------------------------------------------
           
           CALL ReNorm(PsiA,PsiB,gamma,gamma2,IWidth)
           
           IF(IWriteFlag.GE.MAXWriteFLAG+1) THEN   
              PRINT*,"DBG: Orth,iL,PsiA", iLayer, PsiA
              PRINT*,"DBG: Orth,iL,PsiB", iLayer, PsiB; !PAUSE
           ENDIF

           !-------------------------------------------------------
           ! sort the eigenvalues by LARGEST first AND also
           ! resort the eigenvectors accordingly
           !-------------------------------------------------------           
           
           IF (ISortFlag.EQ.1) THEN
               CALL ReSort(PsiA,PsiB,gamma,gamma2,IWidth)          
           ENDIF
           
           !--------------------------------------------------------
           ! "Iter1" counts the number of actual renormalizations
           ! of the transfer matrix.
           !--------------------------------------------------------
           
           !-----------------------------------------------------------
           ! do the gamma computations
           !-----------------------------------------------------------         

              DO iG=1, IWidth
                
                 nGamma(IWidth+1-iG)= gamma(iG)/REAL(NOfOrtho*Iter1)
                
                 acc_variance(IWidth+1-iG)=             &
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
4210        FORMAT(I9.1, G15.7, G15.7)
           ENDIF
           
           !-----------------------------------------------------------
           ! check accuracy and dump the result
           !-----------------------------------------------------------
           
           IF(Iter1.GE.IWidth .AND. &
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

              SELECT CASE (IEdgeFlag)

              CASE (0) !ZigZag

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

             CASE (1) !ArmChair

                 SELECT CASE (IFluxFlag)

                 CASE (0) !Disorder flux    

                    SELECT CASE (IBCFlag)

                    CASE (0) !Hardwall
           
                    IF(IRestart .GE. 1) THEN
                      CALL RestartRoutineACDISHW(IErr)
                      IF( IErr.NE.0 ) THEN
                        PRINT*,"main: error in RestartRoutineACDISHW()"
                        STOP
                      ENDIF
                    ENDIF
                    
                    CASE (1) !Periodic
           
                    IF(IRestart .GE. 1) THEN
                      CALL RestartRoutineACDISPBC(IErr)
                      IF( IErr.NE.0 ) THEN
                        PRINT*,"main: error in RestartRoutineACDISPBC()"
                        STOP
                      ENDIF
                    ENDIF

                    END SELECT

                 CASE (1) !Energy flux

                    SELECT CASE (IBCFlag)

                    CASE (0) !Hardwall

                    IF(IRestart .GE. 1) THEN
                      CALL RestartRoutineACENEHW(IErr)
                      IF( IErr.NE.0 ) THEN
                         PRINT*,"main: error in RestartRoutineACENEHW()"
                         STOP
                      ENDIF
                    ENDIF
                   
                    CASE (1) !Periodic
           
                    IF(IRestart .GE. 1) THEN
                      CALL RestartRoutineACENEPBC(IErr)
                      IF( IErr.NE.0 ) THEN
                        PRINT*,"main: error in RestartRoutineACENEPBC()"
                        STOP
                      ENDIF
                    ENDIF
                
                    END SELECT

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
        
        CALL WriteOutputAvg( IWidth, &
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
     
     CALL DeleteCurrentParameters(IWidth, IErr)
     IF(IErr.NE.0) THEN
        PRINT*,"main: ERR in DeleteCurrentParameters()"
     ENDIF
     
     !--------------------------------------------------------------
     ! DEallocate the memory for the TMM, gamma and RND vectors
     !--------------------------------------------------------------
     
     DEALLOCATE(PsiA,PsiB,nGamma,gamma,gamma2,acc_variance)
     
     ! --------------------------------------------------
     ! get time at the end of the process
     ! --------------------------------------------------
     
     T2 = ETIME(ENDTIME) 
     
     TIME(1,1) = T2 -T
     TIME(2,1) = ENDTIME(1)-STARTTIME(1)
     TIME(3,1) = ENDTIME(2)-STARTTIME(2)
     
     !IF(IWriteFlag.GE.2) THEN
        T  = TIME(1,1)
        T2 = TIME(2,1)
        
        WRITE(*,'(A39,6F12.4)') &
             "SINGLE PROC --> TIME(USR,SYS,DIFF): ", &
             & TIME(1,1), TIME(2,1), TIME(3,1)
        
     !ENDIF
     
     !-----------------------------------------------------------------
     ! end of width loop
     !-----------------------------------------------------------------
     
  ENDDO width_loop
  
  STOP "TMSEGR $Revision: 1.1 $"
  
END PROGRAM TMSEGR
