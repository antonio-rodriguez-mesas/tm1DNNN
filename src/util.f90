! ********************************************************************
!       
! TM1DNNN - Transfer matrix method for the Anderson model
! in 1D with finite-range hopping
!
! ********************************************************************
       
! ********************************************************************
!     
! $Header: /home/cvs/phsht/GrapheneTMM/src/Restart/ALL/utilC.f90,v 1.1 2011/07/22 17:49:19 ccspam Exp $
!
! ********************************************************************

!--------------------------------------------------------------------
! TMMultNNN:
!
! Multiplication of the transfer matrix onto the vector (PSI_A,PSI_B), 
! giving (PSI_B,PSI_A) so that the structure of the transfer matrix 
! can be exploited
!-----------------------------------------------------------------------

SUBROUTINE TMMultNNN(PSI_A,PSI_B, Ilayer, En, DiagDis, KappaA, KappaB, M )

  USE MyNumbers
  USE IPara
  USE RNG
  USE DPara
  
  ! wave functions:
  !       
  ! (PSI_A, PSI_B) on input, (PSI_B,PSI_A) on output
  
  IMPLICIT NONE
  
  INTEGER Ilayer,           &! current # TM multiplications
       M                     ! strip width
  
  REAL(KIND=RKIND)  DiagDis,&! diagonal disorder
       En,                  &! energy
       KappaA,              &! inter-layer hopping odd
       KappaB                ! inter-layer hopping even

  
  COMPLEX(KIND=CKIND) PSI_A(M,M), PSI_B(M,M)
  
  INTEGER iSite, jState
  REAL(KIND=RKIND) OnsitePot
  COMPLEX(KIND=CKIND) new, PsiLeft, PsiRight

  ! first TM 
  
  DO iSite=1,M
     
     ! create the new onsite potential
     SELECT CASE(IRNGFlag)
     CASE(0)
        OnsitePot= En - DiagDis*(DRANDOM()-0.5D0)
     CASE(1)
        OnsitePot= En - DiagDis*(DRANDOM()-0.5D0)*SQRT(12.0D0)
     !CASE(2)
        !OnsitePot= En - GRANDOM(ISeedDummy,0.0D0,DiagDis)
     END SELECT
     
     DO jState=1,M

        ! Boundary conditions for the first TM
        
        IF (iSite.EQ.1) THEN
           
           IF (IBCFlag.EQ.0) THEN
              PsiLeft= CZERO            ! hard wall BC
           ELSE IF (IBCFlag.EQ.1) THEN
              PsiLeft= PSI_A(jState,M)  ! periodic BC
           ELSE IF (IBCFlag.EQ.2) THEN
              PsiLeft= -PSI_A(jState,M) ! antiperiodic BC
           ENDIF

        ELSE IF (iSite.EQ.M) THEN
              
           IF (IBCFlag.EQ.0) THEN
              PsiLeft= CZERO            ! hard wall BC
           ELSE IF (IBCFlag.EQ.1) THEN
              PsiLeft= PSI_A(jState,1)  ! periodic BC
           ELSE IF (IBCFlag.EQ.2) THEN
              PsiLeft= -PSI_A(jState,1) ! antiperiodic BC
           ENDIF

        ELSE 
           PsiLeft= PSI_A(jState,iSite+(-1)**(iSite)) 
        ENDIF

        PsiRight= PSI_A(jState,iSite-(-1)**(iSite))
 
        new= ( CMPLX(OnsitePot,0.0D0,CKIND) * PSI_A(jState,iSite) &
             - CMPLX(KappaA,0.0D0,CKIND) * PsiLeft &
             - CMPLX(KappaB,0.0D0,CKIND) * PsiRight &
             - PSI_B(jState,iSite) )
        
        PSI_B(jState,iSite)= new

     ENDDO ! jState
  ENDDO ! iSite
  
  !PRINT*,"PSIA(1,1),(2,1),(3,1),(4,1)",&
   !    PSI_A(1,1),PSI_A(2,1),PSI_A(3,1),PSI_A(4,1)
  !PRINT*,"PSIB(1,1),(2,1),(3,1),(4,1)",&
   !    PSI_B(1,1),PSI_B(2,1),PSI_B(3,1),PSI_B(4,1)
  
  RETURN
END SUBROUTINE TMMultNNN

!--------------------------------------------------------------------
! TMMult2DZZold:
!
! Multiplication of the transfer matrix onto the vector (PSI_A,PSI_B):
! 1TM: from (PSI_A_n,PSI_B_n-1) to (PSI_B_n,PSI_A)
! 2TM: from (PSI_A_n,PSI_B_n) to (PSI_A_n+1,PSI_B_n) 
! whith which the transfer matrix starts again
!----------------------------------------------------------------------

SUBROUTINE TMMult2DZZold(PSI_A,PSI_B, Ilayer, En, DiagDis, M )

  USE MyNumbers
  USE IPara
  USE RNG
  USE DPara
  
  ! wave functions:
  !       
  ! (PSI_A, PSI_B) on input, (PSI_B,PSI_A) on output
  
  IMPLICIT NONE
  
  INTEGER Ilayer,           &! current # TM multiplications
       M                     ! strip width
  
  REAL(KIND=RKIND)  DiagDis,&! diagonal disorder
       En                    ! energy
  
  COMPLEX(KIND=CKIND) PSI_A(M,M), PSI_B(M,M)
  
  INTEGER iSite, jState
  REAL(KIND=RKIND) OnsitePot
  COMPLEX(KIND=CKIND) newB, newA, PsiLeft, PsiRight

  ! first TM 
  
  DO iSite=1,M
     
     ! create the new onsite potential
     SELECT CASE(IRNGFlag)
     CASE(0)
        OnsitePot= En - DiagDis*(DRANDOM()-0.5D0)
     CASE(1)
        OnsitePot= En - DiagDis*(DRANDOM()-0.5D0)*SQRT(12.0D0)
     !CASE(2)
        !OnsitePot= En - GRANDOM(ISeedDummy,0.0D0,DiagDis)
     END SELECT
     
     !PRINT*,"iS,pL,RndVec", iSite,pLevel,RndVec((pLevel-1)*M+iSite)
     
     DO jState=1,M
        
        !PRINT*,"jState, iSite", jState, iSite,

        ! Boundary conditions for the first TM
        
        IF (iSite.EQ.1) THEN
           
           IF (IBCFlag.EQ.0) THEN
              PsiLeft= CZERO            ! hard wall BC
           ELSE IF (IBCFlag.EQ.1) THEN
              PsiLeft= PSI_A(jState,1)  ! periodic BC
           ELSE IF (IBCFlag.EQ.2) THEN
              PsiLeft= -PSI_A(jState,M) ! antiperiodic BC
           ENDIF

        ELSE

           IF (jState.EQ.M) THEN
              PsiLeft= PSI_A(jState,M)
           ELSE 
              PsiLeft= PSI_A(jState,iSite+(-1)**(iSite))
           ENDIF
 
       ENDIF
 
        newB= ( CMPLX(OnsitePot,0.0D0,CKIND) * PSI_A(jState,iSite) &
             - CMPLX(Kappa,0.0D0,CKIND) * PsiLeft &
             - PSI_B(jState,iSite) )
        
        PSI_B(jState,iSite)= newB

     ENDDO ! jState
  ENDDO ! iSite

  ! second TM
 
  DO iSite=1,M
     
     ! create the new onsite potential
     SELECT CASE(IRNGFlag)
     CASE(0)
        OnsitePot= En - DiagDis*(DRANDOM()-0.5D0)
     CASE(1)
        OnsitePot= En - DiagDis*(DRANDOM()-0.5D0)*SQRT(12.0D0)
     !CASE(2)
        !OnsitePot= En - GRANDOM(ISeedDummy,0.0D0,DiagDis)
     END SELECT
     
     DO jState=1,M

        PsiRight= PSI_B(jState,iSite-(-1)**(iSite))

        newA= ( CMPLX(OnsitePot,0.0D0,CKIND) * PSI_A(jState,iSite) &
             - CMPLX(Kappa,0.0D0,CKIND) * PsiRight &
             - PSI_B(jState,iSite) )

       PSI_A(jState,iSite)= newA
        
     ENDDO ! jState
  ENDDO ! iSite
  
  !PRINT*,"PSIA(1,1),(2,1),(1,2),(2,2)",&
       !PSI_A(1,1),PSI_A(2,1),PSI_A(1,2),PSI_A(2,2)
  
  RETURN
END SUBROUTINE TMMult2DZZold

!--------------------------------------------------------------------
! TMMult2DZZ:
!
! Multiplication of the transfer matrix onto the vector (PSI_A,PSI_B), 
! giving (PSI_B,PSI_A) so that the structure of the transfer matrix 
! can be exploited
!-----------------------------------------------------------------------

SUBROUTINE TMMult2DZZ(PSI_A,PSI_B, Ilayer, En, DiagDis, KappaA, KappaB, M )

  USE MyNumbers
  USE IPara
  USE RNG
  USE DPara
  
  ! wave functions:
  !       
  ! (PSI_A, PSI_B) on input, (PSI_B,PSI_A) on output
  
  IMPLICIT NONE
  
  INTEGER Ilayer,           &! current # TM multiplications
       M                     ! strip width
  
  REAL(KIND=RKIND)  DiagDis,&! diagonal disorder
       En,                  &! energy
       KappaA,              &! inter-layer hopping odd
       KappaB                ! inter-layer hopping even

  
  COMPLEX(KIND=CKIND) PSI_A(M,M), PSI_B(M,M)
  
  INTEGER iSite, jState
  REAL(KIND=RKIND) OnsitePot
  COMPLEX(KIND=CKIND) new, PsiLeft, PsiRight

  ! first TM 
  
  DO iSite=1,M
     
     ! create the new onsite potential
     SELECT CASE(IRNGFlag)
     CASE(0)
        OnsitePot= En - DiagDis*(DRANDOM()-0.5D0)
     CASE(1)
        OnsitePot= En - DiagDis*(DRANDOM()-0.5D0)*SQRT(12.0D0)
     !CASE(2)
        !OnsitePot= En - GRANDOM(ISeedDummy,0.0D0,DiagDis)
     END SELECT
     
     DO jState=1,M

        ! Boundary conditions for the first TM
        
        IF (iSite.EQ.1) THEN
           
           IF (IBCFlag.EQ.0) THEN
              PsiLeft= CZERO            ! hard wall BC
           ELSE IF (IBCFlag.EQ.1) THEN
              PsiLeft= PSI_A(jState,M)  ! periodic BC
           ELSE IF (IBCFlag.EQ.2) THEN
              PsiLeft= -PSI_A(jState,M) ! antiperiodic BC
           ENDIF

        ELSE IF (iSite.EQ.M) THEN
              
           IF (IBCFlag.EQ.0) THEN
              PsiLeft= CZERO            ! hard wall BC
           ELSE IF (IBCFlag.EQ.1) THEN
              PsiLeft= PSI_A(jState,1)  ! periodic BC
           ELSE IF (IBCFlag.EQ.2) THEN
              PsiLeft= -PSI_A(jState,1) ! antiperiodic BC
           ENDIF

        ELSE 
           PsiLeft= PSI_A(jState,iSite+(-1)**(iSite)) 
        ENDIF

        PsiRight= PSI_A(jState,iSite-(-1)**(iSite))
 
        new= ( CMPLX(OnsitePot,0.0D0,CKIND) * PSI_A(jState,iSite) &
             - CMPLX(KappaA,0.0D0,CKIND) * PsiLeft &
             - CMPLX(KappaB,0.0D0,CKIND) * PsiRight &
             - PSI_B(jState,iSite) )
        
        PSI_B(jState,iSite)= new

     ENDDO ! jState
  ENDDO ! iSite
  
  !PRINT*,"PSIA(1,1),(2,1),(3,1),(4,1)",&
   !    PSI_A(1,1),PSI_A(2,1),PSI_A(3,1),PSI_A(4,1)
  !PRINT*,"PSIB(1,1),(2,1),(3,1),(4,1)",&
   !    PSI_B(1,1),PSI_B(2,1),PSI_B(3,1),PSI_B(4,1)
  
  RETURN
END SUBROUTINE TMMult2DZZ

!--------------------------------------------------------------------
! TMMult2DACI:
!
! Multiplication of the transfer matrix onto the vector (PSI_A,PSI_B), 
! giving (PSI_A,PSI_B) so that the structure of the transfer matrix 
! can be exploited
!-----------------------------------------------------------------------

SUBROUTINE TMMult2DACI(PSI_A,PSI_B, Ilayer, En, DiagDis, M )

  USE MyNumbers
  USE IPara
  USE RNG
  USE DPara
  
  ! wave functions:
  !       
  ! (PSI_A, PSI_B) on input, (PSI_B,PSI_A) on output
  
  IMPLICIT NONE
  
  INTEGER Ilayer,           &! current # TM multiplications
       M                     ! strip width
  
  REAL(KIND=RKIND)  DiagDis, &! diagonal disorder
       En                  ! energy

  
  COMPLEX(KIND=CKIND) PSI_A(M,M), PSI_B(M,M)
  
  INTEGER iSite, jState, kIndex, lIndex, IErr, index
  REAL(KIND=RKIND) ax

  ! define a new matrix where store data

  COMPLEX(KIND=CKIND), DIMENSION(:,:), ALLOCATABLE :: PsiC
  COMPLEX(KIND=CKIND), DIMENSION(:,:), ALLOCATABLE :: PsiD
  REAL(KIND=RKIND), DIMENSION(:), ALLOCATABLE :: RndVec
  REAL(KIND=RKIND), DIMENSION(:), ALLOCATABLE :: OnsitePot


  ALLOCATE(PsiC(M,M), STAT = IErr)
  ALLOCATE(PsiD(M,M), STAT = IErr)
  ALLOCATE(RndVec(M), STAT = IErr)
  ALLOCATE(OnsitePot(M), STAT = IErr)

  ! reset to zero
  PsiC= ZERO
  PsiD= ZERO

  ! create the random numbers

  DO index=1,M
     ax=DRANDOM()
     RndVec(index)=ax
  ENDDO

  ! create the new onsite potential
 
  DO index=1,M
     SELECT CASE(IRNGFlag)
     CASE(0)
        OnsitePot(index)= En - DiagDis*(RndVec(index)-0.5D0)
     CASE(1)
        OnsitePot(index)= En - DiagDis*(RndVec(index)-0.5D0)*SQRT(12.0D0)
     !CASE(2)
        !OnsitePot(index)= En - GRANDOM(ISeedDummy,0.0D0,DiagDis)
     END SELECT
  ENDDO

  DO iSite=1,M
     
     DO jState=1,M

        ! First term of the equation is stored in PsiC and second one in PsiB
        
        IF (iSite.EQ.1) THEN
           
           DO kIndex=1,M

              PsiC(jState,iSite)= PsiC(jState,iSite) &
                                  + (OnsitePot(kIndex)*0.5D0) * (-(-1)**kIndex) &
                                  * PSI_A(jState,kIndex)

              PsiD(jState,iSite)= PsiD(jState,iSite) &
                                  + (0.5D0) * (-(-1)**kIndex) * PSI_B(jState,kIndex)

           ENDDO ! kIndex

        ELSE IF (MOD(iSite,2).EQ.1) THEN

             DO kIndex=iSite,M
        
                PsiC(jState,iSite)= PsiC(jState,iSite) &
                                   + (OnsitePot(kIndex)*0.5D0) * (-(-1)**kIndex) &
                                   * PSI_A(jState,kIndex)

                PsiD(jState,iSite)= PsiD(jState,iSite) &
                                    + (0.5D0) * (-(-1)**kIndex) * PSI_B(jState,kIndex)

             ENDDO !kIndex

             DO lIndex=1,iSite-1

                PsiC(jState,iSite)= PsiC(jState,iSite) &
                                    + (OnsitePot(lIndex)*0.5D0) * ((-1)**lIndex) &
                                    * PSI_A(jState,lIndex)

                PsiD(jState,iSite)= PsiD(jState,iSite) &
                                    + (0.5D0) * ((-1)**lIndex) * PSI_B(jState,lIndex)

             ENDDO !lIndex

             ELSE 
 
             DO kIndex=iSite,M
        
                PsiC(jState,iSite)= PsiC(jState,iSite) &
                                    + (OnsitePot(kIndex)*0.5D0) * ((-1)**kIndex) &
                                    * PSI_A(jState,kIndex)

                PsiD(jState,iSite)= PsiD(jState,iSite) &
                                    + (0.5D0) * ((-1)**kIndex) * PSI_B(jState,kIndex)

             ENDDO !kIndex

             DO lIndex=1,iSite-1

                PsiC(jState,iSite)= PsiC(jState,iSite) &
                                    + (OnsitePot(lIndex)*0.5D0) * (-(-1)**lIndex) &
                                    * PSI_A(jState,lIndex)

                PsiD(jState,iSite)= PsiD(jState,iSite) &
                                    + (0.5D0) * (-(-1)**lIndex) * PSI_B(jState,lIndex)

             ENDDO !lIndex
 
        ENDIF

     ENDDO ! jState
  ENDDO ! iSite

  ! Copy the result in PSI_B(M,M)

  DO iSite=1,M
     DO jState=1,M
        PSI_B(jState,iSite)=PsiC(jState,iSite)-PsiD(jState,iSite)
     ENDDO
  ENDDO

  DEALLOCATE(PsiC)
  DEALLOCATE(PsiD)
  DEALLOCATE(RndVec)
  DEALLOCATE(OnsitePot)
  
  !PRINT*,"PSIA(1,1),(1,2),(1,3)",&
       !PSI_A(1,1),PSI_A(1,2),PSI_A(1,3)
  !PRINT*,"PSIB(1,1),(1,2),(1,3)",&
       !PSI_B(1,1),PSI_b(1,2),PSI_B(1,3)
  
  RETURN
END SUBROUTINE TMMult2DACI

!--------------------------------------------------------------------
! TMMult2DACII:
!
! Multiplication of the transfer matrix onto the vector (PSI_A,PSI_B), 
! giving (PSI_A,PSI_B) so that the structure of the transfer matrix 
! can be exploited
!-----------------------------------------------------------------------

SUBROUTINE TMMult2DACII(PSI_A,PSI_B, Ilayer, En, DiagDis, M )

  USE MyNumbers
  USE IPara
  USE RNG
  USE DPara
  
  ! wave functions:
  !       
  ! (PSI_A, PSI_B) on input, (PSI_B,PSI_A) on output
  
  IMPLICIT NONE
  
  INTEGER Ilayer,           &! current # TM multiplications
       M                     ! strip width
  
  REAL(KIND=RKIND)  DiagDis, &! diagonal disorder
       En                  ! energy

  
  COMPLEX(KIND=CKIND) PSI_A(M,M), PSI_B(M,M)

  INTEGER iSite, jState, IErr
  REAL(KIND=RKIND) OnsitePot
  COMPLEX(KIND=CKIND) new, PsiLeft

  ! define a new matrix where store data

  COMPLEX(KIND=CKIND), DIMENSION(:,:), ALLOCATABLE :: PsiC

  ALLOCATE(PsiC(M,M), STAT = IErr)

  ! reset to zero
  PsiC= ZERO
  
  DO iSite=1,M
     
     ! create the new onsite potential
     SELECT CASE(IRNGFlag)
     CASE(0)
        OnsitePot= En - DiagDis*(DRANDOM()-0.5D0)
     CASE(1)
        OnsitePot= En - DiagDis*(DRANDOM()-0.5D0)*SQRT(12.0D0)
     !CASE(2)
        !OnsitePot= En - GRANDOM(ISeedDummy,0.0D0,DiagDis)
     END SELECT
     
     DO jState=1,M

        ! Evaluation of PsiLeft
        
        IF (iSite.EQ.1) THEN
           PsiLeft=PSI_A(jState,M)
        ELSE
           PsiLeft=PSI_A(jState,iSite-1)
        ENDIF
        
        new=( CMPLX(OnsitePot,0.0D0,CKIND) * PSI_B(jState,iSite) &
             - PsiLeft - PSI_A(jState,iSite) )

        PsiC(jState,iSite)=new

     ENDDO ! jState
  ENDDO ! iSite

  ! Copy the result in PSI_A(M,M)

  DO iSite=1,M
     DO jState=1,M
        PSI_A(jState,iSite)=PsiC(jState,iSite)
     ENDDO
  ENDDO

  DEALLOCATE(PsiC)
  
  !PRINT*,"PSIA(1,1),(1,2),(1,3)",&
       !PSI_A(1,1),PSI_A(1,2),PSI_A(1,3)
  !PRINT*,"PSIB(1,1),(1,2),(1,3)",&
       !PSI_B(1,1),PSI_b(1,2),PSI_B(1,3)
  
  RETURN
END SUBROUTINE TMMult2DACII

!--------------------------------------------------------------------
! TMMult2DACIII:
!
! Multiplication of the transfer matrix onto the vector (PSI_A,PSI_B), 
! giving (PSI_A,PSI_B) so that the structure of the transfer matrix 
! can be exploited
!-----------------------------------------------------------------------

SUBROUTINE TMMult2DACIII(PSI_A,PSI_B, Ilayer, En, DiagDis, M )

  USE MyNumbers
  USE IPara
  USE RNG
  USE DPara
  
  ! wave functions:
  !       
  ! (PSI_A, PSI_B) on input, (PSI_B,PSI_A) on output
  
  IMPLICIT NONE
  
  INTEGER Ilayer,           &! current # TM multiplications
       M                     ! strip width
  
  REAL(KIND=RKIND)  DiagDis, &! diagonal disorder
       En                 ! energy

  
  COMPLEX(KIND=CKIND) PSI_A(M,M), PSI_B(M,M)
  
  INTEGER iSite, jState, kIndex, lIndex, IErr, index
  REAL(KIND=RKIND) ax

  ! define a new matrix where store data

  COMPLEX(KIND=CKIND), DIMENSION(:,:), ALLOCATABLE :: PsiC
  COMPLEX(KIND=CKIND), DIMENSION(:,:), ALLOCATABLE :: PsiD
  REAL(KIND=RKIND), DIMENSION(:), ALLOCATABLE :: RndVec
  REAL(KIND=RKIND), DIMENSION(:), ALLOCATABLE :: OnsitePot


  ALLOCATE(PsiC(M,M), STAT = IErr)
  ALLOCATE(PsiD(M,M), STAT = IErr)
  ALLOCATE(RndVec(M), STAT = IErr)
  ALLOCATE(OnsitePot(M), STAT = IErr)

  ! reset to zero
  PsiC= ZERO
  PsiD= ZERO

  ! create the random numbers

  DO index=1,M
     ax=DRANDOM()
     RndVec(index)=ax
  ENDDO

  ! create the new onsite potential
 
  DO index=1,M
     SELECT CASE(IRNGFlag)
     CASE(0)
        OnsitePot(index)= En - DiagDis*(RndVec(index)-0.5D0)
     CASE(1)
        OnsitePot(index)= En - DiagDis*(RndVec(index)-0.5D0)*SQRT(12.0D0)
     !CASE(2)
        !OnsitePot(index)= En - GRANDOM(ISeedDummy,0.0D0,DiagDis)
     END SELECT
  ENDDO
  
  DO iSite=1,M
     
     DO jState=1,M

        ! First term of the equation is stored in PsiC and second one in PsiD
        
        IF (iSite.EQ.M) THEN
           
           DO kIndex=1,M

              PsiC(jState,iSite)= PsiC(jState,iSite) &
                                  + (OnsitePot(kIndex)*0.5D0) * (-(-1)**kIndex) &
                                  * PSI_A(jState,kIndex)

              PsiD(jState,iSite)= PsiD(jState,iSite) &
                                  + (0.5D0) * (-(-1)**kIndex) * PSI_B(jState,kIndex)

           ENDDO ! kIndex

        ELSE IF (MOD(iSite,2).EQ.1) THEN !iSite is odd

             DO kIndex=1,iSite
        
                PsiC(jState,iSite)= PsiC(jState,iSite) &
                                    + (OnsitePot(kIndex)*0.5D0) * (-(-1)**kIndex) &
                                    * PSI_A(jState,kIndex)

                PsiD(jState,iSite)= PsiD(jState,iSite) &
                                    + (0.5D0) * (-(-1)**kIndex) * PSI_B(jState,kIndex)

             ENDDO !kIndex

             DO lIndex=iSite+1,M

                PsiC(jState,iSite)= PsiC(jState,iSite) &
                                    + (OnsitePot(lIndex)*0.5D0) * ((-1)**lIndex) &
                                    * PSI_A(jState,lIndex)

                PsiD(jState,iSite)= PsiD(jState,iSite) &
                                    + (0.5D0) * ((-1)**lIndex) * PSI_B(jState,lIndex)

             ENDDO !lIndex

             ELSE !iSite is even
 
             DO kIndex=1,iSite
        
                PsiC(jState,iSite)= PsiC(jState,iSite) &
                                   + (OnsitePot(kIndex)*0.5D0) * ((-1)**kIndex) &
                                   * PSI_A(jState,kIndex)

                PsiD(jState,iSite)= PsiD(jState,iSite) &
                                    + (0.5D0) * ((-1)**kIndex) * PSI_B(jState,kIndex)

             ENDDO !kIndex

             DO lIndex=1+iSite,M

                PsiC(jState,iSite)= PsiC(jState,iSite) &
                                    + (OnsitePot(lIndex)*0.5D0) * (-(-1)**lIndex) &
                                    * PSI_A(jState,lIndex)

                PsiD(jState,iSite)= PsiD(jState,iSite) &
                                    + (0.5D0) * (-(-1)**lIndex) * PSI_B(jState,lIndex)

             ENDDO !lIndex
 
        ENDIF

     ENDDO ! jState
  ENDDO ! iSite

 ! Copy the result in PSI_B(M,M)

  DO iSite=1,M
     DO jState=1,M
        PSI_B(jState,iSite)=PsiC(jState,iSite)-PsiD(jState,iSite)
     ENDDO
  ENDDO

  DEALLOCATE(PsiC)
  DEALLOCATE(PsiD)
  DEALLOCATE(RndVec)
  DEALLOCATE(OnsitePot)
  
  !PRINT*,"PSIA(1,1),(1,2),(1,3)",&
       !PSI_A(1,1),PSI_A(1,2),PSI_A(1,3)
  !PRINT*,"PSIB(1,1),(1,2),(1,3)",&
       !PSI_B(1,1),PSI_b(1,2),PSI_B(1,3)
  
  RETURN
END SUBROUTINE TMMult2DACIII

!--------------------------------------------------------------------
! TMMult2DACIV:
!
! Multiplication of the transfer matrix onto the vector (PSI_A,PSI_B), 
! giving (PSI_A,PSI_B) so that the structure of the transfer matrix 
! can be exploited
!-----------------------------------------------------------------------

SUBROUTINE TMMult2DACIV(PSI_A,PSI_B, Ilayer, En, DiagDis, M )

  USE MyNumbers
  USE IPara
  USE RNG
  USE DPara
  
  ! wave functions:
  !       
  ! (PSI_A, PSI_B) on input, (PSI_B,PSI_A) on output
  
  IMPLICIT NONE
  
  INTEGER Ilayer,           &! current # TM multiplications
       M                     ! strip width
  
  REAL(KIND=RKIND)  DiagDis, &! diagonal disorder
       En                  ! energy

  
  COMPLEX(KIND=CKIND) PSI_A(M,M), PSI_B(M,M)

  INTEGER iSite, jState, IErr
  REAL(KIND=RKIND) OnsitePot
  COMPLEX(KIND=CKIND) new, PsiRight

  ! define a new matrix where store data

  COMPLEX(KIND=CKIND), DIMENSION(:,:), ALLOCATABLE :: PsiC

  ALLOCATE(PsiC(M,M), STAT = IErr)

  ! reset to zero
  PsiC= ZERO
  
  DO iSite=1,M
     
     ! create the new onsite potential
     SELECT CASE(IRNGFlag)
     CASE(0)
        OnsitePot= En - DiagDis*(DRANDOM()-0.5D0)
     CASE(1)
        OnsitePot= En - DiagDis*(DRANDOM()-0.5D0)*SQRT(12.0D0)
     !CASE(2)
        !OnsitePot= En - GRANDOM(ISeedDummy,0.0D0,DiagDis)
     END SELECT
     
     DO jState=1,M

        ! Evaluation of PsiLeft
        
        IF (iSite.EQ.M) THEN
           PsiRight=PSI_A(jState,1)
        ELSE
           PsiRight=PSI_A(jState,iSite+1)
        ENDIF
        
        new=( CMPLX(OnsitePot,0.0D0,CKIND) * PSI_B(jState,iSite) &
             - PsiRight - PSI_A(jState,iSite) )

        PsiC(jState,iSite)=new

     ENDDO ! jState
  ENDDO ! iSite

  ! Copy the result in PSI_A(M,M)

  DO iSite=1,M
     DO jState=1,M
        PSI_A(jState,iSite)=PsiC(jState,iSite)
     ENDDO
  ENDDO

  DEALLOCATE(PsiC)
  
  !PRINT*,"PSIA(1,1),(1,2),(1,3)",&
       !PSI_A(1,1),PSI_A(1,2),PSI_A(1,3)
  !PRINT*,"PSIB(1,1),(1,2),(1,3)",&
       !PSI_B(1,1),PSI_b(1,2),PSI_B(1,3)
  
  RETURN
END SUBROUTINE TMMult2DACIV

!--------------------------------------------------------------------
! TMMultKACI:
!
! Multiplication of the transfer matrix onto the vector (PSI_A,PSI_B), 
! giving (PSI_A,PSI_B) so that the structure of the transfer matrix 
! can be exploited
!-----------------------------------------------------------------------

SUBROUTINE TMMultKACI(PSI_A,PSI_B, Ilayer, En, DiagDis, Kappa, M )

  USE MyNumbers
  USE IPara
  USE RNG
  USE DPara
  
  ! wave functions:
  !       
  ! (PSI_A, PSI_B) on input, (PSI_B,PSI_A) on output
  
  IMPLICIT NONE
  
  INTEGER Ilayer,           &! current # TM multiplications
       M                     ! strip width
  
  REAL(KIND=RKIND)  DiagDis, &! diagonal disorder
       En,                 & ! energy
       Kappa                 ! diagonal hopping

  
  COMPLEX(KIND=CKIND) PSI_A(M,M), PSI_B(M,M)
  
  INTEGER iSite, jState, kIndex, lIndex, IErr
  REAL(KIND=RKIND) OnsitePot

  ! define a new matrix where store data

  COMPLEX(KIND=CKIND), DIMENSION(:,:), ALLOCATABLE :: PsiC
  COMPLEX(KIND=CKIND), DIMENSION(:,:), ALLOCATABLE :: PsiD

  ALLOCATE(PsiC(M,M), STAT = IErr)
  ALLOCATE(PsiD(M,M), STAT = IErr)

  ! reset to zero
  PsiC= ZERO
  PsiD= ZERO

  DO iSite=1,M
  
     ! create the new onsite potential
     SELECT CASE(IRNGFlag)
     CASE(0)
        OnsitePot= En - DiagDis*(DRANDOM()-0.5D0)
     CASE(1)
        OnsitePot= En - DiagDis*(DRANDOM()-0.5D0)*SQRT(12.0D0)
     !CASE(2)
        !OnsitePot= En - GRANDOM(ISeedDummy,0.0D0,DiagDis)
     END SELECT
     
     DO jState=2,M

        ! First term of the equation is stored in PsiC and second one in PsiB
        
        IF (iSite.EQ.1) THEN
           
           DO kIndex=1,M

              PsiC(jState,iSite)= PsiC(jState,iSite) &
                                  + (CMPLX(OnsitePot,0.0D0,CKIND)) * (1.0D0/(1.0D0+(Kappa)**M)) &
                                  * (-(-1)**kIndex) * (Kappa)**(kIndex-1) * PSI_A(jState,kIndex)

              PsiD(jState,iSite)= PsiD(jState,iSite) &
                                  + (1.0D0/(1.0D0+(Kappa)**M)) * (-(-1)**kIndex)  * (Kappa)**(kIndex-1) &
                                  * PSI_B(jState,kIndex)

           ENDDO ! kIndex

        ELSE IF (MOD(iSite,2).EQ.1) THEN

             DO kIndex=iSite+1,M
        
                PsiC(jState,iSite)= PsiC(jState,iSite) &
                                   + (CMPLX(OnsitePot,0.0D0,CKIND)) * (1.0D0/(1.0D0+(Kappa)**M)) &
                                   * (-(-1)**kIndex)  * ((Kappa)**(kIndex-iSite)) * PSI_A(jState,kIndex)

                PsiD(jState,iSite)= PsiD(jState,iSite) &
                                    + (1.0D0/(1.0D0+(Kappa)**M)) * (-(-1)**kIndex) * ((Kappa)**(kIndex-iSite)) &
                                    * PSI_B(jState,kIndex)

             ENDDO !kIndex

             DO lIndex=1,iSite-1

                PsiC(jState,iSite)= PsiC(jState,iSite) &
                                    + (CMPLX(OnsitePot,0.0D0,CKIND)) * (1.0D0/(1.0D0+(Kappa)**M)) &
                                    * ((-1)**lIndex) * ((Kappa)**(M-iSite+lIndex)) * PSI_A(jState,lIndex)

                PsiD(jState,iSite)= PsiD(jState,iSite) &
                                    + (1.0D0/(1.0D0+(Kappa)**M)) * ((-1)**lIndex) * ((Kappa)**(M-iSite+lIndex))  &
                                    * PSI_B(jState,lIndex)

             ENDDO !lIndex

             ELSE 
 
             DO kIndex=iSite+1,M
        
                PsiC(jState,iSite)= PsiC(jState,iSite) &
                                    + (CMPLX(OnsitePot,0.0D0,CKIND)) * (1.0D0/(1.0D0+(Kappa)**M)) &
                                    * ((-1)**kIndex) * ((Kappa)**(kIndex-iSite)) * PSI_A(jState,kIndex)

                PsiD(jState,iSite)= PsiD(jState,iSite) &
                                    + (1.0D0/(1.0D0+(Kappa)**M)) * ((-1)**kIndex) * ((Kappa)**(kIndex-iSite)) &
                                    * PSI_B(jState,kIndex)

             ENDDO !kIndex

             DO lIndex=1,iSite-1

                PsiC(jState,iSite)= PsiC(jState,iSite) &
                                    + (CMPLX(OnsitePot,0.0D0,CKIND)) * (1.0D0/(1.0D0+(Kappa)**M)) &
                                    * (-(-1)**lIndex) * ((Kappa)**(M-iSite+lIndex)) * PSI_A(jState,lIndex)

                PsiD(jState,iSite)= PsiD(jState,iSite) &
                                    + (1.0D0/(1.0D0+(Kappa)**M)) * (-(-1)**lIndex) * ((Kappa)**(M-iSite+lIndex)) &
                                    * PSI_B(jState,lIndex)

             ENDDO !lIndex
 
        ENDIF

        PsiC(jState,iSite)= PsiC(jState,iSite) &
                            + (1.0D0/(1.0D0+(Kappa)**M)) * (CMPLX(OnsitePot,0.0D0,CKIND)) * PSI_A(jState,iSite)
        PsiD(jState,iSite)= PsiD(jState,iSite) &
                            + (1.0D0/(1.0D0+(Kappa)**M)) * PSI_B(jState,iSite)

     ENDDO ! jState
  ENDDO ! iSite

  ! Copy the result in PSI_B(M,M)

  DO iSite=1,M
     DO jState=1,M
        PSI_B(jState,iSite)=PsiC(jState,iSite)-PsiD(jState,iSite)
     ENDDO
  ENDDO

  DEALLOCATE(PsiC)
  DEALLOCATE(PsiD)
  
  !PRINT*,"PSIA(1,1),(1,2),(1,3)",&
   !    PSI_A(1,1),PSI_A(1,2),PSI_A(1,3)
  !PRINT*,"PSIB(1,1),(1,2),(1,3)",&
   !    PSI_B(1,1),PSI_b(1,2),PSI_B(1,3)
  
  RETURN
END SUBROUTINE TMMultKACI

!--------------------------------------------------------------------
! TMMultKACII:
!
! Multiplication of the transfer matrix onto the vector (PSI_A,PSI_B), 
! giving (PSI_A,PSI_B) so that the structure of the transfer matrix 
! can be exploited
!-----------------------------------------------------------------------

SUBROUTINE TMMultKACII(PSI_A,PSI_B, Ilayer, En, DiagDis, Kappa, M )

  USE MyNumbers
  USE IPara
  USE RNG
  USE DPara
  
  ! wave functions:
  !       
  ! (PSI_A, PSI_B) on input, (PSI_B,PSI_A) on output
  
  IMPLICIT NONE
  
  INTEGER Ilayer,           &! current # TM multiplications
       M                     ! strip width
  
  REAL(KIND=RKIND)  DiagDis, &! diagonal disorder
       En,                  &! energy
       Kappa                 ! diagonal hopping
  
  COMPLEX(KIND=CKIND) PSI_A(M,M), PSI_B(M,M)

  INTEGER iSite, jState, IErr
  REAL(KIND=RKIND) OnsitePot
  COMPLEX(KIND=CKIND) new, PsiLeft

  ! define a new matrix where store data

  COMPLEX(KIND=CKIND), DIMENSION(:,:), ALLOCATABLE :: PsiC

  ALLOCATE(PsiC(M,M), STAT = IErr)

  ! reset to zero
  PsiC= ZERO
  
  DO iSite=1,M
     
     ! create the new onsite potential
     SELECT CASE(IRNGFlag)
     CASE(0)
        OnsitePot= En - DiagDis*(DRANDOM()-0.5D0)
     CASE(1)
        OnsitePot= En - DiagDis*(DRANDOM()-0.5D0)*SQRT(12.0D0)
     !CASE(2)
        !OnsitePot= En - GRANDOM(ISeedDummy,0.0D0,DiagDis)
     END SELECT
     

     DO jState=1,M

        ! Evaluation of PsiLeft
        
        IF (iSite.EQ.1) THEN
           PsiLeft=PSI_A(jState,M)
        ELSE
           PsiLeft=PSI_A(jState,iSite-1)
        ENDIF
        
        new=( CMPLX(OnsitePot,0.0D0,CKIND) * PSI_B(jState,iSite) &
             - Kappa*PsiLeft - PSI_A(jState,iSite) )

        PsiC(jState,iSite)=new

     ENDDO ! jState
  ENDDO ! iSite

  ! Copy the result in PSI_A(M,M)

  DO iSite=1,M
     DO jState=1,M
        PSI_A(jState,iSite)=PsiC(jState,iSite)
     ENDDO
  ENDDO

  DEALLOCATE(PsiC)
  
  !PRINT*,"PSIA(1,1),(1,2),(1,3)",&
   !    PSI_A(1,1),PSI_A(1,2),PSI_A(1,3)
  !PRINT*,"PSIB(1,1),(1,2),(1,3)",&
   !    PSI_B(1,1),PSI_b(1,2),PSI_B(1,3)
  
  RETURN
END SUBROUTINE TMMultKACII

!--------------------------------------------------------------------
! TMMultKACIII:
!
! Multiplication of the transfer matrix onto the vector (PSI_A,PSI_B), 
! giving (PSI_A,PSI_B) so that the structure of the transfer matrix 
! can be exploited
!-----------------------------------------------------------------------

SUBROUTINE TMMultKACIII(PSI_A,PSI_B, Ilayer, En, DiagDis, Kappa, M )

  USE MyNumbers
  USE IPara
  USE RNG
  USE DPara
  
  ! wave functions:
  !       
  ! (PSI_A, PSI_B) on input, (PSI_B,PSI_A) on output
  
  IMPLICIT NONE
  
  INTEGER Ilayer,           &! current # TM multiplications
       M                     ! strip width
  
  REAL(KIND=RKIND)  DiagDis, &! diagonal disorder
       En,                & ! energy
       Kappa               !diagonal hopping

  
  COMPLEX(KIND=CKIND) PSI_A(M,M), PSI_B(M,M)
  
  INTEGER iSite, jState, kIndex, lIndex, IErr
  REAL(KIND=RKIND) OnsitePot

  ! define a new matrix where store data

  COMPLEX(KIND=CKIND), DIMENSION(:,:), ALLOCATABLE :: PsiC
  COMPLEX(KIND=CKIND), DIMENSION(:,:), ALLOCATABLE :: PsiD

  ALLOCATE(PsiC(M,M), STAT = IErr)
  ALLOCATE(PsiD(M,M), STAT = IErr)

  ! reset to zero
  PsiC= ZERO
  PsiD= ZERO
  
  DO iSite=1,M
     
     ! create the new onsite potential
     SELECT CASE(IRNGFlag)
     CASE(0)
        OnsitePot= En - DiagDis*(DRANDOM()-0.5D0)
     CASE(1)
        OnsitePot= En - DiagDis*(DRANDOM()-0.5D0)*SQRT(12.0D0)
     !CASE(2)
        !OnsitePot= En - GRANDOM(ISeedDummy,0.0D0,DiagDis)
     END SELECT
     
     DO jState=1,M

        ! First term of the equation is stored in PsiC and second one in PsiD
        
        IF (iSite.EQ.M) THEN
           
           DO kIndex=1,M-1

              PsiC(jState,iSite)= PsiC(jState,iSite) &
                                  + (CMPLX(OnsitePot,0.0D0,CKIND)) * (1.0D0/(1.0D0+(Kappa)**M)) &
                                  * (-(-1)**kIndex) * ((Kappa)**(M-kIndex)) * PSI_A(jState,kIndex)

              PsiD(jState,iSite)= PsiD(jState,iSite) &
                                  + (1.0D0/(1.0D0+(Kappa)**M)) * (-(-1)**kIndex) * ((Kappa)**(M-kIndex)) &
                                  * PSI_B(jState,kIndex)

           ENDDO ! kIndex

        ELSE IF (MOD(iSite,2).EQ.1) THEN !iSite is odd

             DO kIndex=1,iSite-1
        
                PsiC(jState,iSite)= PsiC(jState,iSite) &
                                    + (CMPLX(OnsitePot,0.0D0,CKIND)) * (1.0D0/(1.0D0+(Kappa)**M)) &
                                    * (-(-1)**kIndex) * ((Kappa)**(iSite-kIndex)) * PSI_A(jState,kIndex)

                PsiD(jState,iSite)= PsiD(jState,iSite) &
                                    + (1.0D0/(1.0D0+(Kappa)**M)) * (-(-1)**kIndex) * ((Kappa)**(iSite-kIndex)) &
                                    * PSI_B(jState,kIndex)

             ENDDO !kIndex

             DO lIndex=iSite+1,M

                PsiC(jState,iSite)= PsiC(jState,iSite) &
                                    + (CMPLX(OnsitePot,0.0D0,CKIND)) * (1.0D0/(1.0D0+(Kappa)**M)) &
                                    * ((-1)**lIndex) * ((Kappa)**(M+iSite-lIndex)) * PSI_A(jState,lIndex)

                PsiD(jState,iSite)= PsiD(jState,iSite) &
                                    + (1.0D0/(1.0D0+(Kappa)**M)) * ((-1)**lIndex) * ((Kappa)**(M+iSite-lIndex)) & 
                                    * PSI_B(jState,lIndex)

             ENDDO !lIndex

             ELSE !iSite is even
 
             DO kIndex=1,iSite-1
        
                PsiC(jState,iSite)= PsiC(jState,iSite) &
                                   + (CMPLX(OnsitePot,0.0D0,CKIND)) * (1.0D0/(1.0D0+(Kappa)**M)) &
                                   * ((-1)**kIndex) * ((Kappa)**(iSite-kIndex)) * PSI_A(jState,kIndex)

                PsiD(jState,iSite)= PsiD(jState,iSite) &
                                    + (1.0D0/(1.0D0+(Kappa)**M)) * ((-1)**kIndex) * ((Kappa)**(iSite-kIndex)) &
                                    * PSI_B(jState,kIndex)

             ENDDO !kIndex

             DO lIndex=1+iSite,M

                PsiC(jState,iSite)= PsiC(jState,iSite) &
                                    + (CMPLX(OnsitePot,0.0D0,CKIND)) * (1.0D0/(1.0D0+(Kappa)**M)) &
                                    * (-(-1)**lIndex) * ((Kappa)**(M+iSite-lIndex)) * PSI_A(jState,lIndex)

                PsiD(jState,iSite)= PsiD(jState,iSite) &
                                    + (1.0D0/(1.0D0+(Kappa)**M)) * (-(-1)**lIndex) * ((Kappa)**(M+iSite-lIndex)) &
                                    * PSI_B(jState,lIndex)

             ENDDO !lIndex
 
        ENDIF

        PsiC(jState,iSite)= PsiC(jState,iSite) & 
                            + (1.0D0/(1.0D0+(Kappa)**M)) * (CMPLX(OnsitePot,0.0D0,CKIND))*PSI_A(jState,iSite)
        PsiD(jState,iSite)= PsiD(jState,iSite) &
                            + (1.0D0/(1.0D0+(Kappa)**M)) * PSI_B(jState,iSite)

     ENDDO ! jState
  ENDDO ! iSite

 ! Copy the result in PSI_B(M,M)

  DO iSite=1,M
     DO jState=1,M
        PSI_B(jState,iSite)=PsiC(jState,iSite)-PsiD(jState,iSite)
     ENDDO
  ENDDO

  DEALLOCATE(PsiC)
  DEALLOCATE(PsiD)
  
  !PRINT*,"PSIA(1,1),(1,2),(1,3)",&
   !    PSI_A(1,1),PSI_A(1,2),PSI_A(1,3)
  !PRINT*,"PSIB(1,1),(1,2),(1,3)",&
   !    PSI_B(1,1),PSI_b(1,2),PSI_B(1,3)
  
  RETURN
END SUBROUTINE TMMultKACIII

!--------------------------------------------------------------------
! TMMultKACIV:
!
! Multiplication of the transfer matrix onto the vector (PSI_A,PSI_B), 
! giving (PSI_A,PSI_B) so that the structure of the transfer matrix 
! can be exploited
!-----------------------------------------------------------------------

SUBROUTINE TMMultKACIV(PSI_A,PSI_B, Ilayer, En, DiagDis, Kappa, M )

  USE MyNumbers
  USE IPara
  USE RNG
  USE DPara
  
  ! wave functions:
  !       
  ! (PSI_A, PSI_B) on input, (PSI_B,PSI_A) on output
  
  IMPLICIT NONE
  
  INTEGER Ilayer,           &! current # TM multiplications
       M                     ! strip width
  
  REAL(KIND=RKIND)  DiagDis, &! diagonal disorder
       En,                 & ! energy
       Kappa                ! diagonal hopping

  
  COMPLEX(KIND=CKIND) PSI_A(M,M), PSI_B(M,M)

  INTEGER iSite, jState, IErr
  REAL(KIND=RKIND) OnsitePot
  COMPLEX(KIND=CKIND) new, PsiRight

  ! define a new matrix where store data

  COMPLEX(KIND=CKIND), DIMENSION(:,:), ALLOCATABLE :: PsiC

  ALLOCATE(PsiC(M,M), STAT = IErr)

  ! reset to zero
  PsiC= ZERO
  
  DO iSite=1,M
     
     ! create the new onsite potential
     SELECT CASE(IRNGFlag)
     CASE(0)
        OnsitePot= En - DiagDis*(DRANDOM()-0.5D0)
     CASE(1)
        OnsitePot= En - DiagDis*(DRANDOM()-0.5D0)*SQRT(12.0D0)
     !CASE(2)
        !OnsitePot= En - GRANDOM(ISeedDummy,0.0D0,DiagDis)
     END SELECT

     DO jState=1,M

        ! Evaluation of PsiLeft
        
        IF (iSite.EQ.M) THEN
           PsiRight=PSI_A(jState,1)
        ELSE
           PsiRight=PSI_A(jState,iSite+1)
        ENDIF
        
        new=( CMPLX(OnsitePot,0.0D0,CKIND) * PSI_B(jState,iSite) &
             - Kappa*PsiRight - PSI_A(jState,iSite) )

        PsiC(jState,iSite)=new

     ENDDO ! jState
  ENDDO ! iSite

  ! Copy the result in PSI_A(M,M)

  DO iSite=1,M
     DO jState=1,M
        PSI_A(jState,iSite)=PsiC(jState,iSite)
     ENDDO
  ENDDO

  DEALLOCATE(PsiC)
  
  !PRINT*,"PSIA(1,1),(1,2),(1,3)",&
   !    PSI_A(1,1),PSI_A(1,2),PSI_A(1,3)
  !PRINT*,"PSIB(1,1),(1,2),(1,3)",&
   !    PSI_B(1,1),PSI_b(1,2),PSI_B(1,3)
  
  RETURN
END SUBROUTINE TMMultKACIV


!--------------------------------------------------------------------
! TMMultACHWI:
!
! Multiplication of the transfer matrix onto the vector (PSI_A,PSI_B), 
! giving (PSI_A,PSI_B) so that the structure of the transfer matrix 
! can be exploited
!-----------------------------------------------------------------------

SUBROUTINE TMMultACHWI(PSI_A,PSI_B, Ilayer, En, DiagDis, M )

  USE MyNumbers
  USE IPara
  USE RNG
  USE DPara
  
  ! wave functions:
  !       
  ! (PSI_A, PSI_B) on input, (PSI_B,PSI_A) on output
  
  IMPLICIT NONE
  
  INTEGER Ilayer,           &! current # TM multiplications
       M                     ! strip width
  
  REAL(KIND=RKIND)  DiagDis, &! diagonal disorder
       En                  ! energy

  
  COMPLEX(KIND=CKIND) PSI_A(M,M), PSI_B(M,M)
  
  INTEGER iSite, jState, kIndex, lIndex, IErr, index
  REAL(KIND=RKIND) ax

  ! define a new matrix where store data

  COMPLEX(KIND=CKIND), DIMENSION(:,:), ALLOCATABLE :: PsiC
  COMPLEX(KIND=CKIND), DIMENSION(:,:), ALLOCATABLE :: PsiD
  REAL(KIND=RKIND), DIMENSION(:), ALLOCATABLE :: RndVec
  REAL(KIND=RKIND), DIMENSION(:), ALLOCATABLE :: OnsitePot

  ALLOCATE(PsiC(M,M), STAT = IErr)
  ALLOCATE(PsiD(M,M), STAT = IErr)
  ALLOCATE(RndVec(M), STAT = IErr)
  ALLOCATE(OnsitePot(M), STAT = IErr)

  ! reset to zero
  PsiC= ZERO
  PsiD= ZERO

  ! create the random numbers

  DO index=1,M
     ax=DRANDOM()
     RndVec(index)=ax
  ENDDO

  ! create the new onsite potential
 
  DO index=1,M
     SELECT CASE(IRNGFlag)
     CASE(0)
        OnsitePot(index)= En - DiagDis*(RndVec(index)-0.5D0)
     CASE(1)
        OnsitePot(index)= En - DiagDis*(RndVec(index)-0.5D0)*SQRT(12.0D0)
     !CASE(2)
        !OnsitePot(index)= En - GRANDOM(ISeedDummy,0.0D0,DiagDis)
     END SELECT
  ENDDO

  DO iSite=1,M
     
     DO jState=1,M

        ! First term of the equation is stored in PsiC and second one in PsiB
        
        IF (iSite.EQ.1) THEN
           
           DO kIndex=1,M

              PsiC(jState,iSite)= PsiC(jState,iSite) &
                                  + (OnsitePot(kIndex)) * (-(-1)**kIndex) &
                                  * PSI_A(jState,kIndex)

              PsiD(jState,iSite)= PsiD(jState,iSite) &
                                  + (-(-1)**kIndex) * PSI_B(jState,kIndex)

           ENDDO ! kIndex

        ELSE IF (MOD(iSite,2).EQ.1) THEN !iSite Odd

             DO kIndex=iSite,M
        
                PsiC(jState,iSite)= PsiC(jState,iSite) &
                                   + (OnsitePot(kIndex)) * (-(-1)**kIndex) &
                                   * PSI_A(jState,kIndex)

                PsiD(jState,iSite)= PsiD(jState,iSite) &
                                    + (-(-1)**kIndex) * PSI_B(jState,kIndex)

             ENDDO !kIndex

             ELSE !iSite Even
 
             DO kIndex=iSite,M 
        
                PsiC(jState,iSite)= PsiC(jState,iSite) &
                                    + (OnsitePot(kIndex)) * ((-1)**kIndex) &
                                    * PSI_A(jState,kIndex)

                PsiD(jState,iSite)= PsiD(jState,iSite) &
                                    + ((-1)**kIndex) * PSI_B(jState,kIndex)

             ENDDO !kIndex
 
        ENDIF

     ENDDO ! jState
  ENDDO ! iSite

  ! Copy the result in PSI_B(M,M)

  DO iSite=1,M
     DO jState=1,M
        PSI_B(jState,iSite)=PsiC(jState,iSite)-PsiD(jState,iSite)
     ENDDO
  ENDDO

  DEALLOCATE(PsiC)
  DEALLOCATE(PsiD)
  DEALLOCATE(RndVec)
  DEALLOCATE(OnsitePot)
  
  !PRINT*,"First step"
  !PRINT*,"PSIA(1,1),(2,1),(3,1),(4,1)",&
   !    PSI_A(1,1),PSI_A(2,1),PSI_A(3,1),PSI_A(4,1)
  !PRINT*,"PSIB(1,1),(2,1),(3,1),(4,1)",&
   !    PSI_B(1,1),PSI_B(2,1),PSI_B(3,1),PSI_B(4,1)
  
  RETURN
END SUBROUTINE TMMultACHWI

!--------------------------------------------------------------------
! TMMultACHWII:
!
! Multiplication of the transfer matrix onto the vector (PSI_A,PSI_B), 
! giving (PSI_A,PSI_B) so that the structure of the transfer matrix 
! can be exploited
!-----------------------------------------------------------------------

SUBROUTINE TMMultACHWII(PSI_A,PSI_B, Ilayer, En, DiagDis, M )

  USE MyNumbers
  USE IPara
  USE RNG
  USE DPara
  
  ! wave functions:
  !       
  ! (PSI_A, PSI_B) on input, (PSI_B,PSI_A) on output
  
  IMPLICIT NONE
  
  INTEGER Ilayer,           &! current # TM multiplications
       M                     ! strip width
  
  REAL(KIND=RKIND)  DiagDis, &! diagonal disorder
       En                  ! energy

  
  COMPLEX(KIND=CKIND) PSI_A(M,M), PSI_B(M,M)

  INTEGER iSite, jState, IErr
  REAL(KIND=RKIND) OnsitePot
  COMPLEX(KIND=CKIND) new, PsiLeft

  ! define a new matrix where store data

  COMPLEX(KIND=CKIND), DIMENSION(:,:), ALLOCATABLE :: PsiC

  ALLOCATE(PsiC(M,M), STAT = IErr)

  ! reset to zero
  PsiC= ZERO
  
  DO iSite=1,M
     
     ! create the new onsite potential
     SELECT CASE(IRNGFlag)
     CASE(0)
        OnsitePot= En - DiagDis*(DRANDOM()-0.5D0)
     CASE(1)
        OnsitePot= En - DiagDis*(DRANDOM()-0.5D0)*SQRT(12.0D0)
     !CASE(2)
        !OnsitePot= En - GRANDOM(ISeedDummy,0.0D0,DiagDis)
     END SELECT
     
     DO jState=1,M

        ! Evaluation of PsiLeft
        
        IF (iSite.EQ.1) THEN
           PsiLeft=ZERO
        ELSE
           PsiLeft=PSI_A(jState,iSite-1)
        ENDIF
        
        new=( CMPLX(OnsitePot,0.0D0,CKIND) * PSI_B(jState,iSite) &
             - PsiLeft - PSI_A(jState,iSite) )

        PsiC(jState,iSite)=new

     ENDDO ! jState
  ENDDO ! iSite

  ! Copy the result in PSI_A(M,M)

  DO iSite=1,M
     DO jState=1,M
        PSI_A(jState,iSite)=PsiC(jState,iSite)
     ENDDO
  ENDDO

  DEALLOCATE(PsiC)
  
 !PRINT*,"Second step"
 !PRINT*,"PSIA(1,1),(2,1),(3,1),(4,1)",&
  !     PSI_A(1,1),PSI_A(2,1),PSI_A(3,1),PSI_A(4,1)
  !PRINT*,"PSIB(1,1),(2,1),(3,1),(4,1)",&
   !    PSI_B(1,1),PSI_B(2,1),PSI_B(3,1),PSI_B(4,1)
  
  RETURN
END SUBROUTINE TMMultACHWII

!--------------------------------------------------------------------
! TMMultACHWIII:
!
! Multiplication of the transfer matrix onto the vector (PSI_A,PSI_B), 
! giving (PSI_A,PSI_B) so that the structure of the transfer matrix 
! can be exploited
!-----------------------------------------------------------------------

SUBROUTINE TMMultACHWIII(PSI_A,PSI_B, Ilayer, En, DiagDis, M )

  USE MyNumbers
  USE IPara
  USE RNG
  USE DPara
  
  ! wave functions:
  !       
  ! (PSI_A, PSI_B) on input, (PSI_B,PSI_A) on output
  
  IMPLICIT NONE
  
  INTEGER Ilayer,           &! current # TM multiplications
       M                     ! strip width
  
  REAL(KIND=RKIND)  DiagDis, &! diagonal disorder
       En                 ! energy

  
  COMPLEX(KIND=CKIND) PSI_A(M,M), PSI_B(M,M)
  
  INTEGER iSite, jState, kIndex, lIndex, IErr, index
  REAL(KIND=RKIND) ax

  ! define a new matrix where store data

  COMPLEX(KIND=CKIND), DIMENSION(:,:), ALLOCATABLE :: PsiC
  COMPLEX(KIND=CKIND), DIMENSION(:,:), ALLOCATABLE :: PsiD
  REAL(KIND=RKIND), DIMENSION(:), ALLOCATABLE :: RndVec
  REAL(KIND=RKIND), DIMENSION(:), ALLOCATABLE :: OnsitePot

  ALLOCATE(PsiC(M,M), STAT = IErr)
  ALLOCATE(PsiD(M,M), STAT = IErr)
  ALLOCATE(RndVec(M), STAT = IErr)
  ALLOCATE(OnsitePot(M), STAT = IErr)

  ! reset to zero
  PsiC= ZERO
  PsiD= ZERO

  ! create the random numbers

  DO index=1,M
     ax=DRANDOM()
     RndVec(index)=ax
  ENDDO

  ! create the new onsite potential
 
  DO index=1,M
     SELECT CASE(IRNGFlag)
     CASE(0)
        OnsitePot(index)= En - DiagDis*(RndVec(index)-0.5D0)
     CASE(1)
        OnsitePot(index)= En - DiagDis*(RndVec(index)-0.5D0)*SQRT(12.0D0)
     !CASE(2)
        !OnsitePot(index)= En - GRANDOM(ISeedDummy,0.0D0,DiagDis)
     END SELECT
  ENDDO
  
  DO iSite=1,M

     DO jState=1,M

        ! First term of the equation is stored in PsiC and second one in PsiD
        
        IF (iSite.EQ.M) THEN 
           
           DO kIndex=1,M

              PsiC(jState,iSite)= PsiC(jState,iSite) &
                                  + (OnsitePot(kIndex)) * ((-1)**(M+kIndex)) &
                                  * PSI_A(jState,kIndex)

              PsiD(jState,iSite)= PsiD(jState,iSite) &
                                  + ((-1)**(kIndex+M)) * PSI_B(jState,kIndex)

           ENDDO ! kIndex

        ELSE IF (MOD(iSite,2).EQ.1) THEN !iSite Odd

             DO kIndex=1,iSite
        
                PsiC(jState,iSite)= PsiC(jState,iSite) &
                                    + (OnsitePot(kIndex)) * (-(-1)**kIndex) &
                                    * PSI_A(jState,kIndex)

                PsiD(jState,iSite)= PsiD(jState,iSite) &
                                    + (-(-1)**kIndex) * PSI_B(jState,kIndex)

             ENDDO !kIndex

             ELSE !iSite Even
 
             DO kIndex=1,iSite
        
                PsiC(jState,iSite)= PsiC(jState,iSite) &
                                   + (OnsitePot(kIndex)) * ((-1)**kIndex) &
                                   * PSI_A(jState,kIndex)

                PsiD(jState,iSite)= PsiD(jState,iSite) &
                                    + ((-1)**kIndex) * PSI_B(jState,kIndex)

             ENDDO !kIndex
 
        ENDIF

     ENDDO ! jState
  ENDDO ! iSite

 ! Copy the result in PSI_B(M,M)

  DO iSite=1,M
     DO jState=1,M
        PSI_B(jState,iSite)=PsiC(jState,iSite)-PsiD(jState,iSite)
     ENDDO
  ENDDO

  DEALLOCATE(PsiC)
  DEALLOCATE(PsiD)
  DEALLOCATE(RndVec)
  DEALLOCATE(OnsitePot)

  !PRINT*, "Third step"
  !PRINT*,"PSIA(1,1),(2,1),(3,1),(4,1)",&
   !    PSI_A(1,1),PSI_A(2,1),PSI_A(3,1),PSI_A(4,1)
  !PRINT*,"PSIB(1,1),(2,1),(3,1),(4,1)",&
   !    PSI_B(1,1),PSI_B(2,1),PSI_B(3,1),PSI_B(4,1)
  
  RETURN
END SUBROUTINE TMMultACHWIII

!--------------------------------------------------------------------
! TMMultACHWIV:
!
! Multiplication of the transfer matrix onto the vector (PSI_A,PSI_B), 
! giving (PSI_A,PSI_B) so that the structure of the transfer matrix 
! can be exploited
!-----------------------------------------------------------------------

SUBROUTINE TMMultACHWIV(PSI_A,PSI_B, Ilayer, En, DiagDis, M )

  USE MyNumbers
  USE IPara
  USE RNG
  USE DPara
  
  ! wave functions:
  !       
  ! (PSI_A, PSI_B) on input, (PSI_B,PSI_A) on output
  
  IMPLICIT NONE
  
  INTEGER Ilayer,           &! current # TM multiplications
       M                     ! strip width
  
  REAL(KIND=RKIND)  DiagDis, &! diagonal disorder
       En                  ! energy

  
  COMPLEX(KIND=CKIND) PSI_A(M,M), PSI_B(M,M)

  INTEGER iSite, jState, IErr
  REAL(KIND=RKIND) OnsitePot
  COMPLEX(KIND=CKIND) new, PsiRight

  ! define a new matrix where store data

  COMPLEX(KIND=CKIND), DIMENSION(:,:), ALLOCATABLE :: PsiC

  ALLOCATE(PsiC(M,M), STAT = IErr)

  ! reset to zero
  PsiC= ZERO
  
  DO iSite=1,M
     
     ! create the new onsite potential
     SELECT CASE(IRNGFlag)
     CASE(0)
        OnsitePot= En - DiagDis*(DRANDOM()-0.5D0)
     CASE(1)
        OnsitePot= En - DiagDis*(DRANDOM()-0.5D0)*SQRT(12.0D0)
     !CASE(2)
        !OnsitePot= En - GRANDOM(ISeedDummy,0.0D0,DiagDis)
     END SELECT
     
     DO jState=1,M

        ! Evaluation of PsiLeft
        
        IF (iSite.EQ.M) THEN
           PsiRight=ZERO
        ELSE
           PsiRight=PSI_A(jState,iSite+1)
        ENDIF
        
        new=( CMPLX(OnsitePot,0.0D0,CKIND) * PSI_B(jState,iSite) &
             - PsiRight - PSI_A(jState,iSite) )

        PsiC(jState,iSite)=new

     ENDDO ! jState
  ENDDO ! iSite

  ! Copy the result in PSI_A(M,M)

  DO iSite=1,M
     DO jState=1,M
        PSI_A(jState,iSite)=PsiC(jState,iSite)
     ENDDO
  ENDDO

  DEALLOCATE(PsiC)
  
  !PRINT*, "Fourth step"
  !PRINT*,"PSIA(1,1),(2,1),(3,1),(4,1)",&
   !    PSI_A(1,1),PSI_A(2,1),PSI_A(3,1),PSI_A(4,1)
  !PRINT*,"PSIB(1,1),(2,1),(3,1),(4,1)",&
   !    PSI_B(1,1),PSI_B(2,1),PSI_B(3,1),PSI_B(4,1)
  
  RETURN
END SUBROUTINE TMMultACHWIV

!------------------------------------------------------------------------
!	ReNorm:
!
!	Gram-Schmidt orthonormalization of the wave vectors (PSI_A,PSI_B),
!	see e.g. Horn/Johnson, "Matrix Analysis", pp. 15. 
!
!	PLUS reading off the Lyapunov exponents
!
!	(PSI_A,PSI_B) is the incoming vector of eigenvectors,
!	GAMMA, GAMMA2 are Lyapunov exponent and its square,
!	M is the width of the strip,
!
!	IVec/JVec label different vectors, 
!	KIndex is the component index of a SINGLE vector,
!
!	dummy, sum and norm are just different names for the same
!	double precision workspace
!-------------------------------------------------------------------------

SUBROUTINE ReNorm(PSI_A,PSI_B,GAMMA,GAMMA2,M)
  
  USE MyNumbers
  USE IConstants
  USE IPara
  
  INTEGER M
  
  COMPLEX(KIND=CKIND) PSI_A(M,M), PSI_B(M,M)
  REAL(KIND=RKIND) GAMMA(M), GAMMA2(M)
  
  INTEGER IVec,JVec,KIndex
  
  COMPLEX(KIND=CKIND) sum
  REAL(KIND=RKIND) dummy,norm
  EQUIVALENCE (dummy,norm)

  !make the local variables static
  !SAVE
  
  !PRINT*,"DBG: ReNorm()"
  
  DO 100 IVec=1,M
     
     DO 200 JVec=1,IVec-1
        
        sum= CZERO
        
        DO 300 KIndex=1,M
           
           sum= sum + CONJG(PSI_A(JVec,KIndex))*PSI_A(IVec,KIndex) &
                + CONJG(PSI_B(JVec,KIndex))*PSI_B(IVec,KIndex)
300     ENDDO
        
        DO 400 KIndex=1,M
           
           PSI_A(IVec,KIndex)= PSI_A(IVec,KIndex) - &
                sum * PSI_A(JVec,KIndex)
           PSI_B(IVec,KIndex)= PSI_B(IVec,KIndex) - &
                sum * PSI_B(JVec,KIndex)
           
400     ENDDO
        
200  ENDDO
     
     ! calculation of norm
     
     norm= REAL(0.D0,RKIND)
     DO 500 KIndex=1,M                      
        norm= norm + CONJG(PSI_A(IVec,KIndex)) * PSI_A(IVec,KIndex) &
             + CONJG(PSI_B(IVec,KIndex)) * PSI_B(IVec,KIndex)
500  ENDDO
     dummy= 1.D0/SQRT(norm)
     DO 600 KIndex=1,M
        PSI_A(IVec,KIndex)= CMPLX(dummy,0.0D0,CKIND) * PSI_A(IVec,KIndex)
        PSI_B(IVec,KIndex)= CMPLX(dummy,0.0D0,CKIND) * PSI_B(IVec,KIndex)
600  ENDDO
     
     !	----------------------------------------------------------------
     !	gammadummy is ordered s.t. the vector with the smallest overall 
     !	norm is last and the vector with the largest overall norm first. 
     !	We sum this so that the smallest value is used for GAMMA(M)
     !	and the largest for GAMMA(1). Same for GAMMA2(). Therefore, the
     !	inverse of the smallest Lyapunov exponent is obtained by taking
     !	LAMBDA(1) = N/GAMMA(M)
     dummy       = LOG(dummy)
     GAMMA(IVec) = GAMMA(IVec) - dummy
     GAMMA2(IVec)= GAMMA2(IVec) + dummy*dummy
          
     ! ----------------------------------------------------------------
     ! check orthogonality if desired
     GOTO 100
     IF(IWriteFlag.GE.MAXWriteFlag) THEN

        DO JVec=1,IVec-1
           sum= REAL(0.D0,RKIND)
           DO KIndex=1,M
              sum= sum + CONJG(PSI_A(JVec,KIndex))*PSI_A(IVec,KIndex) &
                   + CONJG(PSI_B(JVec,KIndex))*PSI_B(IVec,KIndex)
           ENDDO
           PRINT*,"Renorm: <",JVec,"|",IVec,">=",sum
        ENDDO
        
     ENDIF
     
100 ENDDO
  
  RETURN
  
END SUBROUTINE ReNorm

!---------------------------------------------------------------------
!	Swap:
!
!	(PSI_A,PSI_B)= (old,new) is the incoming vector, this is swapped
!	into (PSI_A,PSI_B)= (new,old)
!-------------------------------------------------------------------

SUBROUTINE Swap( PSI_A, PSI_B, M)

  USE MyNumbers
  
  INTEGER M
  REAL(KIND=RKIND) PSI_A(M,M), PSI_B(M,M)
  
  INTEGER jState, index
  REAL(KIND=RKIND) dummy
  
  !	PRINT*,"DBG: Swap()"
  
  DO jState=1,M
     DO index=1,M
        
        dummy              = PSI_B(index,jState)
        PSI_B(index,jState)= PSI_A(index,jState)
        PSI_A(index,jState)= dummy
        
     ENDDO
  ENDDO
  
  RETURN
  
END SUBROUTINE Swap


!--------------------------------------------------------------------
!	ReSort:
!
!	sort the Lyapunov eigenvalues s.t. the largest comes first. RESORT()
!	is ShellSort taken from NumRec, SHELL().
!---------------------------------------------------------------------

SUBROUTINE ReSort( PSI_A, PSI_B, array0, array1, N )

  USE MyNumbers
  
  INTEGER N
  COMPLEX(KIND=CKIND) PSI_A(N,N),PSI_B(N,N)
  REAL(KIND=RKIND) array0(N), array1(N)
  
  REAL(KIND=RKIND) ALN2I, LocalTINY
  PARAMETER (ALN2I=1.4426950D0, LocalTINY=1.D-5)
  
  INTEGER NN,M,L,K,J,I,LOGNB2, index
  REAL(KIND=RKIND) dummyA, dummyB
  
  !	PRINT*,"DBG: ReSort()"
  
  !	PRINT*,"array0(1),array0(N)",array0(1),array0(N)
  
  LOGNB2=INT(LOG(REAL(N))*ALN2I+LocalTINY)
  M=N
  DO 12 NN=1,LOGNB2
     M=M/2
     K=N-M
     DO 11 J=1,K
        I=J
3       CONTINUE
        L=I+M
        IF(array0(L).GT.array0(I)) THEN
           
           dummyA   = array0(I)
           array0(I)= array0(L)
           array0(L)= dummyA
           
           dummyB   = array1(I)
           array1(I)= array1(L)
           array1(L)= dummyB
           
           DO 100 index=1,N
              dummyA        = PSI_A(index,I)
              dummyB        = PSI_B(index,I)
              
              PSI_A(index,I)= PSI_A(index,L)
              PSI_B(index,I)= PSI_B(index,L)
              
              PSI_A(index,L)= dummyA
              PSI_B(index,L)= dummyB
100        ENDDO
        
           I=I-M
           IF(I.GE.1) GOTO 3
        ENDIF
11   ENDDO
12 ENDDO
  
  !	PRINT*,"array0(1),array0(N)",array0(1),array0(N)
  RETURN

END SUBROUTINE ReSort

! --------------------------------------------------------------------
SUBROUTINE CheckUnitarity(Umat,Wmat,M,IErr)

  USE MyNumbers
  USE IPara
  USE IConstants
  USE DPara

  IMPLICIT NONE

  INTEGER(KIND=IKIND) M, IErr, index,jndex
  COMPLEX(KIND=CKIND) Umat(M,M), Wmat(M,M)
  REAL(KIND=RKIND) sum

  IErr=0

  !PRINT*,"Umat=", Umat

  Wmat=CONJG(TRANSPOSE(Umat))
  !Wmat=CONJG(Wmat)
  Wmat=MATMUL(Umat,Wmat)

  sum= 0.0D0
  DO index=1,M
     DO jndex=1,M
        sum= sum+ABS(Wmat(index,jndex))
     ENDDO
  ENDDO
  sum= (sum-REAL(M,RKIND))/REAL(M*M,RKIND)

  !PRINT*,"Wmat=", Wmat
  IF(IWriteFlag.GE.MAXWriteFlag) THEN
     PRINT*,"CheckUnit: sum=",sum!; PAUSE
  ENDIF

  !IF(ABS(sum).GE.EPSILON(0.0D0)) IErr=1
  IF(ABS(sum).GE.TINY) THEN
     PRINT*,"CheckUnit: sum=",sum," > TINY=",TINY," !"
     IErr=1
  ENDIF

  ! check symmetrization if no magnetic field
  IF(MagFlux.LT.TINY) THEN
     Wmat=TRANSPOSE(Umat)
     
     sum= 0.0D0
     DO index=1,M
        DO jndex=1,M
           sum= sum+ABS(Umat(index,jndex)-Wmat(index,jndex))
        ENDDO
     ENDDO
     sum= (sum)/REAL(M*M,RKIND)
     
     !PRINT*,"Wmat=", Wmat
     IF(IWriteFlag.GE.MAXWriteFlag) THEN
        PRINT*,"CheckSym: sum=",sum!; PAUSE
     ENDIF
  ENDIF

  RETURN

END SUBROUTINE CheckUnitarity
