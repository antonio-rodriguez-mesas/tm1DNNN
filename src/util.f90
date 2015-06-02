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

SUBROUTINE TMMultNNN(PSI_A,PSI_B, Ilayer, &
     HopMatiLR, HopMatiL, EMat, dummyMat, OnsitePotVec, &
     En, Dis, M )

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
  
  REAL(KIND=RKIND)  Dis,    &! diagonal disorder
       En                    ! energy
  
  COMPLEX(KIND=CKIND) PSI_A(M,M), PSI_B(M,M)
  REAL(KIND=RKIND) HopMatiLR(M,M),HopMatiL(M,M), EMat(M,M),dummyMat(M,M)
  REAL(KIND=RKIND) OnsitePotVec(M)
  
  INTEGER iSite, jState
  REAL(KIND=RKIND) OnsitePot
  COMPLEX(KIND=CKIND) new, PsiLeft, PsiRight

  PRINT*,"TMMultNNN()"

  ! create the new onsite matrix 
  DO iSite=1,M
     SELECT CASE(IRNGFlag)
     CASE(0)
        EMat(iSite,iSite)= En - Dis*(DRANDOM()-0.5D0)
     CASE(1)
        EMat(iSite,iSite)= En - Dis*(DRANDOM()-0.5D0)*SQRT(12.0D0)
     !CASE(2)
        !EMat(iSite,iSite)= En - GRANDOM(ISeedDummy,0.0D0,Dis)
     END SELECT
  ENDDO ! iSite

  dummyMat= MATMUL(HopMatiL,EMat)
  dummyMat= MATMUL(dummyMat,PSI_A)

  PSI_A= dummyMat + MATMUL(HopMatiLR,PSI_B)

  !PRINT*,"PSIA(1,1),(2,1),(3,1),(4,1)",&
   !    PSI_A(1,1),PSI_A(2,1),PSI_A(3,1),PSI_A(4,1)
  !PRINT*,"PSIB(1,1),(2,1),(3,1),(4,1)",&
   !    PSI_B(1,1),PSI_B(2,1),PSI_B(3,1),PSI_B(4,1)
  
  RETURN
END SUBROUTINE TMMultNNN

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
