!********************************************************************
!
! TM1DNNN - Transfer matrix method for the Anderson model
! in 1D with finite-range hopping
!
!********************************************************************

SUBROUTINE INVERT(MatrixSize,Matrix,InvertedMatrix,IErr)  
!
!   Invert an M*M Complex Matrix
!   Matrix: the Matrix (Destroyed)
!   InvertedMatrix: the Inverse
!  USE WriteToScreen
  USE MyNumbers

  USE IConstants
  USE IPara

!  USE MyMPI
  
  IMPLICIT NONE
  
  INTEGER :: MatrixSize, LWORK, INFO, I, IErr
  REAL(KIND=RKIND), DIMENSION(1:MatrixSize,1:MatrixSize) :: Matrix
  REAL(KIND=RKIND), DIMENSION(1:MatrixSize,1:MatrixSize) :: InvertedMatrix
  
  INTEGER, DIMENSION(:), ALLOCATABLE :: IPIV
  REAL(KIND=RKIND), DIMENSION(:), ALLOCATABLE :: WORK
  
  !PRINT*,"DBG: Invert()"
!!$
!!$  B = CZERO 
!!$  DO I=1,M
!!$     B(I,I) = CONE
!!$  END DO
  !INFO=0

  ALLOCATE(IPIV(MatrixSize),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Invert(): ERR in ALLOCATE(IPIV(MatrixSize)) statement, MatrixSize=", &
          MatrixSize
     RETURN
  ENDIF
  
  CALL DGETRF(MatrixSize,MatrixSize,Matrix,MatrixSize,IPIV,IErr)
  LWORK = MatrixSize*MatrixSize
  !LWORK = 0
  IF ( IErr.NE.0 ) THEN
     PRINT *,'Invert() : Datatype Error: IFAIL=',INFO
     RETURN
  END IF
  ALLOCATE(WORK(LWORK),STAT=IErr)   
  IF( IErr.NE.0 ) THEN
     PRINT*,"Invert(): ERR in ALLOCATE(WORK(LWORK)) statement, LWORK=", LWORK
     RETURN
  ENDIF
  
  CALL DGETRI(MatrixSize,Matrix,MatrixSize,IPIV,WORK,LWORK,IErr)
  IF ( IErr.NE.0 ) THEN
     PRINT *,'Inversion Error: IFAIL=',INFO
     RETURN
  END IF
  DEALLOCATE(IPIV,WORK,STAT=IErr)
  IF ( IErr.NE.0 ) THEN
     PRINT *,'Invert : Deallocation Error',INFO
     RETURN
  END IF
  !DEALLOCATE(IPIV,WORK)
  InvertedMatrix = Matrix  
  RETURN
END SUBROUTINE INVERT


