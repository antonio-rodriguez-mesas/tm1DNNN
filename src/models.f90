! ********************************************************************
!       
! TM1DNNN - Transfer matrix method for the Anderson model
! in 1D with finite-range hopping
!
! ********************************************************************
       
MODULE Models

USE MyNumbers
USE DPara

CONTAINS
  
  REAL (RKIND) FUNCTION HoppingModel(&
       IFlag, IDistance )
    
    INTEGER(IKIND) IFlag, IDistance

    SELECT CASE(IFlag)
    CASE(-1)
       HoppingModel=REAL(IDistance,RKIND)
    CASE(0)
       HoppingModel=1.0_RKIND
    CASE(1)
       HoppingModel=1.0_RKIND/REAL(IDistance,RKIND)
    CASE(2)
       HoppingModel=EXP(-REAL(IDistance,RKIND)/Kappa)
    CASE(3)
       HoppingModel=1.0_RKIND/(REAL(IDistance,RKIND)**Kappa)
    CASE DEFAULT
       PRINT*,"HoppingModel(): IFlag=", IFlag, " not implemented --- aborting!"
       STOP
    END SELECT
    
    RETURN
  END FUNCTION HoppingModel

END MODULE Models
