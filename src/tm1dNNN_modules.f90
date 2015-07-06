!********************************************************************
!
! TM1DNNN - Transfer matrix method for the Anderson model
! in 1D with finite-range hopping
!
!********************************************************************
       
!********************************************************************
!
! $Header: /home/cvs/phsht/GrapheneTMM/src/Restart/ALL/tmseGR_modules.f90,v 1.1 2011/07/22 17:49:19 ccspam Exp $
!
!********************************************************************

!--------------------------------------------------------------------
MODULE MyNumbers     
  IMPLICIT NONE

  INTEGER, PARAMETER :: IKIND = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: RKIND = SELECTED_REAL_KIND(15,307)
  INTEGER, PARAMETER :: CKIND = SELECTED_REAL_KIND(15,307)

  REAL(KIND=RKIND) :: PI, TWOPI, ONEPLS, ONEMNS

  REAL(KIND=RKIND), PARAMETER :: ZERO = 0.0, ONE = 1.0 ,TWO = 2.0, THREE = 3.0, FOUR = 4.0
  COMPLEX(KIND=RKIND), PARAMETER :: CZERO = (0.0d0,0.0d0), CONE = (1.0d0,0.0d0), &
       CIMAGONE= (0.0d0,1.0d0)            

  REAL (KIND=RKIND), PARAMETER :: HALF = 0.5D0, QUARTER = 0.25D0, EIGHTH = 0.125D0

  REAL(KIND=RKIND) :: TINY= 1.0D-9
  
CONTAINS
  SUBROUTINE INIT_NUMBERS
    PI     = 4.0D0* ATAN(1.0D0)
    TWOPI  = 8.0D0* ATAN(1.0D0)
    ONEMNS = SQRT(EPSILON(ONEMNS))
    ONEPLS = ONE + ONEMNS
    ONEMNS = ONE - ONEMNS
  END SUBROUTINE INIT_NUMBERS

  FUNCTION ARG(X,Y)
    
    REAL(KIND=RKIND) ARG, X, Y
    
    IF( X > 0. ) THEN 
       ARG= ATAN(Y/X)
    ELSE IF ( (X == 0.) .and. (Y > 0. )) THEN 
       ARG = PI/2.0D0
    ELSE IF ( (X == 0.) .and. (Y < 0. )) THEN 
       ARG = -PI/2.0D0
    ELSE IF ( (X < 0. ) .and. (Y >= 0.)) THEN 
       ARG = PI + ATAN(Y/X)
    ELSE IF ( (X < 0. ) .and. (Y < 0. )) THEN 
       ARG = -PI + ATAN(Y/X)
    ENDIF
    
    RETURN
  END FUNCTION ARG

END MODULE MyNumbers

!--------------------------------------------------------------------------------------------
MODULE CConstants
  CHARACTER*18, PARAMETER :: RStr= "$Revision: 0.1 $ "
  CHARACTER*30, PARAMETER :: DStr= "$Date: 2015/06/01 17:49:19 $ "
  CHARACTER*16, PARAMETER :: AStr= "$Author: phsht $ "
END MODULE CConstants

!------------------------------------------------------------------------------------------
! MAXGamma needs to be equal to MAXRange, as we need to find ALL
! Lyapunov exponents, so do not change!
MODULE IConstants
  INTEGER, PARAMETER :: MAXRange= 1000, MAXGamma= MAXRange, MAXIter=2147483646
  INTEGER, PARAMETER :: MAXKeepFlag= 3, MAXWriteFlag= 5, MAXFluxFlag= 3, MAXRNGFlag=2
  INTEGER, PARAMETER :: MAXSortFlag=1, MAXModelFlag=2, MAXBCFlag=2
  INTEGER, PARAMETER :: MINDimenFlag=2, MAXDimenFlag=3
  INTEGER, PARAMETER :: MAXFiles= 5, MINIter=3
END MODULE IConstants

!--------------------------------------------------------------------
MODULE IPara
  INTEGER :: ISeed, NOfIter, NOfOrtho, NOfPrint, NOfGamma
  INTEGER :: IRange0, IRange1, dIRange, IKeepFlag, IWriteFlag, ISortFlag
  INTEGER :: IFluxFlag, IBCFlag, IModelFlag, IRNGFlag, IKappaFlag
  INTEGER :: IRestart, IMidDump, IWalltime
  INTEGER :: index1, Iter1, IRange
END MODULE IPara

!--------------------------------------------------------------------
MODULE DPara
  USE MyNumbers
  REAL(KIND=RKIND) :: DiagDis0,DiagDis1,dDiagDis,DiagDis
  REAL(KIND=RKIND) :: Energy0,Energy1,dEnergy,Energy
  REAL(KIND=RKIND) :: Kappa,MagFlux, Epsilon
  
  !REAL(KIND=RKIND), DIMENSION(:,:), ALLOCATABLE :: PsiA, PsiB
END MODULE DPara

!--------------------------------------------------------------------
!      Input- and Outputchannels
MODULE IChannels
  INTEGER, PARAMETER :: IChInp= 40, IChOut= 41, IChOutGam= 42, &
       IChOutPsi= 43, IChOutRHO=44,IChOutWrite= 39, IChOutAvgRHO= 45, &
       IChLoadLP= 46, IChSaveLP= 47, IChRes= 48, &
       IChOutAvgRHO1= 50, IChOutAvgRHOL= 51, &
       ICHtmp= 52, IChOutHAV= 53
END MODULE IChannels

!------------------------------------------------------------------
MODULE Randoms
  USE MyNumbers
  INTEGER(KIND=IKIND) :: z1, z2, z3, z4
END MODULE Randoms

!------------------------------------------------------------------
MODULE Extensions
  USE MyNumbers
  REAL(KIND=RKIND), ALLOCATABLE, DIMENSION(:,:) :: PsiA, PsiB
END MODULE Extensions

!-------------------------------------------------------------------------
MODULE Gammas
  USE MyNumbers
  REAL(KIND=RKIND), ALLOCATABLE, DIMENSION(:) :: nGamma, gamma, gamma2, acc_variance
END MODULE Gammas





