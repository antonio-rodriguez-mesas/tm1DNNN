!********************************************************************
!
! TMSEXD - Transfer matrix method for the Anderson
! model with diagonal disorder in X dimensions
!
!********************************************************************
       
!********************************************************************
!
! $Header: /home/cvs/phsht/GrapheneTMM/src/Restart/ALL/tmseGR_modules.f90,v 1.1 2011/07/22 17:49:19 ccspam Exp $
!
!********************************************************************

!**************************************************************************
!$Log: tmseGR_modules.f90,v $
!Revision 1.1  2011/07/22 17:49:19  ccspam
!Programs for ZZ and AC with the restart routine for francesca
!
!Revision 1.2  2011/05/31 13:53:43  ccspam
!*** empty log message ***
!
!Revision 1.1  2011/05/06 08:13:09  phsht
!1st installement
!
!Revision 1.1  2010/11/11 11:16:25  phrkaj
!Renamed files to get rid of any mention of SB.
!
!
!**************************************************************************

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
  CHARACTER*18, PARAMETER :: RStr= "$Revision: 1.1 $ "
  CHARACTER*30, PARAMETER :: DStr= "$Date: 2011/07/22 17:49:19 $ "
  CHARACTER*16, PARAMETER :: AStr= "$Author: ccspam $ "
END MODULE CConstants

!------------------------------------------------------------------------------------------
! MAXGamma needs to be equal to MAXWidth, as we need to find ALL
! Lyapunov exponents, so do not change!
MODULE IConstants
  INTEGER, PARAMETER :: MAXWidth= 1000, MAXGamma= MAXWidth, MAXIter=2147483646
  INTEGER, PARAMETER :: MAXKeepFlag= 3, MAXWriteFlag= 4, MAXFluxFlag= 3, MAXRNGFlag=2
  INTEGER, PARAMETER :: MAXSortFlag=1, MAXBCFlag=2
  INTEGER, PARAMETER :: MINDimenFlag=2, MAXDimenFlag=3
  INTEGER, PARAMETER :: MAXFiles= 5, MINIter=3
END MODULE IConstants

!--------------------------------------------------------------------
MODULE IPara
  INTEGER :: ISeed, NOfIter, NOfOrtho, NOfPrint, NOfGamma
  INTEGER :: Width0, Width1, dWidth, IKeepFlag, IWriteFlag, ISortFlag
  INTEGER :: IFluxFlag, IBCFlag, IRNGFlag, IKappaFlag, IEdgeFlag
  INTEGER :: IRestart, IMidDump, IWalltime
  INTEGER :: index1, Iter1, IWidth
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
  COMPLEX(KIND=CKIND), ALLOCATABLE, DIMENSION(:,:) :: PsiA, PsiB
END MODULE Extensions

!-------------------------------------------------------------------------
MODULE Gammas
  USE MyNumbers
  REAL(KIND=RKIND), ALLOCATABLE, DIMENSION(:) :: nGamma, gamma, gamma2, acc_variance
END MODULE Gammas





