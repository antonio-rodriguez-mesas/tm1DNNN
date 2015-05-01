C       ********************************************************************
C       
C       TMSE3DAO - Transfer matrix method for the 3D Anderson model with 
C                  automatic reorthogonalization
C
C       ********************************************************************
       
C       ********************************************************************
C       
C       $Header: /home/kerekou1/rar/f77/tmse3dAO/RCS/common.h,v 2.2 1998/06/04 07:06:12 rar Exp rar $
C
C       ********************************************************************

C **************************************************************************
C
C $Log: common.h,v $
C * Revision 2.2  1998/06/04  07:06:12  rar
C * included "MINNOrtho= 5" to provide the minimal starting value for IOrtho,
C * the input parameter NOfOrtho is now the maximum allowed value of IOrtho
C *
C * Revision 2.1  1998/04/29  13:50:48  rar
C * BLAS version
C *
C * Revision 1.3  1998/04/02  08:40:42  rar
C * removed tmOD remnants
C *
C * Revision 1.2  1998/03/30  14:19:26  rar
C * new comments
C *
C * Revision 1.1  1998/03/30  14:15:27  rar
C * Initial revision
C *
C * Revision 1.4  96/12/10  12:39:18  12:39:18  rar (Rudolf Roemer)
C * DBLE->REAL*8
C * 
C * Revision 1.3  1996/12/04  11:19:47  rar
C * included a diagonal disorder potential (DiagDis0,DiagDis1,dDiagDis)
C *
C * Revision 1.2  96/11/04  11:33:47  11:33:47  rar (Rudolf Roemer)
C * changed the header statement
C * 
C * Revision 1.1  96/11/01  16:39:56  16:39:56  rar (Rudolf Roemer)
C * Initial revision
C * 
C * Revision 3.2  96/04/22  14:28:59  14:28:59  rar (Rudolf Roemer)
C * made TINY global, introduced R/D/AStr to make logging in output files
C * with version number for later comparison, changed MAXInitFlag to 2
C * 
c Revision 3.1  96/04/16  11:36:06  11:36:06  rar (Rudolf Roemer)
c initial RCS version of V2.5
c 
C
C **************************************************************************


C       ********************************************************************
C
C       History prior to version 3.0:
C
C       15/04/96 RAR: adjustments for V 2.5
C       16/02/96 RAR: NOfConfig0/NOfConfig1/NOfIter
C       06/02/96 RAR: taken from TMSE2D
C       24/01/96 RAR: perfectly new installation
C
C       ********************************************************************

C	--------------------------------------------------------------------
        CHARACTER*17 RStr
        CHARACTER*30 DStr
        CHARACTER*15 AStr

        PARAMETER(
     +     RStr= "$Revision: 2.2 $ ",
     +     DStr= "$Date: 1998/06/04 07:06:12 $ ",
     +     AStr= "$Author: rar $ "
     +  )

        INTEGER MAXWidth, MAXGamma, MINNOrtho

C       MAXGamma needs to be equal to MAXWidth, as we need to find ALL
C       Lyapunov exponents, so do not change!

        PARAMETER (MAXWidth= 50, MAXGamma= MAXWidth*MAXWidth,
     +             MINNOrtho= 5)

C	--------------------------------------------------------------------
        INTEGER MAXKeepFlag, MAXWriteFlag, MAXSortFlag,
     +     MAXBCXFlag, MAXBCYFlag, MAXFiles, MINIter, MAXIMWrite
        PARAMETER (MAXKeepFlag=1,MAXWriteFlag= 2,
     +     MAXSortFlag=1, MAXBCXFlag=1, MAXBCYFlag=1, MAXFiles= 2, MINIter=3,
     +     MAXIMWrite=2)

C	--------------------------------------------------------------------
        INTEGER ISeed, NOfIter, NOfOrtho, NOfGamma,
     +     WidthX0, WidthX1, dWidthX, WidthY, IKeepFlag,
     +     IWriteFlag, ISortFlag, IBCXFlag, IBCYFlag, IMatWrite

        COMMON/IPara/ISeed, NOfIter, NOfOrtho, NOfGamma,
     +     WidthX0, WidthX1, dWidthX, WidthY, IKeepFlag,
     +     IWriteFlag, ISortFlag,IBCXFlag, IBCYFlag, IMatWrite

C	--------------------------------------------------------------------
        REAL*8
     +     DiagDis0,DiagDis1,dDiagDis,
     +     Energy, Epsilon

        REAL*8 PsiA(MAXWidth*MAXWidth*MAXWidth*MAXWidth),
     +         PsiB(MAXWidth*MAXWidth*MAXWidth*MAXWidth)

        COMMON/DPara/
     +     DiagDis0,DiagDis1,dDiagDis,
     +     Energy,  Epsilon,
     +     PsiA, PsiB

C	--------------------------------------------------------------------
C       Input- and Outputchannels

        INTEGER IChInp, IChOut, IChOutGam, IChOutPsi
        PARAMETER (IChInp= 40, IChOut= 41, IChOutGam= 42, IChOutPsi= 43)

C	--------------------------------------------------------------------
C       constants
C	--------------------------------------------------------------------

        REAL*8 TINY
        PARAMETER( TINY= 1.0D-20 )

















