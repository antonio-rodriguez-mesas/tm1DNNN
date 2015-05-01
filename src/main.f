C       ********************************************************************
C       
C       dfcc3d - Transfer matrix method for the 3D Anderson model with
C                  automatic reorthogonalization
c                  FCC lattice, TMM along the main diagonal                 
C
C       ********************************************************************
       
C       ********************************************************************
C
C       $Header: /home/kerekou1/rar/f77/tmse3dAO/RCS/main.f,v 2.2 1998/06/04 07:04:46 rar Exp $
C
C       ********************************************************************

C **************************************************************************
C
C $Log: main.f,v $
c Revision 2.2  1998/06/04  07:04:46  rar
c made NOfOrtho the maximum number of IOrtho allowed
c
c Revision 2.1  1998/04/29  13:50:58  rar
c BLAS version
c
c Revision 1.8  1998/04/02  08:43:02  rar
c added a "D"
c
c Revision 1.7  1998/04/02  08:39:32  rar
c slight changes, removed old variables
c
c Revision 1.6  1998/04/01  16:00:47  crv
c removed DisCenter, DisWidth
c
c Revision 1.5  1998/04/01  09:43:44  crv
c remove PresentVec, PastVec and HopFuture
c
c Revision 1.2  1998/03/30  14:19:35  rar
c new comments
c
c Revision 1.1  1998/03/30  14:16:03  rar
c Initial revision
c
c Revision 1.6  1996/12/10  12:39:50  rar
c DBLE->REAL*8
c
c Revision 1.5  1996/12/04  11:20:29  rar
c included a diagonal disorder potential (DiagDis0,DiagDis1,dDiagDis)
c
c Revision 1.4  1996/11/14  11:13:19  rar
c TINY is now passed on to TMMult() to check for division by zero;
c
c Revision 1.3  96/11/04  11:42:27  11:42:27  rar (Rudolf Roemer)
c removed pLevel
c 
c Revision 1.2  96/11/04  11:28:31  11:28:31  rar (Rudolf Roemer)
c changed some comments in the header
c 
c Revision 1.1  96/11/01  16:40:20  16:40:20  rar (Rudolf Roemer)
c Initial revision
c 
C
C **************************************************************************

C       ********************************************************************
C
C       Comments:
C
C       IKeepFlag       0        overwrite old .raw files
C                       1        keep old .raw files and skip to next
C                                configuration
C
C       IWriteFlag      0        no protocolling
C                       1        write out final wave function for each
C                                strip width and flux (disorder/U)
C                       2        additionally, write out each nGamma
C                                value after reorthonormalisation
C
C       IMatWrite       0        no protocolling
C                       1        write the connectivity matrix and its inverse
C                                to '*.mat'
C                       2        as 2, stop after writing
C
C
C       ISortFlag       0        NO sorting of eigenvalues/vectors after
C                                each reorthonormalization
C                       1        YES, sorting is done.
C
C       boundary conditions along the X direction:
C
C       IBCXFlag        0        hard wall boundary conditions 
C                       1        helical boundary conditions
C
C
C       boundary conditions along the Y direction
C
C       IBCYFlag        0        hard wall boundary conditions
C                       1        periodic boundary conditions
C
C
C       Notation:
C
C       see B. Kramer and M. Schreiber, "Transfer-Matrix Methods and Finite-
C	Size Scaling for Disordered Systems" and K. Frahm et al., EPL 31,
C       169 (1995).
C
C       Also, Martin Hennecke, "Anderson Uebergang im Magnetfeld", PTB
C       Bericht, PTB-PG-6, Braunschweig, Dec 94.
C
C       ********************************************************************

        PROGRAM dfcc3d

        IMPLICIT NONE

C	--------------------------------------------------------------------
C       parameter and global variable definitions
C	--------------------------------------------------------------------

	INCLUDE "common.h"

C	--------------------------------------------------------------------
C       local variable definitions
C	--------------------------------------------------------------------

        INTEGER IWidthX,IWidthY,MXY, index, Iter1,Iter2, NOfG,iG,
     +         IterNum, IOrtho

        REAL*8 flux,flux0,flux1,dflux, DiagDis, DIter1

        REAL*8 nGamma(MAXGamma), gamma(MAXGamma),
     +         gamma2(MAXGamma), acc_variance(MAXGamma)

        INTEGER IErr
        integer MAX4,MAX2

c       variables needed to invert the connectivity matrix

        parameter (MAX4=MAXWidth*MAXWidth*MAXWidth*MAXWidth)
	parameter (MAX2=MAXWidth*MAXWidth)
	real*8 con(MAX4),coninv(MAX4),righth(MAX4),af(MAX4)
	real*8 berr(MAX2),ferr(MAX2),c(MAX2),r(MAX2)
	real*8 work(3*MAX2)
	integer  ipiv(MAX2),iwork(MAX2)




C	--------------------------------------------------------------------
C       constants
C	--------------------------------------------------------------------


C	--------------------------------------------------------------------
C	protocal feature startup
C	--------------------------------------------------------------------

	PRINT*,"TMSE3DAO ASYMP (double precision)"
C       PRINT*,"$Revision: 2.2 $ $Date: 1998/06/04 07:04:46 $ $Author: rar $"
        PRINT*,RStr,DStr,AStr
C       PRINT*,"$Revision: 2.2 $ $Date: 1998/06/04 07:04:46 $ $Author: rar $"
        PRINT*,"----------------------------------------------------------"

C	--------------------------------------------------------------------
C	input handling
C	--------------------------------------------------------------------

        CALL Input( IErr )
D	PRINT*, "DBG: IErr=", IErr
	IF( IErr.NE.0 ) THEN
	  PRINT*,"main: error in Input()"
	  STOP
	ENDIF

        IF (ISortFlag.EQ.1) THEN
           PRINT*,"ReSort() is done after each renormalization !"
           PRINT*,"-------------------------------------------------------"
        ENDIF

C	--------------------------------------------------------------------
C       select the quantity to be varied
C	--------------------------------------------------------------------

        flux0     = DiagDis0
        flux1     = DiagDis1
        dflux     = dDiagDis

C	--------------------------------------------------------------------
C       main parameter sweep
C	--------------------------------------------------------------------

        DO 1000 IWidthX= WidthX0,WidthX1,dWidthX

C          --------------------------------------------------------------
C          decide on MxM bars or MXxMY bars
C          --------------------------------------------------------------
           IF (WidthY.NE.WidthX0) THEN
              IWidthY= WidthY
           ELSE
              IWidthY= IWidthX
           ENDIF

c          construct and invert the connectivity matrix


           CALL dcopy(MAX4, 0.0D0, 0, con, 1)
	   CALL dcopy(MAX4, 0.0D0, 0, righth, 1)
	   call cmatinver(IWidthX,IWidthY,con,coninv,righth,af,
     $                    berr,ferr,c,r,work,ipiv,iwork,IBCXFlag,
     $                    IBCYFlag,IMatWrite)

     	   if (IMatWrite.eq.2) goto 1000


C	   --------------------------------------------------------------
C          the # Lyapunov exponents is maximally .EQ. to IWidthX*IWidthY
C	   --------------------------------------------------------------

           NOfG= MIN( NOfGamma, IWidthX*IWidthY )

C	   --------------------------------------------------------------
C	   open the AVG file
C	   --------------------------------------------------------------

           IF (IKeepFlag.EQ.1 ) THEN
              CALL CheckOutputAvg( IWidthX,IWidthY, IErr )
D             PRINT*, "DBG: IErr=", IErr
              IF( IErr.EQ.1 ) THEN
                 PRINT*,"main: error in CheckOutputAvg()"
                 STOP
              ELSE IF (IErr.EQ.2 ) THEN
                 GOTO 1000
              ENDIF
           ENDIF

           CALL OpenOutputAvg( IWidthX,IWidthY, IErr )
D          PRINT*, "DBG: IErr=", IErr
           IF( IErr.EQ.1 ) THEN
              PRINT*,"main: error in OpenOutputAvg()"
              STOP
           ENDIF

C	   --------------------------------------------------------------
C	   flux loop
C	   --------------------------------------------------------------

           DO 2000 flux= flux0,flux1,dflux

C	    --------------------------------------------------------------
C	    set values for the physical quantities
C	    --------------------------------------------------------------

              DiagDis   = flux

C	      -----------------------------------------------------------
C	      open the nGamma file
C	      -----------------------------------------------------------

              IF(IWriteFlag.GE.1) THEN
                 CALL OpenOutputGamma( IWidthX,IWidthY,
     +                DiagDis, IErr )
D                PRINT*, "DBG: IErr=", IErr
                 IF( IErr.NE.0 ) THEN
                    PRINT*,"main: error in OpenOutputGamma()"
                    STOP
                 ENDIF
              ENDIF

 2500         WRITE(*,2510) IWidthX,IWidthY, DiagDis
 2510         FORMAT("START @ IWidthX= ",I4.1,
     +             ", IWidthY= ",I4.1,
     +             ", DiagDis= ", G10.3)

C	      -----------------------------------------------------------
C             initialize the wave vectors and the gamme sums
C	      -----------------------------------------------------------

              MXY = IWidthX*IWidthY

C             reset the gamma sum approximants for a SINGLE
C             configuration

              CALL dcopy(MXY, 0.0D0, 0, gamma , 1)
              CALL dcopy(MXY, 0.0D0, 0, gamma2, 1)

C             reset the wave vectors

              CALL dcopy(MAX4, 0.0D0, 0, PsiA, 1)
              CALL dcopy(MAX4, 0.0D0, 0, PsiB, 1)

              CALL dcopy(MXY, 1.0D0, 0, PsiA, MXY+1)

C             -----------------------------------------------------------
C             initialize the random number generator. NOTE that this is
C             done with the same ISeed for each value of IWidth and
C             Disorder in order to make these runs INDIVIDUALLY reprodu-
C             cible
C	      -----------------------------------------------------------

              CALL SRANDOM(ISeed)

C             -----------------------------------------------------------
C             initialize IOrtho
C             -----------------------------------------------------------

              IOrtho= MIN(NOfOrtho,MINNOrtho)

C             -----------------------------------------------------------
C	      iteration loop
C	      -----------------------------------------------------------

C                 DO 4000 Iter1= 1,NOfIter, 1

                 Iter1   = 0
                 IterNum = 0

 4000            Iter1   = Iter1 +1

                 DO 4100 Iter2= 1, IOrtho-1, 2

C                   do the 2 TM multiplications

                    CALL TMMult( PsiA,PsiB,Righth,
     +                   Energy, DiagDis, IWidthX,IWidthY,
     +                   IBCXFlag, IBCYFlag, coninv)

                    CALL TMMult( PsiB,PsiA,Righth,
     +                   Energy, DiagDis, IWidthX,IWidthY,
     +                   IBCXFlag, IBCYFlag, coninv)

                    IterNum = IterNum + 2

 4100            CONTINUE

C                -------------------------------------------------------
C                renormalize via Gram-Schmidt
C                -------------------------------------------------------

D                PRINT *,"Iter: ",Iter1," IOrtho: ", IOrtho

                 CALL ReOrtho(PsiA,PsiB,gamma,gamma2,
     +                        IWidthX,IWidthY,IOrtho)

C                -------------------------------------------------------
C                check if IOrtho is too large, then reset
C                -------------------------------------------------------

                 IF(IOrtho.GT.NOfOrtho) THEN
                    IOrtho= NOfOrtho
                 ENDIF

C                -------------------------------------------------------
C                sort the eigenvalues by LARGEST first AND also
C                resort the eigenvectors accordingly
C                -------------------------------------------------------

                 IF (ISortFlag.EQ.1) THEN
                    CALL ReSort(PsiA,PsiB,gamma,gamma2,IWidthX,IWidthY)
                 ENDIF

C	         --------------------------------------------------------
C                "Iter1" counts the number of actual renormalizaions
C                of the transfer matrix.
C	         --------------------------------------------------------

                 DIter1= DBLE(Iter1)

                 DO 4200 iG=1, MXY

                    nGamma(MXY+1-iG)= gamma(iG)/REAL(IterNum)

                    acc_variance(MXY+1-iG)=
     +                   SQRT( ABS(
     +                   (gamma2(iG)/DIter1 -
     +                   (gamma(iG)/DIter1)**2 )
     +                   / MAX(DIter1-1.D0,1.D0)
     +                   )) / ABS( gamma(iG)/DIter1 )

 4200            CONTINUE

C	         -----------------------------------------------------------
C	         write the nGamma data to file
C	         -----------------------------------------------------------

                 IF(IWriteFlag.GE.1) THEN
                    CALL WriteOutputGamma( Iter1, nGamma,
     +                   acc_variance,NOfG,
     +                   IErr )
D                   PRINT*, "DBG: IErr=", IErr
                    IF( IErr.NE.0 ) THEN
                       PRINT*,"main: error in WriteOutputGamma()"
                       STOP
                    ENDIF
                 ENDIF

C	         --------------------------------------------------------
C	         write the nGamma data to stdout
C	         --------------------------------------------------------

                 IF(IWriteFlag.GE.2) THEN
                    WRITE(*,4210) Iter1, IOrtho, nGamma(1), acc_variance(1)
 4210               FORMAT(I7.1, I7.1, G15.7, G15.7)
                 ENDIF

C                this system dependent subroutine flushes all open
C                buffers, note that small letters are needed for the
C                HP F77 compiler
                 CALL f77fflush()

C	         -----------------------------------------------------------
C                check accuracy and dump the result
C	         -----------------------------------------------------------

                 IF( acc_variance(1).LE.Epsilon  .AND. 
     +                acc_variance(1).GE.TINY    .AND.
     +                Iter1          .GE.MXY     .AND.
     +                Iter1          .GE.MINIter ) THEN
D                   PRINT*,"convergence in NOfIter-loop:"
                    GOTO 5000
                 ENDIF

C 4000            CONTINUE
                 IF(Iter1.LT.NOfIter) THEN
                    GOTO 4000
                 ENDIF

C             -------------------------------------------------------------
C             continue through here if convergence for a single configura-
C             tion is NOT achieved, reset Iter1
C             -------------------------------------------------------------

D             PRINT*,"NO convergence in NOfIter-loop:"
              Iter1= Iter1-1

C             -------------------------------------------------------------
C             jump to this label if convergence for a single configura-
C             tion is achieved. 
C             -------------------------------------------------------------

 5000         WRITE(*,5010) Iter1, IterNum, DiagDis,
     +             nGamma(1), acc_variance(1)
 5010         FORMAT("END @ ", I7.1, ",",I7.1,
     +             ",", G15.7, ",", G15.7, ",", G15.7, 
     +             ",", G15.7, ",", G15.7)

C	      --------------------------------------------------------------
C	      close the nGamma file
C	      --------------------------------------------------------------

              IF(IWriteFlag.GE.1) THEN
                 CALL CloseOutputGamma( IErr )
D                PRINT*, "DBG: IErr=", IErr
                 IF( IErr.NE.0 ) THEN
                    PRINT*,"main: error in CloseOutputGamma()"
                    STOP
                 ENDIF
              ENDIF

C	      --------------------------------------------------------------
C	      write the AVG data
C	      --------------------------------------------------------------

              CALL WriteOutputAvg( IWidthX,IWidthY, DiagDis, 
     +             nGamma, acc_variance, NOfG, PsiA, IErr )
D             PRINT*, "DBG: IErr=", IErr
              IF( IErr.NE.0 ) THEN
                 PRINT*,"main: error in WriteOutputAvg()"
                 STOP
              ENDIF

C             --------------------------------------------------------------
C             end of flux loop
C             --------------------------------------------------------------

 2000      CONTINUE

C	   -----------------------------------------------------------------
C	   close the AVG files
C	   -----------------------------------------------------------------

           CALL CloseOutputAvg( IErr )
D          PRINT*, "DBG: IErr=", IErr
           IF( IErr.NE.0 ) THEN
              PRINT*,"main: error in CloseOutputAvg()"
              STOP
           ENDIF

C          -----------------------------------------------------------------
C          end of width loop
C          -----------------------------------------------------------------

 1000   CONTINUE

        STOP "TMSE3DAO $Revision: 2.2 $"
        END
