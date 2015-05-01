C       ********************************************************************
C       
C       TM2E3DAO - Transfer matrix method for the 3D Anderson model with 
C                  automatic reorthogonalization
C
C       ********************************************************************
       
       
C       ********************************************************************
C       
C       $Header: /home/kerekou1/rar/f77/tmse3dAO/RCS/inout.f,v 2.1 1998/04/29 13:51:13 rar Exp $
C
C       ********************************************************************

C **************************************************************************
C
C $Log: inout.f,v $
c Revision 2.1  1998/04/29  13:51:13  rar
c BLAS version
c
c Revision 1.5  1998/04/03  13:49:55  rar
c swapped PSI_A/B(x,y) -> PSI_A/B(y,x)
c
c Revision 1.4  1998/04/02  09:10:42  rar
c made DEC compatible
c
c Revision 1.3  1998/04/02  08:39:57  rar
c removed old parameters, changed filenames
c
c Revision 1.2  1998/04/01  16:03:15  crv
c removed DisCenter & DisWidth, added IBCXFlag & IBCYFlag
c
c Revision 1.1  1998/03/30  14:17:35  rar
c Initial revision
c
c Revision 1.5  1996/12/10  12:39:56  rar
c DBLE->REAL*8
c
c Revision 1.4  1996/12/04  11:20:44  rar
c included a diagonal disorder potential (DiagDis0,DiagDis1,dDiagDis)
c
c Revision 1.3  1996/11/13  10:38:10  rar
c F90-linux compatible
c
c Revision 1.2  1996/11/04  11:33:32  rar
c changed the help output
c
c Revision 1.1  96/11/01  16:40:37  16:40:37  rar (Rudolf Roemer)
c Initial revision
c 
C
C **************************************************************************

C	--------------------------------------------------------------------
C	Input:
C
C	IErr	error code

	SUBROUTINE Input( IErr )

        INCLUDE "common.h"

        INTEGER IErr, ILine

D	PRINT*,"DBG: Input()"

	IErr = 0
        ILine= 0

        OPEN(UNIT= IChInp, ERR= 120, FILE= "dfcc3d.inp", STATUS= 'OLD')

        ILine= ILine+1
 	READ(IChInp,10,ERR=20,END=30) ISeed
D	PRINT*,"ISeed        = ",ISeed

        ILine= ILine+1
	READ(IChInp,10,ERR=20,END=30) NOfIter
D	PRINT*,"NOfIter      = ", NOfIter

        ILine= ILine+1
	READ(IChInp,10,ERR=20,END=30) NOfOrtho
D	PRINT*,"NOfOrtho     = ", NOfOrtho

        ILine= ILine+1
	READ(IChInp,10,ERR=20,END=30) NOfGamma
D	PRINT*,"NOfGamma     = ", NOfGamma

        ILine= ILine+1
	READ(IChInp,10,ERR=20,END=30) IKeepFlag
D	PRINT*,"IKeepFlag    = ", IKeepFlag

        ILine= ILine+1
	READ(IChInp,10,ERR=20,END=30) IWriteFlag
D	PRINT*,"IWriteFlag   = ", IWriteFlag

        ILine= ILine+1
	READ(IChInp,10,ERR=20,END=30) IMatWrite
D	PRINT*,"IMatWrite   = ", IMatWrite

        ILine= ILine+1
	READ(IChInp,10,ERR=20,END=30) ISortFlag
D	PRINT*,"ISortFlag    = ", ISortFlag

        ILine= ILine+1
	READ(IChInp,10,ERR=20,END=30) IBCXFlag
D	PRINT*,"IBCXFlag     = ", IBCXFlag

        ILine= ILine+1
	READ(IChInp,10,ERR=20,END=30) IBCYFlag
D	PRINT*,"IBCYFlag     =",  IBCYFlag

        ILine= ILine+1
	READ(IChInp,10,ERR=20,END=30) WidthX0
D	PRINT*,"WidthX0       = ",WidthX0

        ILine= ILine+1
	READ(IChInp,10,ERR=20,END=30) WidthX1
D	PRINT*,"WidthX1       = ", WidthX1

        ILine= ILine+1
	READ(IChInp,10,ERR=20,END=30) dWidthX
D	PRINT*,"dWidthX       = ", dWidthX

        ILine= ILine+1
	READ(IChInp,10,ERR=20,END=30) WidthY
D	PRINT*,"WidthY       = ",WidthY

        ILine= ILine+1
	READ(IChInp,15,ERR=20,END=30) DiagDis0
D	PRINT*,"DiagDis0     = ", DiagDis0

        ILine= ILine+1
	READ(IChInp,15,ERR=20,END=30) DiagDis1
D	PRINT*,"DiagDis1     = ", DiagDis1

        ILine= ILine+1
	READ(IChInp,15,ERR=20,END=30) dDiagDis
D	PRINT*,"dDiagDis     = ", dDiagDis

        ILine= ILine+1
	READ(IChInp,15,ERR=20,END=30) Energy
D	PRINT*,"Energy       = ", Energy

        ILine= ILine+1
	READ(IChInp,15,ERR=20,END=30) Epsilon
D	PRINT*,"Epsilon      = ", Epsilon

 10	FORMAT(16X,I10.1)
C 10	FORMAT("IMAXIteration= ",I10.1)
 15     FORMAT(16X,F18.9)
C 15     FORMAT("IMAXIteration= ",F18.9)

        CLOSE(IChInp, ERR=130)

C	check the parameters for validity

   	IF( NOfIter.LE.0 ) THEN
           PRINT*,"Input(): NOfIter <= 0"
           IErr= 1
	ENDIF

   	IF( NOfGamma.GT.MAXGamma ) THEN
           PRINT*,"Input(): NOfGamma > MAXGamma (=",MAXGamma,")"
           IErr= 1
	ENDIF

   	IF( WidthX0.LE.0 ) THEN
           PRINT*,"Input(): WidthX0 <= 0"
           IErr= 1
	ENDIF

   	IF( WidthX0.GT.MAXWidth ) THEN
           PRINT*,"Input(): WidthX0 > MAXWidth (=",MAXWidth,")"
           IErr= 1
	ENDIF

   	IF( WidthX1.GT.MAXWidth ) THEN
           PRINT*,"Input(): WidthX1 > MAXWidth (=",MAXWidth,")"
           IErr= 1
	ENDIF

   	IF( (WidthX0.GT.WidthX1) .AND. (dWidthX.GT.0) ) THEN
           PRINT*,"Input(): WidthX0 > WidthX1 and dWidthX>0"
           IErr= 1
	ENDIF
        IF( WidthY.LE.0 ) THEN
           PRINT*,"Input(): WidthY <= 0"
           IErr= 1
	ENDIF

   	IF( WidthY.GT.MAXWidth ) THEN
           PRINT*,"Input(): WidthY > MAXWidth (=",MAXWidth,")"
           IErr= 1
	ENDIF

   	IF( IKeepFlag.GT.MAXKeepFlag ) THEN
           PRINT*,"Input(): IKeepFlag > MAXKeepFlag (=",MAXKeepFlag,")"
           IErr= 1
	ENDIF

   	IF( IWriteFlag.GT.MAXWriteFlag ) THEN
           PRINT*,"Input(): IWriteFlag > MAXWriteFlag (=",MAXWriteFlag,")"
           IErr= 1
	ENDIF

	IF( IMatWrite.GT.MAXIMWrite ) THEN
           PRINT*,"Input(): IMatWrite > MAXIMWrite (=",MAXIMWrite,")"
           IErr= 1
	ENDIF

   	IF( ISortFlag.GT.MAXSortFlag ) THEN
           PRINT*,"Input(): ISortFlag > MAXSortFlag (=",MAXSortFlag,")"
           IErr= 1
	ENDIF

   	IF( IBCXFlag.GT.MAXBCXFlag ) THEN
           PRINT*,"Input(): IBCXFlag > MAXBCXFlag (=",MAXBCXFlag,")"
           IErr= 1
	ENDIF

   	IF( IBCYFlag.GT.MAXBCYFlag ) THEN
           PRINT*,"Input(): IBCYFlag > MAXBCYFlag (=",MAXBCYFlag,")"
           IErr= 1
	ENDIF

   	IF( IBCXFlag.EQ.IBCYFlag ) THEN
           PRINT*,"Input(): IBCXFlag = IBCYFlag =",IBCXFlag
           IErr= 0
	ENDIF

        IF( Epsilon.LE.0.0D0) THEN
           PRINT*,"Input(): nonpositive Epsilon"
           IErr= 1
        ENDIF

	RETURN

C	error in OPEN detected
 120	PRINT*,"Input(): ERR in OPEN"
        GOTO 1000

C	error in CLOSE detected
 130	PRINT*,"Input(): ERR in CLOSE"
        GOTO 1000

C	error in READ detected
 20	PRINT*,"Input(): ERR in READ at line", ILine
        GOTO 1000

C	EOF in READ occured prematurely
 30	PRINT*,"Input(): EOF in READ at line", ILine

C       dump the input help

 1000   PRINT*,"Input parameters:          ; explanation:"
        PRINT*,"--------------------------------------------------------------------"
        PRINT*,"ISeed         = 123456     ; seed for random number generator"

        PRINT*,"NOfIter       = 1          ; # steps for convergence of SINGLE config."
        PRINT*,"NOfOrtho      = 10         ; # steps until reorthonormalization"
        PRINT*,"NOfGamma      = 1          ; # smallest Lyapunov exponents"

        PRINT*,"IKeepFlag     = 0          ; 0/1 = yes/no data overwrite"
        PRINT*,"IWriteFlag    = 1          ; 0/2/1 = no/wave fcn/Gamma output"
	PRINT*,"IMatWrite     = 0          ; 0/1/2 no/yes write connectivity matrix/stop"
        PRINT*,"ISortFlag     = 0          ; 0/1 = no/yes ReSort()"
        PRINT*,"IBCXFlag      = 0          ; 0/1 = hard wall/helical b.c. along X-direction"
 	PRINT*,"IBCYFlag      = 0          ; 0/1 = hard wall/periodic b.c. along Y-direction"

        PRINT*,"WidthX0       = 0          ; minimal X width"
        PRINT*,"WidthX1       = 0          ; maximal X width"
        PRINT*,"dWidthX       = 0          ; X width increment"

        PRINT*,"WidthY        = 0          ; Y width"

        PRINT*,"DiagDis0      = 0.         ; minimal diagonal disorder"
        PRINT*,"DiagDis1      = 5.         ; maximal  diagonal disorder"
        PRINT*,"dDiagDis      = 0.5        ; increment of diagonal disorder"

        PRINT*,"Energy        = 0.0        ; energy goal"

        PRINT*,"Epsilon       = 5.0E-2     ; accuracy goal of iteration"

	IErr= 1
	RETURN

	END


C	--------------------------------------------------------------------
C	CheckOutputAvg:
C
C	IErr	error code

	SUBROUTINE CheckOutputAvg( IWidthX,IWidthY, IErr )

        INCLUDE "common.h"

        INTEGER IWidthX,IWidthY, IErr

        CHARACTER*12 CName

D	PRINT*,"DBG: CheckOutputAvg()"

	IErr= 0

C       WRITE out the input parameter
        
        WRITE(CName, '(I4.4,I4.4,A4)') IWidthX, IWidthY, '.raw'

        OPEN(UNIT= IChOut, ERR= 10, STATUS= 'NEW', FILE=CName)

        IErr= 0

 20     CLOSE(UNIT= IChOut, ERR= 100)

        RETURN

 10     WRITE(*,15) CName
 15     FORMAT(" OpenOutputAvg(): ", A12, 
     +       " exists -- skipped!")

        IErr= 2
        GOTO 20

C       ERR in CLOSE detected
 100	PRINT*,"CloseOutputAvg(): ERR in CLOSE()"
	IErr= 1
	RETURN
        
        END


C	--------------------------------------------------------------------
C	OpenOutputAvg:
C
C	IErr	error code

	SUBROUTINE OpenOutputAvg( IWidthX,IWidthY, IErr )

        INCLUDE "common.h"

        INTEGER IWidthX,IWidthY, IErr

        INTEGER ICh, IChList(MAXFiles)

        CHARACTER*12 CName

D	PRINT*,"DBG: OpenOutputAvg()"

	IErr= 0

C       WRITE out the input parameter
        
        WRITE(CName, '(I4.4,I4.4,A4)') IWidthX, IWidthY, '.raw'
        OPEN(UNIT= IChOut, ERR= 10, STATUS= 'UNKNOWN', FILE=CName)

        WRITE(CName, '(I4.4,I4.4,A4)') IWidthX, IWidthY, '.psi'
        OPEN(UNIT= IChOutPsi, ERR= 10, STATUS= 'UNKNOWN', FILE=CName)

        IChList(1)= IChOut
C        IChList(2)= IChOutPsi

        DO 1000 ICh= 1, 1

           WRITE(IChList(ICh),90,ERR=20) RStr,DStr,AStr
 90        FORMAT("(* ",3A," *)")

           WRITE(IChList(ICh),100,ERR=20) ISeed
 100       FORMAT("ISeed        = ", I10.1, "; (* + Config *)")

           WRITE(IChList(ICh),120,ERR=20) NOfIter
 120       FORMAT("NOfIter      = ", I10.1, ";")

           WRITE(IChList(ICh),140,ERR=20) NOfOrtho
 140       FORMAT("NOfOrtho     = ", I10.1, ";")

           WRITE(IChList(ICh),150,ERR=20) NOfGamma
 150       FORMAT("NOfGamma     = ", I10.1, ";")

           WRITE(IChList(ICh),170,ERR=20) IKeepFlag
 170       FORMAT("IKeepFlag    = ", I10.1, ";")

           WRITE(IChList(ICh),180,ERR=20) IWriteFlag
 180       FORMAT("IWriteFlag   = ", I10.1, ";")

           WRITE(IChList(ICh),185,ERR=20) IMatWrite
 185       FORMAT("IMatWrite   = ", I10.1, ";")

           WRITE(IChList(ICh),190,ERR=20) ISortFlag
 190       FORMAT("ISortFlag    = ", I10.1, ";")

           WRITE(IChList(ICh),192,ERR=20) IBCXFlag
 192       FORMAT("IBCXFlag     = ", I10.1, ";")

           WRITE(IChList(ICh),194,ERR=20) IBCYFlag
 194       FORMAT("IBCYFlag     = ", I10.1, ";")

           WRITE(IChList(ICh),200,ERR=20) WidthX0
 200       FORMAT("WidthX0      = ", I10.1, ";")

           WRITE(IChList(ICh),210,ERR=20) WidthX1
 210       FORMAT("WidthX1      = ", I10.1, ";")

           WRITE(IChList(ICh),220,ERR=20) dWidthX
 220       FORMAT("dWidthX      = ", I10.1, ";")

           WRITE(IChList(ICh),222,ERR=20) WidthY
 222       FORMAT("WidthY       = ", I10.1, ";")

           WRITE(IChList(ICh),290,ERR=20) DiagDis0
 290       FORMAT("DiagDis0     = ", G18.9, ";")

           WRITE(IChList(ICh),293,ERR=20) DiagDis1
 293       FORMAT("DiagDis1     = ", G18.9, ";")

           WRITE(IChList(ICh),296,ERR=20) dDiagDis
 296       FORMAT("dDiagDis     = ", G18.9, ";")

           WRITE(IChList(ICh),300,ERR=20) Energy
 300       FORMAT("energy       = ", G18.9, ";")

           WRITE(IChList(ICh),310,ERR=20) Epsilon
 310       FORMAT("epsilon      = ", G18.9, ";")

           WRITE(IChList(ICh),400,ERR=20) IWidthX
 400       FORMAT("WidthX       = ", I10.1, ";")
 
           WRITE(IChList(ICh),402,ERR=20) IWidthY
 402       FORMAT("WidthY       = ", I10.1, ";")

           WRITE(IChList(ICh),500,ERR=20)
 500       FORMAT("data= {")

 1000   CONTINUE

        RETURN

C	error in OPEN detected
 10	PRINT*,"OpenOutputAvg(): ERR in OPEN()"
        IErr= 1
	RETURN

C	error in WRITE detected
 20	PRINT*,"OpenOutputAvg(): ERR in WRITE()"
	IErr= 1
	RETURN

        END


C	--------------------------------------------------------------------
C	WriteOutputAvg:
C
C	IErr	error code

	SUBROUTINE WriteOutputAvg(  
     +     IWidthX, IWidthY, DiagDis,
     +     gam, var, NOfL,
     +     psi,
     +     IErr )
        
        INCLUDE "common.h"

	INTEGER IWidthX, IWidthY, NOfL, IErr
        
        REAL*8 DiagDis, psi(IWidthX*IWidthY,IWidthX*IWidthY)

        REAL*8 gam(NOfL), var(NOfL) 

        INTEGER iState, jSite, iL

D	PRINT*,"DBG: WriteOutputAvg()"

	IErr= 0

C       average Lyapunov exponent Gamma

        DO 100 iL= 1,NOfL

           WRITE(IChOut,410,ERR=10) 
     +          iL, DiagDis, 
     +          gam(iL), var(iL)

 410       FORMAT("{ ", I7.1, ", ", G15.6, 
     +          ", ", G25.16, ", ", G25.16, " },")

 100    CONTINUE

C       wave functions

        IF( IWriteFlag.GE.2 ) THEN

           DO 1000 iState= 1, IWidthX*IWidthY

C              WRITE(IChOutPsi,510,ERR=10) 
C 510          FORMAT("{ ")

              DO 2000 jSite= 1, IWidthX*IWidthY
                 WRITE(IChOutPsi,550,ERR=10) iState,jSite,
     +              psi(jSite,iState)**2
 550             FORMAT( I4.1, " ", I4.1, " ", 1(G25.16) )
 2000         CONTINUE

C              WRITE(IChOutPsi,570,ERR=10) 
C 570          FORMAT("}, ")

 1000      CONTINUE

        ENDIF

        RETURN

C       ERR in Write detected
 10	PRINT*,"WriteOutputAvg(): ERR in WRITE()"
	IErr= 1
	RETURN

        END


C	--------------------------------------------------------------------
C	CloseOutputAvg:
C
C	IErr	error code

	SUBROUTINE CloseOutputAvg( IErr )

        INCLUDE "common.h"

        INTEGER IErr

        INTEGER ICh, IChList(MAXFiles)

D	PRINT*,"DBG: CloseOutputAvg()"

	IErr= 0

        IChList(1)= IChOut
C        IChList(2)= IChOutPsi

        DO 1000 ICh= 1, 1

           WRITE(IChList(ICh),100,ERR=20)
 100       FORMAT("}")

           CLOSE(UNIT= IChList(ICh), ERR= 10)

 1000   CONTINUE

        RETURN

C       ERR in CLOSE detected
 10	PRINT*,"CloseOutputAvg(): ERR in CLOSE()"
	IErr= 1
	RETURN

C	error in WRITE detected
 20	PRINT*,"CloseOutputAvg(): ERR in WRITE()"
	IErr= 1
	RETURN

        END


C	--------------------------------------------------------------------
C	OpenOutputGamma:
C
C	IErr	error code

	SUBROUTINE OpenOutputGamma( IWidthX, IWidthY, 
     +     DiagDis,
     +     IErr )

        INCLUDE "common.h"

        INTEGER IWidthX,IWidthY, IErr
        REAL*8 DiagDis

        INTEGER ICh

        CHARACTER*13 CName

D	PRINT*,"DBG: OpenOutputGamma()"

	IErr= 0
        ICh = IChOutGam

        WRITE(CName, '(I4.4,I4.4,A1,I4.4)') IWidthX,IWidthY,
     +        '.', NINT(100.0D0*ABS(DiagDis))

C       the filename is different for this gamma logging

        OPEN(UNIT= ICh, ERR= 10, STATUS= 'UNKNOWN', FILE=CName)

        RETURN

C	error in OPEN detected
 10	PRINT*,"OpenOutputGamma(): ERR in OPEN()"
	IErr= 1
	RETURN

        END


C	--------------------------------------------------------------------
C	WriteOutputGamma:
C
C	IErr	error code

	SUBROUTINE WriteOutputGamma( index, gam, var, NOfL, IErr )
        
        INCLUDE "common.h"

	INTEGER index, NOfL, IErr
        REAL*8 gam(NOfL), var(NOfL)

        INTEGER ICh, iL

D	PRINT*,"DBG: WriteOutputGamma()"

	IErr= 0
        ICh = IChOutGam

C       Lyapunov exponent Gamma

        DO 100 iL=1,NOfL
           WRITE(ICh,410,ERR=10) index, iL, gam(iL), var(iL)
 410       FORMAT(" ",I7.1, " ", I4.1, " ", 
     +        G25.16, " ", G25.16)
 100    CONTINUE

        RETURN

C       ERR in Write detected
 10	PRINT*,"WriteOutputGamma(): ERR in WRITE()"
	IErr= 1
	RETURN

        END


C	--------------------------------------------------------------------
C	CloseOutputGamma:
C
C	IErr	error code

	SUBROUTINE CloseOutputGamma( IErr )

        INCLUDE "common.h"

        INTEGER IErr

        INTEGER ICh

D	PRINT*,"DBG: CloseOutputGamma()"

	IErr= 0
        ICh = IChOutGam

        CLOSE(UNIT= ICh, ERR= 10)

        RETURN

C       ERR in CLOSE detected
 10	PRINT*,"CloseOutputGamma(): ERR in CLOSE()"
	IErr= 1
	RETURN

        END


c --------------------------------------------------------------------

	subroutine matwrit(a,n,iflag)

	integer n,i,j,iflag
	real*8 a(n,n)

	write(91,*) ' '

	if (iflag.eq.0) then
	do i=1,n
	  write(91,*) 'row ',i
	  write(91,'(20f4.0)') (a(i,j),j=1,n)
	  write(91,*) ' '
	enddo
	endif
	if (iflag.eq.1) then
	do i=1,n
	  write(91,*) 'row ',i
	  write(91,'(10g18.8)') (a(i,j),j=1,n)
	  write(91,*) ' '
	enddo
	endif
	if (iflag.eq.2) then
	do i=1,n
	  write(91,*) 'row ',i
	  write(91,'(10g14.6)') (a(i,j),j=1,n)
	  write(91,*) ' '
	enddo
	endif
	write(91,*) ' '
	return
	end
















