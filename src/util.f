C       ********************************************************************
C       
C       TMSE3DAO - Transfer matrix method for the 3D Anderson model with
C                  automatic reorthogonalization
C
C       ********************************************************************
       
C       ********************************************************************
C       
C       $Header: /home/kerekou1/rar/f77/tmse3dAO/RCS/util.f,v 2.2 1998/06/04 07:05:07 rar Exp $
C
C       ********************************************************************

C **************************************************************************
C
C $Log: util.f,v $
c Revision 2.2  1998/06/04  07:05:07  rar
c hardly any change
c
c Revision 2.1  1998/04/29  13:50:38  rar
c BLAS version
c
c Revision 1.6  1998/04/03  13:49:29  rar
c swapped PSI_A/B(x,y) -> PSI_A/B(y,x)
c
c Revision 1.5  1998/04/02  08:40:24  rar
c removd commented old tmOD lines
c
c Revision 1.4  1998/04/01  16:01:14  crv
c removed DisCenter, DisWidth
c
c Revision 1.3  1998/04/01  09:45:02  crv
c removed PresentVec, PastVec and HopFuture
c
c Revision 1.1  1998/03/30  14:17:48  rar
c Initial revision
c
c Revision 1.7  1996/12/10  12:40:02  rar
c DBLE->REAL*8
c
c Revision 1.6  1996/12/04  11:20:35  rar
c included a diagonal disorder potential (DiagDis0,DiagDis1,dDiagDis)
c
c Revision 1.5  1996/11/14  11:13:51  rar
c TMMult() receives a TINY to check for division by a zero random number;
c
c Revision 1.4  96/11/13  10:38:23  10:38:23  rar (Rudolf Roemer)
c F90-linux compatible;
c recorrected a mistake in the coding of the hopping,
c replaced for PsiLeft with the correct version
c 
c Revision 1.3  1996/11/04  11:42:40  rar
c removed pLevel
c
c Revision 1.2  96/11/04  11:30:06  11:30:06  rar (Rudolf Roemer)
c replaced the wrong PastVec with the right PresentVec in the BC-IF's
c also replaced a wrong 1 in the BC-IF's
c 
c Revision 1.1  96/11/01  16:40:25  16:40:25  rar (Rudolf Roemer)
c Initial revision
c 
C
C **************************************************************************

C	--------------------------------------------------------------------
C       TMMult:
C
C       Multiplication of the transfer matrix onto the vector (PSI_A,PSI_B),
C       giving (PSI_B,PSI_A) so that the structure of the transfer matrix 
C       can be exploited

        SUBROUTINE TMMult(PSI_A,PSI_B,Righth,
     +                    En,DiagDis, MX,MY, IBCXFlag, IBCYFlag, coninv)

C       wave functions:
C
C       (PSI_A, PSI_B) on input, (PSI_B,PSI_A) on output


        INTEGER MX,MY,                    ! plane width
     +          IBCXFlag, IBCYFlag        ! boundary condition flags
                                          ! along the X and Y axes

        REAL*8            DiagDis,        ! diagonal disorder
     +                    En              ! energy

        REAL*8 PSI_A(MX,MY,MX*MY), PSI_B(MX,MY,MX*MY)
	real*8 Righth(MX,MY,MX*MY),coninv(MX*MY,MX*MY)

        INTEGER iSite,jSite,iState,n
        REAL*8  new, OnsitePot
	real*8  PsiPast,PsiPastLeft,PsiPastRight
	real*8  PsiUpRight,PsiRight,PsiDownRight
	real*8  PsiDownLeft,PsiLeft,PsiUpLeft


        REAL*8 DRANDOM
        EXTERNAL DRANDOM

D	PRINT*,"DBG: TMMult()"

        DO 100 iSite=1,MX

           DO 120 jSite=1,MY

C             create the new onsite potential

              OnsitePot= En - (DRANDOM(iSite)-0.5d0)*DiagDis

              DO 200 iState=1,MX*MY

CD             PRINT*,"iState, iSite, jSite",
CD    +            iState,iSite,jSite


	      PsiPast=PSI_B(iSite,jSite,iState)

	      if ((IBCXFlag.eq.0.or.IBCYFlag.eq.0).and.
     +	      (iSite.eq.1.or.iSite.eq.MX.or.jSite.eq.1.or.jSite.eq.MY)) then
	       PsiPastLeft=0.d0
	       PsiPastRight=0.d0
	       PsiUpRight=0.d0
	       PsiRight=0.d0
	       PsiDownRight=0.0d0
	       PsiDownLeft=0.d0
	       PsiLeft=0.d0
	       PsiUpLeft=0.d0
	      endif           ! if we are at the boundary we set these variables
	                      ! to 0, thus we don't have to check hard wall
			      ! boundary condition

c             setting PsiPastLeft

	      if (iSite.eq.1) then
	        if (IBCXFlag.eq.1) then
	           if (jSite.eq.1) then
		     if (IBCYFlag.eq.1) PsiPastLeft=PSI_B(MX,MY,iState)
		   else
		     PsiPastLeft=PSI_B(MX,jSite-1,iState)
		   endif
		 endif
	      else
		   PsiPastLeft=PSI_B(iSite-1,jSite,iState)
              endif


c            setting PsiPastRight

	      if (iSite.eq.1) then
	        if (IBCXFlag.eq.1) then
		     PsiPastRight=PSI_B(MX,jSite,iState)
                endif
	      else
	        if (jSite.eq.MY) then
		   if (IBCYFlag.eq.1) PsiPastRight=PSI_B(iSite-1,1,iState)
		else
		   PsiPastRight=PSI_B(iSite-1,jSite+1,iState)
		endif
              endif


c             setting PsiUpRight

	      if (iSite.eq.1) then
	         if (IBCXFlag.eq.1) then
		    PsiUpRight=PSI_A(MX,jSite,iState)
		 endif
	      else
		 if (jSite.eq.MY) then
		    if (IBCYFlag.eq.1) PsiUpRight=PSI_A(iSite-1,1,iState)
		 else
		    PsiUpRight=PSI_A(iSite-1,jSite+1,iState)
		 endif
	      endif


c             setting PsiRight

	      if (jSite.eq.MY) then
		 if (IBCYFlag.eq.1) PsiRight=PSI_A(iSite,1,iState)
	      else
	         PsiRight=PSI_A(iSite,jSite+1,iState)
	      endif


c            setting PsiDownRight

	     if (iSite.eq.MX) then
		 if (IBCXFlag.eq.1) then
	           if (jSite.eq.MY) then
		     if (IBCYFlag.eq.1) PsiDownRight=PSI_A(1,1,iState)
		   else
		     PsiDownRight=PSI_A(1,jSite+1,iState)
		   endif
		 endif
	     else
	         PsiDownRight=PSI_A(iSite+1,jSite,iState)
	     endif


c            setting PsiDownLeft

	     if (iSite.eq.MX) then
	         if (IBCXFlag.eq.1) PsiDownLeft=PSI_A(1,jSite,iState)
	     else
	         if (jSite.eq.1) then
		    if (IBCYFlag.eq.1) PsiDownLeft=PSI_A(iSite+1,MY,iState)
		 else
		    PsiDownLeft=PSI_A(iSite+1,jSite-1,iState)
		 endif
	     endif


c            setting PsiLeft

             if (jSite.eq.1) then
	       if (IBCYFlag.eq.1) PsiLeft=PSI_A(iSite,MY,iState)
	     else
	       PsiLeft=PSI_A(iSite,jSite-1,iState)
	     endif


c            setting PsiUpLeft

	     if (iSite.eq.1) then
		if (IBCXFlag.eq.1) then
		   if (jSite.eq.1) then
		     if (IBCYFlag.eq.1) PsiUpLeft=PSI_A(MX,MY,iState)
		   else
		     PsiUpLeft=PSI_A(MX,jSite-1,iState)
		   endif
		endif
	     else
	        PsiUpLeft=PSI_A(iSite-1,jSite,iState)
	     endif



                 new=  OnsitePot * PSI_A(iSite,jSite,iState)
     +           - PsiUpRight - PsiRight - PsiDownRight
     +           - PsiDownLeft - PsiLeft - PsiUpLeft
     +           - PsiPast - PsiPastLeft - PsiPastRight

CD                  PRINT*,"En, PL, PR, PU, PD, PA,PB, PN"
CD                  PRINT*, En, PsiLeft, PsiRight,PsiUp,PsiDown
CD    +             PSI_A(iState,iSite,jSite), PSI_B(iState,iSite,jSite),
CD    +             new

                 Righth(iSite,jSite,iState)= new
 200          CONTINUE

 120       CONTINUE
 100    CONTINUE

 	   n=MX*MY

	   call dgemm ('n','n', n, n, n, 1.0d0, coninv, n, Righth,
     $	               n, 0.d0, PSI_B, n)



        RETURN
        END

C	--------------------------------------------------------------------
C       ReOrtho:
C
C       Gram-Schmidt orthonormalization of the wave vectors (PSI_A,PSI_B),
C       see e.g. Horn/Johnson, "Matrix Analysis", pp. 15.
C
C       PLUS reading off the Luapunov exponents
C
C       (PSI_A,PSI_B) is the incoming vector of eigenvectors,
C       GAMMA, GAMMA2 are Luapunov exponent and its square,
C       M is the width of the strip,
C
C       IVec/JVec label different vectors, 
C       KIndex is the component index of a SINGLE vector,
C
C       dummy, sum and norm are just different names for the same
C       double precision workspace

      SUBROUTINE ReOrtho(PSI_A,PSI_B,GAMMA,GAMMA2,MX,MY,NOrtho)

      INCLUDE "common.h"

      INTEGER MX,MY,NOrtho

      REAL*8 PSI_A(MX*MY,MX*MY), PSI_B(MX*MY,MX*MY)
      REAL*8 GAMMA(MX*MY), GAMMA2(MX*MY)
        
      INTEGER IVec,JVec,KIndex, MXY

      REAL*8 dummy,sum,norm,normbefore,quot
      EQUIVALENCE (dummy,sum)

      REAL*8 MinAccuracy,MaxAccuracy
      PARAMETER (MinAccuracy= 1D-11, MaxAccuracy= 1D-8)

C     BLAS routines

      REAL*8 ddot
      EXTERNAL ddot

C       make the local variables static
      SAVE

D     PRINT*,"DBG: ReOrtho()"
D     PRINT*,"DBG: MinAccuracy=", MinAccuracy,
D    +          ", MaxAccuray=", MaxAccuracy
        
      MXY       = MX*MY

      normbefore= 
     +   ddot(MXY, PSI_A(1,MXY),1, PSI_A(1,MXY),1) +
     +   ddot(MXY, PSI_B(1,MXY),1, PSI_B(1,MXY),1)

      DO 100 IVec=1,MXY

           DO 200 JVec=1,IVec-1

              sum= -(
     +           ddot(MXY, PSI_A(1,JVec),1, PSI_A(1,IVec),1) +
     +           ddot(MXY, PSI_B(1,JVec),1, PSI_B(1,IVec),1)
     +              )

              call daxpy(MXY, sum, PSI_A(1,JVec),1, PSI_A(1,IVec),1)
              call daxpy(MXY, sum, PSI_B(1,JVec),1, PSI_B(1,IVec),1)

 200       CONTINUE

C          calculation of norm

           norm= ddot(MXY, PSI_A(1,IVec),1, PSI_A(1,IVec),1) +
     +           ddot(MXY, PSI_B(1,IVec),1, PSI_B(1,IVec),1)
 
 500       CONTINUE   

           dummy= 1.D0/SQRT(norm)

           call dscal(MXY, dummy, PSI_A(1,IVec), 1)
           call dscal(MXY, dummy, PSI_B(1,IVec), 1)

C          ----------------------------------------------------------------
C          gammadummy is ordered s.t. the vector with the smallest overall 
C          norm is last and the vector with the largest overall norm first. 
C          We sum this so that the smallest value is used for GAMMA(M)
C          and the largest for GAMMA(1). Same for GAMMA2(). Therefore, the
C          inverse of the smallest Lyapunov exponent is obtained by taking
C          LAMBDA(1) = N/GAMMA(M)
C           PRINT *,"DBG: dummy,LOG(dummy)", dummy,LOG(dummy)

           dummy       = LOG(dummy)

           GAMMA(IVec) = GAMMA(IVec) - dummy
           GAMMA2(IVec)= GAMMA2(IVec) + dummy*dummy

C           IF(IVec.EQ.MX*MY) PRINT*, dummy

 100    CONTINUE
        
C       determine whether NOrtho needs to be changed

        quot = SQRT(norm/normbefore)

C       decrease NOrtho if accuracy is bad
        IF (quot.LT.MinAccuracy) THEN
           NOrtho = NOrtho-2
           IF (NOrtho.LT.1) THEN
              NOrtho = 1
           ENDIF
        ENDIF

C       increase NOrtho if accuracy is good
        IF (quot.GT.MaxAccuracy) THEN
           NOrtho = NOrtho+2
        ENDIF

        RETURN
        END


C	--------------------------------------------------------------------
C       Swap:
C
C       (PSI_A,PSI_B)= (old,new) is the incoming vector, this is swapped
C       into (PSI_A,PSI_B)= (new,old)

        SUBROUTINE Swap( PSI_A, PSI_B, MX, MY)

        INTEGER MX,MY
        REAL*8 PSI_A(MX*MY,MX*MY), PSI_B(MX*MY,MX*MY)

        INTEGER jState, index
        REAL*8 dummy

D	PRINT*,"DBG: Swap()"

        DO 100 jState=1,MX*MY
           DO 200 index=1,MX*MY

              dummy              = PSI_B(index,jState)
              PSI_B(index,jState)= PSI_A(index,jState)
              PSI_A(index,jState)= dummy

 200       CONTINUE
 100    CONTINUE

        RETURN
        END


C	--------------------------------------------------------------------
C       ReSort:
C
C       sort the Lyapunov eigenvalues s.t. the largest comes first. RESORT()
C       is ShellSort taken from NumRec, SHELL().

        SUBROUTINE ReSort( PSI_A, PSI_B, array0, array1, MX,MY )

        INTEGER MX,MY,N
        REAL*8 PSI_A(MX*MY,MX*MY),PSI_B(MX*MY,MX*MY), 
     +         array0(MX*MY), array1(MX*MY)

        REAL*8 ALN2I, TINY
        PARAMETER (ALN2I=1.4426950D0, TINY=1.D-5)

        INTEGER NN,M,L,K,J,I,LOGNB2, index
        REAL*8 dummyA, dummyB

D	PRINT*,"DBG: ReSort()"

CD       PRINT*,"array0(1),array0(N)",array0(1),array0(N)

        
        N=MX*MY
        LOGNB2=INT(LOG(REAL(N))*ALN2I+TINY)
        M=N
        DO 12 NN=1,LOGNB2
           M=M/2
           K=N-M
           DO 11 J=1,K
              I=J
 3            CONTINUE
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

                    PSI_A(index,I)= PSI_A(index,I)
                    PSI_B(index,I)= PSI_B(index,I)

                    PSI_A(index,L)= dummyA
                    PSI_B(index,L)= dummyB
 100             CONTINUE

                 I=I-M
                 IF(I.GE.1) GOTO 3
              ENDIF
 11        CONTINUE
 12     CONTINUE

CD       PRINT*,"array0(1),array0(N)",array0(1),array0(N)
        RETURN
        END


c ------------------------------------------------------------------------

	subroutine cmatinver(MX,MY,con,coninv,righth,af,
     $                    berr,ferr,c,r,work,ipiv,iwork,IBCX,
     $                    IBCY,IMatWrite)

     	integer MX,MY,IBCX,IBCY,IMatWrite
	real*8 con(MX,MY,MX,MY),coninv(MX,MY,MX,MY)
	real*8 righth(MX,MY,MX,MY),af(MX*MY,MX*MY)
	real*8 berr(MX*MY),ferr(MX*MY),c(MX*MY),r(MX,MY)
	real*8 work(3*MX*MY),rcond

	integer ipiv(MX*MY),iwork(MX*MY),info
	character*1 equed,fact,trans,transa,transb
	character*12 name

	real*8  dalpha,dbeta

	integer irow,jcolumn,n

	do irow=1,MX
	 do jcolumn=1,MY
	  con(irow,jcolumn,irow,jcolumn)=1.0d0

	  if (irow.eq.MX) then
	     if (IBCX.eq.1) then
	       con(irow,jcolumn,1,jcolumn)=1.0d0
	       if (jcolumn.eq.MY) then
	         if (IBCY.eq.1) con(irow,jcolumn,1,1)=1.0d0
	       else
	         con(irow,jcolumn,1,jcolumn+1)=1.0d0
	       endif
	     endif
	  else
	     con(irow,jcolumn,irow+1,jcolumn)=1.0d0
	     if (jcolumn.eq.1) then
	       if (IBCY.eq.1) con(irow,jcolumn,irow+1,MY)=1.0d0
	     else
	       con(irow,jcolumn,irow+1,jcolumn-1)=1.0d0
	     endif
	  endif

	 enddo
	enddo

	do irow=1,MX
	  do jcolumn=1,MY
	   righth(irow,jcolumn,irow,jcolumn)=1.0d0
	  enddo
	enddo

	fact='e'
	trans='n'

	n=MX*MY

	if (IMatWrite.ne.0) then
	  name(9:12)='.mat'
	  write(name(1:8),'(i4.4,i4.4)') MX,MY
	  open(91,file=name)
	  write(91,*) ' '
	  write(91,*) 'conmat'
	  call matwrit(con,n,1)
	endif


	call dgesvx(fact,trans,n,n,con,n,af,n,ipiv,equed,r,c,righth,
     $              n,coninv,n,rcond,ferr,berr,work,iwork,info)


        if (IMatWrite.ne.0) then
	   write(91,*) 'info= ',info
	   write(91,*) 'rcond= ',rcond
	   write(91,*) 'equed= ',equed
	   write(91,*) ' '

	   if (info.eq.0) then
	      write(91,*) 'conmat(-1)'
	      call matwrit(coninv,n,2)

	      transa='n'
	      transb='n'
	      dalpha=1.0d0
	      dbeta=0.0d0

	   call dgemm ('n','n', n, n, n, 1.0d0, con, n, coninv,
     $	               n, 0.d0, righth, n)

	      write(91,*) 'con*conmat(-1)'
	      call matwrit(righth,n,2)
	      close(91)
	    endif
	 endif


	 if (info.ne.0) then
	   write(*,*) 'singular matrix'
	   write(*,*) 'info= ',info,' rcond= ',rcond
	   stop
	 endif

	 return
	 end


