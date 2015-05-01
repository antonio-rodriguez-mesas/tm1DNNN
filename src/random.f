C       ********************************************************************
C       
C       RANDOM - Standard F77/F90 interface for random number generators
C
C       ********************************************************************
       
C       ********************************************************************
C       
C       $Header: /home/kerekou1/rar/f77/tmODse2d/RCS/random.f,v 1.4 1996/12/10 12:40:41 rar Exp $
C
C       ********************************************************************

C **************************************************************************
C
C $Log: random.f,v $
c Revision 1.4  1996/12/10  12:40:41  rar
c DBLE->REAL*8
c
c Revision 1.3  1996/11/14  11:14:15  rar
c everything is capitalized to be compatible with the HP Compiler;
c
c Revision 1.1  1996/05/21  10:09:02  rar
c Initial revision
c
c Revision 3.1  96/04/16  11:35:37  11:35:37  rar (Rudolf Roemer)
c initial RCS version of V2.5
c 
C
C **************************************************************************


C       ********************************************************************
C
C       History prior to version 2.0:
C
C       24/01/96 RAR: perfectly new installation
C
C       ********************************************************************

C	--------------------------------------------------------------------
C       SRANDOM: 
C
C       Random number generator SEED interface for use with any old RND

        SUBROUTINE SRANDOM( ISeed )

        INTEGER ISeed, IDUM
        REAL*8 Dummy

        COMMON/IFIXED/IDUM

C       change these lines to incorporate different RND generators
        REAL*8 RAN2
        EXTERNAL RAN2

        IDUM = ISeed
        Dummy= RAN2(-IDUM)

        RETURN
        END

C	--------------------------------------------------------------------
C       DRANDOM: 
C
C       Random number generator interface for use with any old RND
C
C       NOTE that ISeed is never used!

        REAL*8 FUNCTION DRANDOM( ISeed )

        INTEGER ISeed, IDUM
        COMMON/IFIXED/IDUM

C       change these lines to incorporate different RND generators
        REAL*8 RAN2
        EXTERNAL RAN2

        DRANDOM= RAN2(IDUM)

        RETURN
        END

C	--------------------------------------------------------------------
C       RAN2: 
C
C       Random number generator of Numerical Recipes,
C       (C) Copr. 1986-92 Numerical Recipes Software 6?6>)AY.

        REAL*8 FUNCTION RAN2(IDUM)

        IMPLICIT REAL*8 (A-H,O-Z)

        INTEGER IDUM,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
        REAL AM,EPS,RNMX
        PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *       IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *       NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2D-7,RNMX=1.-EPS)
        INTEGER IDUM2,J,K,IV(NTAB),IY
        SAVE IV,IY,IDUM2
        DATA IDUM2/123456789/, IV/NTAB*0/, IY/0/

        IF (IDUM.LE.0) THEN
           IDUM=MAX(-IDUM,1)
           IDUM2=IDUM
           DO 11 J=NTAB+8,1,-1
              K=IDUM/IQ1
              IDUM=IA1*(IDUM-K*IQ1)-K*IR1
              IF (IDUM.LT.0) IDUM=IDUM+IM1
              IF (J.LE.NTAB) IV(J)=IDUM
 11        CONTINUE
           IY=IV(1)
        ENDIF
        K=IDUM/IQ1
        IDUM=IA1*(IDUM-K*IQ1)-K*IR1
        IF (IDUM.LT.0) IDUM=IDUM+IM1
        K=IDUM2/IQ2
        IDUM2=IA2*(IDUM2-K*IQ2)-K*IR2
        IF (IDUM2.LT.0) IDUM2=IDUM2+IM2
        J=1+IY/NDIV
        IY=IV(J)-IDUM2
        IV(J)=IDUM
        IF(IY.LT.1)IY=IY+IMM1
        RAN2=MIN(AM*IY,RNMX)
        RETURN

        END

