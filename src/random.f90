!********************************************************************
! $Header: /home/cvs/phsht/GrapheneTMM/src/Restart/ALL/random.f90,v 1.1 2011/07/22 17:49:19 ccspam Exp $
!********************************************************************

!********************************************************************
! $Log: random.f90,v $
! Revision 1.1  2011/07/22 17:49:19  ccspam
! Programs for ZZ and AC with the restart routine for francesca
!
! Revision 1.1  2010/04/19 15:20:46  phrial
! *** empty log message ***
!
! Revision 1.2  2009/06/09 15:59:23  phrial
! New modular form of code.
!
!********************************************************************

module RNG
  
  USE MyNumbers
  USE Randoms

  implicit none
  
  !accessibility
  !private
  !public :: SETVALUES
  !public :: DRANDOM
  
  integer(kind=ikind) :: b, k
  real(kind=rkind) :: dummy, dret
  real(kind=rkind),    parameter :: AM=4.656612873077d-10
  integer(kind=ikind), parameter :: IA=16807, IM=2147483647, IQ=127773, IR=2836
  
contains
  subroutine SETVALUES(ISeed)

    integer(kind=ikind) :: ISeed
    integer(kind=ikind) :: idum
    
    idum=ISeed

    if (idum.le.0) idum=1
    k=(idum)/IQ
    idum=IA*(idum-k*IQ)-IR*k
    if (idum.lt.0) idum = idum + IM
    if (idum.lt.2) then
      z1=idum+2 
    else 
      z1=idum
    endif
    k=(idum)/IQ
    idum=IA*(idum-k*IQ)-IR*k    
    if (idum.lt.0) idum = idum + IM
    if (idum.lt.8) then 
      z2=idum+8 
    else 
      z2=idum
    endif
    k=(idum)/IQ
    idum=IA*(idum-k*IQ)-IR*k
    if (idum.lt.0) idum = idum + IM
    if (idum.lt.16) then
      z3=idum+16 
    else 
      z3=idum
    endif
    k=(idum)/IQ
    idum=IA*(idum-k*IQ)-IR*k
    if (idum.lt.0) idum = idum + IM
    if (idum.lt.128) then
      z4=idum+128 
    else 
      z4=idum
    endif

    dummy=DRANDOM()
  end subroutine SETVALUES

  real function DRANDOM() result(dRet)

    b  = ishft(ieor(ishft(z1,6),z1),-13)
    z1 = ieor(ishft(iand(z1,-2),18),b)
    
    b  = ishft(ieor(ishft(z2,2),z2),-27)
    z2 = ieor(ishft(iand(z2,-8),2),b)
    
    b  = ishft(ieor(ishft(z3,13),z3),-21)
    z3 = ieor(ishft(iand(z3,-16),7),b)
    
    b  = ishft(ieor(ishft(z4,3),z4),-12)
    z4 = ieor(ishft(iand(z4,-128),13),b)
    
    dRet=ishft(ieor(ieor(ieor(z1,z2),z3),z4),-1)*AM

  end function DRANDOM
end module RNG

