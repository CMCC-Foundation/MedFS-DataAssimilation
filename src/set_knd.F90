MODULE SET_KND

!-----------------------------------------------------------------------
!                                                                      !
! THE PRECISION OF REALS AND INTEGERS                                  !
!                                                                      !
! VERSION 1: S.DOBRICIC 2006                                           !
!-----------------------------------------------------------------------


IMPLICIT NONE

PUBLIC

   INTEGER, PARAMETER ::                &
      DP = SELECTED_REAL_KIND(12,307),  &  ! REAL*8
      R4 = SELECTED_REAL_KIND( 6, 37),  &  ! REAL*4
#ifdef NO_DOUBLEPREC
      R8 = R4
#else
      R8 = DP
#endif

   INTEGER, PARAMETER ::                &
      I4 = SELECTED_INT_KIND(9) ,       &  ! INTEGER*4
      I6 = SELECTED_INT_KIND(12) ,      &  !
      I8 = SELECTED_INT_KIND(14)           ! INTEGER*8

END MODULE SET_KND
