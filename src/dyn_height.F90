SUBROUTINE DYN_HEIGHT

!-----------------------------------------------------------------------
!                                                                      !
! DYNAMIC HEIGHT FORMULATION
!                                                                      !
! A.S.
!-----------------------------------------------------------------------


  USE SET_KND
  USE GRD_STR
  USE BAL
  USE RUN, ONLY : LL_SSH_UNBALANCED
  USE LBCLNK
  USE OCEANTOOLS

  IMPLICIT NONE

  INTEGER(I4)    :: K
  REAL   (R8)    :: RHTL(GRD%IM,GRD%JM)

  RHTL = 0._R8
  IF( .NOT. LL_SSH_UNBALANCED ) GRD%ETA = 0._R8
  DO K=NLEVS,1,-1
      RHTL = RHTL  + GRD%DNS(:,:,K)*DHDZ(K)
  ENDDO
  GRD%ETA = GRD%ETA -1._R8 * RHTL/RHO0

END SUBROUTINE DYN_HEIGHT
