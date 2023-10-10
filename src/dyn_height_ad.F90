SUBROUTINE DYN_HEIGHT_AD

!-----------------------------------------------------------------------
!                                                                      !
! CALCULATE VERTICAL INTEGRAL OF BOUYANCY GRADIENT                     !
!                                                                      !
! VERSION 1: S.DOBRICIC 2007                                           !
!-----------------------------------------------------------------------


  USE SET_KND
  USE GRD_STR
  USE BAL
  USE LBCLNK
  USE RUN, ONLY : LL_SSH_UNBALANCED

  IMPLICIT NONE

  INTEGER(I4)    :: K

  GRD%DNS = 0._R8
  DO K=NLEVS,1,-1
     GRD%DNS(:,:,K) = - DHDZ(K)*GRD%ETA_AD/RHO0 
  ENDDO

END SUBROUTINE DYN_HEIGHT_AD
