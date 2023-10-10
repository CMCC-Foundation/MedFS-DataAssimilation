SUBROUTINE INT_PAR

!-----------------------------------------------------------------------
!                                                                      !
! CALCULATE INTERPOLATION PARAMETERS                                   !
!                                                                      !
! VERSION 1: S.DOBRICIC 2006                                           !
!-----------------------------------------------------------------------


 USE SET_KND
 USE OBS_STR

 IMPLICIT NONE

  INTEGER(I4)    ::  K, I

! ----
! LOAD SLA OBSERVATIONS
  WRITE(*,*) 'INTPAR CALLING INTAPARSLA'
  CALL INT_PAR_SLA
  WRITE(*,*) 'INTPAR CALLING INTABACVAL'
  CALL INT_BACVAL_SLA
  WRITE(*,*) 'CALLING 2ND OBSOPER'

! LOAD ARGO OBSERVATIONS
  CALL INT_PAR_ARG

! LOAD XBT OBSERVATIONS
  CALL INT_PAR_XBT

! LOAD GLIDER OBSERVATIONS
  CALL INT_PAR_GLD

! LOAD VELOCITY OBSERVATIONS
  CALL INT_PAR_VEL

END SUBROUTINE INT_PAR
