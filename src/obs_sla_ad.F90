SUBROUTINE OBS_SLA_AD

!-----------------------------------------------------------------------
!                                                                      !
! APPLY OBSERVATIONAL OPERATOR FOR SLA - ADJOINT                       !
!                                                                      !
! VERSION 1: S.DOBRICIC 2006                                           !
!-----------------------------------------------------------------------


 USE SET_KND
 USE GRD_STR
 USE OBS_STR

 IMPLICIT NONE

 INTEGER(I4)   ::  JP, K
 INTEGER(I4)   ::  XIND1

 DO K=1,SLA%NO

  IF(SLA%FLC(K).EQ.1)THEN

    OBS%K = OBS%K + 1

    DO JP=1,NPQ
      GRD%ETA_AD(SLA%IB(K,JP),SLA%JB(K,JP)) = &
      & GRD%ETA_AD(SLA%IB(K,JP),SLA%JB(K,JP)) + SLA%PQ(K,JP) *OBS%GRA(OBS%K)
    ENDDO

  ENDIF

 ENDDO


END SUBROUTINE OBS_SLA_AD
