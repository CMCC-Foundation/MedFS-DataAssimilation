SUBROUTINE OBS_SLA

!-----------------------------------------------------------------------
!                                                                      !
! APPLY OBSERVATIONAL OPERATOR FOR SLA                                 !
!                                                                      !
! VERSION 1: S.DOBRICIC 2006                                           !
!-----------------------------------------------------------------------


 USE SET_KND
 USE GRD_STR
 USE OBS_STR

 IMPLICIT NONE

 INTEGER(I4)   ::  K

 DO K=1,SLA%NO

  IF(SLA%FLC(K).EQ.1)THEN

    SLA%INC(K) = &
    & OSUM( SLA%PQ(K,1:NPQ),GRD%ETA(SLA%IB(K,1:NPQ),SLA%JB(K,1:NPQ)) )

  ENDIF

 ENDDO

END SUBROUTINE OBS_SLA
