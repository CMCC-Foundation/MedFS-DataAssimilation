SUBROUTINE CHECK_SLABG(KITER)

!-----------------------------------------------------------------------
!                                                                      !
! APPLY OBSERVATIONAL OPERATOR FOR SLA                                 !
!                                                                      !
! VERSION 1: S.DOBRICIC 2006                                           !
!-----------------------------------------------------------------------


 USE SET_KND
 USE GRD_STR
 USE OBS_STR
 USE IOUNITS

 IMPLICIT NONE

 INTEGER(I4), INTENT(IN) ::  KITER
 INTEGER(I4)   ::  K, JVL, KBOT
 LOGICAL :: LLOK
 REAL(R8) :: TB, SB

#include "obs_events.h"

IF(SLA%NO.LE.0) RETURN

#ifdef SHARED_MEMORY
!$OMP PARALLEL DEFAULT(SHARED), PRIVATE(K,JVL,TB,SB,KBOT,LLOK)
!$OMP DO SCHEDULE(STATIC,1)
#endif
 CYOBS : DO K=1,SLA%NO

    IF( SLA%FLC(K) .NE. 1 ) CYCLE CYOBS

    KBOT=SLA%BOT(K)

    LLOK=.TRUE.
    DO JVL = KBOT,1,-1

       TB      = SLA%TB(K,JVL) 
       SB      = SLA%SB(K,JVL)

       IF(TB .LT. -10._R8 .OR. TB .GT. 50._R8 .OR. &
       SB .LT. 0._R8 .OR. SB .GT. 50._R8) LLOK=.FALSE.

    ENDDO

    IF(.NOT. LLOK) THEN
       SLA%FLC(K) = 0
       SLA%EVE(K) = keve_BGCH
    ENDIF

 ENDDO CYOBS
#ifdef SHARED_MEMORY
!$OMP END DO
!$OMP END PARALLEL
#endif

WRITE(IOUNLOG,*) ' ITER ',KITER
WRITE(IOUNLOG,*) ' CHECK SLA BG, OBS WRONG :',&
                & COUNT(SLA%EVE(1:SLA%NO).EQ.keve_BGCH)
CALL FLUSH(IOUNLOG)

END SUBROUTINE CHECK_SLABG
