SUBROUTINE OBS_SLA2

!-----------------------------------------------------------------------
!                                                                      !
! APPLY OBSERVATIONAL OPERATOR FOR SLA                                 !
!                                                                      !
! VERSION 1: S.DOBRICIC 2006                                           !
!-----------------------------------------------------------------------


 USE SET_KND
 USE GRD_STR
 USE OBS_STR
 USE RUN, ONLY : RTIMES1950, KLEVNM, RHO0,NTSTEPS, LLVARMDT, NCONF
 USE OCEANTOOLS , ONLY : RHO_UNESCOTL
 USE MYFRTPROF, ONLY : MYFRTPROF_WALL
 USE TLAD_VARS

 IMPLICIT NONE

 INTEGER(I4)   ::  I, J, K, KSTEP, JVL, JTSTEP, I2,KBOT
 INTEGER(I4)   ::  XIND1

 REAL(R8) :: ZC, T, S, TB, SB, RHTL, RTDIFF,ZRHO
 REAL(R8) ::  A2,A2_TL,B2,B2_TL,C2,ROOTS,RHOW_TL

 CALL MYFRTPROF_WALL('OBS_SLA2: SLA OBSERVATION OPERATOR',0)

#ifdef SHARED_MEMORY
!$OMP PARALLEL DEFAULT(SHARED), PRIVATE(K,JVL,KSTEP,RTDIFF,JTSTEP,&
!$OMP & RHTL,T,S,TB,SB,ZRHO,KBOT,A2,A2_TL,B2,B2_TL,C2,ROOTS,RHOW_TL)
!$OMP DO SCHEDULE(DYNAMIC,1)
#endif
CYOBS : DO K=1,SLA%NO

  IF( .NOT. LL_4DVAR .OR. ( LL_4DVAR .AND. SLA%TATS(K) .EQ. TLAD_TS )) THEN

  IF(NCONF.NE.202) THEN

  KBOT=SLA%BOT(K)

    RHTL = 0._R8
    DO JVL = KBOT,1,-1

       T = &
 & OSUM( SLA%PQ(K,1:NPQ),GRD%TEM(SLA%IB(K,1:NPQ),SLA%JB(K,1:NPQ),JVL) )
       S = &
 & OSUM( SLA%PQ(K,1:NPQ),GRD%SAL(SLA%IB(K,1:NPQ),SLA%JB(K,1:NPQ),JVL) )

        TB      = SLA%TB(K,JVL) 
        SB      = SLA%SB(K,JVL)

!       ZRHO = RHO_UNESCOTL(SB,TB,S,T,0._R8,.FALSE.)

#include "rho_unescotl_r8.h"

        RHTL = RHTL + ZRHO*GRD%DZ(JVL)

    ENDDO

    SLA%INC(K) = -1._R8*RHTL/RHO0

  ENDIF  

  IF(LLVARMDT) THEN
      SLA%INC(K) = SLA%INC(K) + &
      & OSUM( SLA%PQ(K,1:NPQ),GRD%CMDT(SLA%IB(K,1:NPQ),SLA%JB(K,1:NPQ)) )
  ENDIF
  IF(NCONF.EQ.202) THEN
      SLA%INC(K) = &
      & OSUM( SLA%PQ(K,1:NPQ),GRD%CMDT(SLA%IB(K,1:NPQ),SLA%JB(K,1:NPQ)) )
  ENDIF

 ENDIF

ENDDO CYOBS
#ifdef SHARED_MEMORY
!$OMP END DO
!$OMP END PARALLEL
#endif

 CALL MYFRTPROF_WALL('OBS_SLA2: SLA OBSERVATION OPERATOR',1)
END SUBROUTINE OBS_SLA2
