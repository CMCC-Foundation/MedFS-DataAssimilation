SUBROUTINE OBS_SLA2E

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

 INTEGER(I4)   ::  K

 CALL MYFRTPROF_WALL('OBS_SLA2E: SLA OBSERVATION OPERATOR',0)

#ifdef SHARED_MEMORY
!$OMP PARALLEL DEFAULT(SHARED), PRIVATE(K)
!$OMP DO SCHEDULE(DYNAMIC,1)
#endif
 CYOBS : DO K=1,SLA%NO

  IF( .NOT. LL_4DVAR .OR. ( LL_4DVAR .AND. SLA%TATS(K) .EQ. TLAD_TS )) THEN

  IF(NCONF.NE.202) THEN

    SLA%INC(K) = &
    & OSUM( SLA%PQ(K,1:NPQ),GRD%ETA(SLA%IB(K,1:NPQ),SLA%JB(K,1:NPQ)) )

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

 CALL MYFRTPROF_WALL('OBS_SLA2E: SLA OBSERVATION OPERATOR',1)
END SUBROUTINE OBS_SLA2E
