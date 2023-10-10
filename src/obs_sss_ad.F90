#undef SHARED_MEMORY
SUBROUTINE OBS_SSS_AD

!-----------------------------------------------------------------------
!                                                                      !
! APPLY OBSERVATIONAL OPERATOR FOR INS DATA, ADJOINT VERSION           !
!                                                                      !
! VERSION 1: A.STORTO 2009                                             !
! VERSION 2: A.STORTO 2009, OPENMP PARALLELIZATION                     !
!-----------------------------------------------------------------------

 USE SET_KND
 USE GRD_STR
 USE OBSDEF
 USE OBS_STR
 USE CTL_STR
 USE IOUNITS , ONLY : IOUNERR
 USE MYFRTPROF, ONLY : MYFRTPROF_WALL
 USE TLAD_VARS
 USE WEAKLY

 IMPLICIT NONE

 INTEGER(I4)   ::  I, JP, D, KK,I2,OBS_K,NO(SSS%NO)
 INTEGER(I4)   ::  XIND1

CALL MYFRTPROF_WALL('OBS_SSS_AD: ADJOINT OF SSS OPERATOR',0)

!... WATCHOUT : KK FOLLOWS SSS/SLA/ ETC ARRAYS
!               K  FOLLOWS GLOBAL OBSERVATION VECTOR

CYOBS : DO KK = 1,SSS%NO

 IF(SSS%FLC(KK).EQ.0)CYCLE CYOBS

  IF( .NOT. LL_4DVAR .OR. ( LL_4DVAR .AND. SSS%TATS(KK) .EQ. TLAD_TS )) THEN

  OBS%K = OBS%K + 1
  OBS_K = OBS%K
  IF( LL_4DVAR ) OBS_K = SSS%OBIN(KK)

    DO JP=1,NPQ
      GRD%SAL_AD(SSS%IB(KK,JP),SSS%JB(KK,JP),1) = &
      & GRD%SAL_AD(SSS%IB(KK,JP),SSS%JB(KK,JP),1) + SSS%PQ(KK,JP)*OBS%GRA(OBS_K)
    ENDDO

    IF(LL_WEAKLY) THEN
     DO JP=1,NPQ
      SERR_AD(SSS%IB(KK,JP),SSS%JB(KK,JP),1) = &
      & SERR_AD(SSS%IB(KK,JP),SSS%JB(KK,JP),1) + SSS%PQ(KK,JP)*OBS%GRA(OBS_K)
     ENDDO
    ENDIF

  ENDIF

ENDDO CYOBS

! CALL FLUSH(1203)
CALL MYFRTPROF_WALL('OBS_SSS_AD: ADJOINT OF SSS OPERATOR',1)
END SUBROUTINE OBS_SSS_AD
