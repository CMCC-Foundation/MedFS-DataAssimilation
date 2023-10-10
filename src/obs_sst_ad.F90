SUBROUTINE OBS_SST_AD

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

 INTEGER(I4)   ::  I, J, D, KK,I2,OBS_K,NO(SST%NO),JP
 INTEGER(I4) :: XIND1

CALL MYFRTPROF_WALL('OBS_SST_AD: ADJOINT OF SST OPERATOR',0)

!... WATCHOUT : KK FOLLOWS SST/SLA/ ETC ARRAYS
!               K  FOLLOWS GLOBAL OBSERVATION VECTOR

 CYOBS : DO KK = 1,SST%NO

  IF(SST%FLC(KK).EQ.0) CYCLE CYOBS
  IF( .NOT. LL_4DVAR .OR. ( LL_4DVAR .AND. SST%TATS(KK) .EQ. TLAD_TS )) THEN

  OBS%K = OBS%K + 1
  OBS_K = OBS%K
  IF( LL_4DVAR ) OBS_K = SST%OBIN(KK)

    DO JP=1,NPQ
      GRD%TEM_AD(SST%IB(KK,JP),SST%JB(KK,JP),1) = &
      & GRD%TEM_AD(SST%IB(KK,JP),SST%JB(KK,JP),1) + SST%PQ(KK,JP)*OBS%GRA(OBS_K)
    ENDDO
    IF(LL_WEAKLY) THEN
     DO JP=1,NPQ
      TERR_AD(SST%IB(KK,JP),SST%JB(KK,JP),1) = &
      & TERR_AD(SST%IB(KK,JP),SST%JB(KK,JP),1) + SST%PQ(KK,JP)*OBS%GRA(OBS_K)
     ENDDO
    ENDIF

  ENDIF

 ENDDO CYOBS

CALL MYFRTPROF_WALL('OBS_SST_AD: ADJOINT OF SST OPERATOR',1)
END SUBROUTINE OBS_SST_AD
