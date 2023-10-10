SUBROUTINE OBS_SST_TL

!-----------------------------------------------------------------------
!                                                                      !
! APPLY OBSERVATIONAL OPERATOR FOR SST DATA, TANGENT-LINEAR VERSION    !
!                                                                      !
! VERSION 1: A.STORTO 2009                                             !
!-----------------------------------------------------------------------


 USE SET_KND
 USE GRD_STR
 USE OBS_STR
 USE CTL_STR
 USE IOUNITS , ONLY : IOUNERR
 USE TLAD_VARS
 USE MYFRTPROF

 IMPLICIT NONE

 INTEGER(I4)   ::  I, J, KK,I2
 INTEGER(I4) :: XIND1

CALL MYFRTPROF_WALL('OBS_SST_TL: SST OBSERVATION OPERATOR',0)

#ifdef SHARED_MEMORY
!$OMP PARALLEL DEFAULT(SHARED), PRIVATE(KK)
!$OMP DO SCHEDULE(STATIC,1)
#endif
CYOBS : DO KK = 1,SST%NO
 
  IF( .NOT. LL_4DVAR .OR. ( LL_4DVAR .AND. SST%TATS(KK) .EQ. TLAD_TS )) THEN

    SST%INC(KK) = &
    OSUM( SST%PQ(KK,1:NPQ),GRD%TEM(SST%IB(KK,1:NPQ),SST%JB(KK,1:NPQ),1) )

  ENDIF

 ENDDO CYOBS
#ifdef SHARED_MEMORY
!$OMP END DO
!$OMP END PARALLEL
#endif

CALL MYFRTPROF_WALL('OBS_SST_TL: SST OBSERVATION OPERATOR',1)

END SUBROUTINE OBS_SST_TL
