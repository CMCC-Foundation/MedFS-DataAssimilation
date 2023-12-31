SUBROUTINE OBSOP

!-----------------------------------------------------------------------
!                                                                      !
! APPLY OBSERVATIONAL OPERATORS
!                                                                      !
! VERSION 1: S.DOBRICIC 2006                                           !
!-----------------------------------------------------------------------


 USE SET_KND
 USE OBS_STR
 USE RUN, ONLY : COBSOPT
 USE GRD_STR, ONLY : LL_SSH
 USE MYFRTPROF, ONLY : MYFRTPROF_WALL

 IMPLICIT NONE

 CALL MYFRTPROF_WALL('OBSOP: OBSERVATION OPERATOR',0)

IF(COBSOPT.EQ.'MFS') THEN

! SATELLITE OBSERVATIONS OF SLA
  CALL OBS_SLA

! OBSERVATIONS BY ARGO FLOATS
  CALL OBS_ARG

! OBSERVATIONS BY XBT PROFILES
  CALL OBS_XBT

! OBSERVATIONS BY GLIDERS
  CALL OBS_GLD

! OBSERVATIONS OF VELOCITY
  CALL OBS_VEL

! OBSERVATIONS OF VELOCITY
  IF(TRJ%NO.GT.0) CALL OBS_TRJ

ELSEIF(COBSOPT.EQ.'GLO' .OR. COBSOPT.EQ.'ASC' ) THEN

! SATELLITE OBSERVATIONS OF SLA
  IF ( LL_SSH ) THEN
    CALL OBS_SLA2E
  ELSE
    CALL OBS_SLA2
  ENDIF

! INS OBSERVATIONS
  CALL OBS_INS_TL

! SST OBSERVATIONS
  IF ( SST%NO .GT. 0 )CALL OBS_SST_TL

! SSS OBSERVATIONS
  IF ( SSS%NO .GT. 0 )CALL OBS_SSS_TL

ELSE

  CALL ABOR1('COBSOPT OPTION IN OBSOP '//COBSOPT//' UNSUPPORTED')

ENDIF

 CALL MYFRTPROF_WALL('OBSOP: OBSERVATION OPERATOR',1)
END SUBROUTINE OBSOP
