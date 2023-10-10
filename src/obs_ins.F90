#undef SHARED_MEMORY
SUBROUTINE OBS_INS

!-----------------------------------------------------------------------
!                                                                      !
! APPLY OBSERVATIONAL OPERATOR FOR INS DATA                            !
!                                                                      !
! VERSION 1: A.STORTO 2009                                             !
!-----------------------------------------------------------------------


 USE SET_KND
 USE GRD_STR
 USE OBSDEF
 USE OBS_STR
 USE CTL_STR
 USE READFG
 USE IOUNITS , ONLY : IOUNERR,IOUNLOG
 USE RUN, ONLY : RTIMES1950,NTSTEPS,LL_DENSLEVS
 USE MYFRTPROF, ONLY : MYFRTPROF_WALL

 IMPLICIT NONE

 INTEGER(I4)   ::  I,J,D,KK,JTSTEP,KSTEP,JP,FIRSTTS,JT
 REAL(R8)      ::  SALB, ZSUM

#include "obs_events.h"

CALL MYFRTPROF_WALL('OBS_INS: COMPUTE INSITU BACKGROUND VALUES',0)

FIRST : DO JTSTEP=1,NTSTEPS
  IF( LL_FGASSIM(JTSTEP) ) THEN
      FIRSTTS=JTSTEP
      EXIT FIRST
  ENDIF
ENDDO FIRST

WRITE(IOUNLOG,*) ' SURFACE OBS :',COUNT( INS%PAR(1:INS%NO) .EQ. KKSST ),&
& COUNT( INS%DSRC(1:INS%NO) .EQ. 'DRF' ),&
& COUNT( INS%DSRC(1:INS%NO) .EQ. 'DRF' .AND. INS%FLC(1:INS%NO) .EQ. 1)

WHERE( INS%PAR(1:INS%NO) .EQ. KKSST )
       INS%ERR(1:INS%NO) = INS%DPT(1:INS%NO)
       INS%DPT(1:INS%NO) = 0._R8
       INS%PAR(1:INS%NO) = KKTEMP
ENDWHERE

#ifdef SHARED_MEMORY
!$OMP PARALLEL DEFAULT(SHARED),PRIVATE(KK,KSTEP,JTSTEP,D,JP)
!$OMP DO SCHEDULE(DYNAMIC,1)
#endif
 CYOBS : DO KK = 1,INS%NO

  IF(INS%FLC(KK).NE.1 ) CYCLE CYOBS

  KSTEP = MINLOC( ABS(RTIMES1950(1:NTSTEPS)-INS%TIM(KK) ), DIM=1 ) - (FIRSTTS-1)
  IF(KSTEP .LT. 1 .OR. KSTEP .GT. NFGTSTEPS) THEN
     WRITE(IOUNOUT,*) ' OBS_INS: FOUND KSTEP EQUAL TO ',KSTEP
     WRITE(IOUNOUT,*) ' OBS_INS: TIME EQUAL TO ',INS%TIM(KK)
     WRITE(IOUNOUT,*) ' LL_FGASSIM : ',LL_FGASSIM
     WRITE(IOUNOUT,*) ' RTIMES1950 : ',RTIMES1950(1:NTSTEPS)
     CALL ABOR1('OBS_INS: KSTEP OUT OF RANGE')
  ENDIF

  D=INS%KB(KK)

  IF(D .LT. 1 .OR. D .GT. GRD%KM-1) THEN
     WRITE(IOUNOUT,*) ' OBS_INS: FOUND D EQUAL TO ',D
     CALL ABOR1('OBS_INS: D OUT OF RANGE')
  ENDIF

  IF( INS%PAR(KK).EQ.KKTEMP)THEN

    IF(LL_DENSLEVS) THEN
      ZSUM = 0._R8
      DO JP=1,NPQ
        IF( GRD%TEMB(INS%IB(KK,JP),INS%JB(KK,JP),D,KSTEP) .NE. 0._R8 ) THEN
            ZSUM = ZSUM + INS%PQ(KK,JP)
        ELSE
            INS%PQ(KK,JP) = 0._R8
        ENDIF
        IF( GRD%TEMB(INS%IB(KK,JP),INS%JB(KK,JP),D+1,KSTEP) .NE. 0._R8 ) THEN
            ZSUM = ZSUM + INS%PQ(KK,JP+NPQ)
        ELSE
            INS%PQ(KK,JP+NPQ) = 0._R8
        ENDIF
      ENDDO
      IF( ZSUM .GT. 0._R8 ) THEN
            INS%PQ(KK,:) = INS%PQ(KK,:) / ZSUM
      ELSE
            INS%FLC(KK) = 0
            INS%EVE(KK) = KEVE_UDEN
      ENDIF
    ENDIF

    INS%BAC(KK) = 0._R8
    DO JP=1,NPQ
       INS%BAC(KK) = INS%BAC(KK) + INS%PQ(KK,JP)*GRD%TEMB(INS%IB(KK,JP),INS%JB(KK,JP),D,KSTEP) + &
                     INS%PQ(KK,JP+NPQ)*GRD%TEMB(INS%IB(KK,JP),INS%JB(KK,JP),D+1,KSTEP)
    ENDDO
    IF( INS%VAL(KK) .GT. 150._R8 ) THEN
       SALB = 0._R8
       DO JP=1,NPQ
          SALB = SALB + INS%PQ(KK,JP)*GRD%SALB(INS%IB(KK,JP),INS%JB(KK,JP),D,KSTEP) + &
          INS%PQ(KK,JP+NPQ)*GRD%SALB(INS%IB(KK,JP),INS%JB(KK,JP),D+1,KSTEP)
       ENDDO
       INS%VAL(KK) = POTEMP(SALB,INS%VAL(KK)-200._R8,DEP_TO_P(INS%DPT(KK),&
       & INS%LAT(KK)),0._R8)
    ENDIF

    INS%RES(KK) = INS%VAL(KK) - INS%BAC(KK)

  ELSE IF(INS%PAR(KK).EQ.KKSAL )THEN

    INS%BAC(KK) = 0._R8

    IF(LL_DENSLEVS) THEN
      ZSUM = 0._R8
      DO JP=1,NPQ
        IF( GRD%SALB(INS%IB(KK,JP),INS%JB(KK,JP),D,KSTEP) .NE. 0._R8 ) THEN
            ZSUM = ZSUM + INS%PQ(KK,JP)
        ELSE
            INS%PQ(KK,JP) = 0._R8
        ENDIF
        IF( GRD%SALB(INS%IB(KK,JP),INS%JB(KK,JP),D+1,KSTEP) .NE. 0._R8 ) THEN
            ZSUM = ZSUM + INS%PQ(KK,JP+NPQ)
        ELSE
            INS%PQ(KK,JP+NPQ) = 0._R8
        ENDIF
      ENDDO
      IF( ZSUM .GT. 0._R8 ) THEN
            INS%PQ(KK,:) = INS%PQ(KK,:) / ZSUM
      ELSE
            INS%FLC(KK) = 0
            INS%EVE(KK) = KEVE_UDEN
      ENDIF
    ENDIF

    DO JP=1,NPQ
       INS%BAC(KK) = INS%BAC(KK) + INS%PQ(KK,JP)*GRD%SALB(INS%IB(KK,JP),INS%JB(KK,JP),D,KSTEP) + &
                     INS%PQ(KK,JP+NPQ)*GRD%SALB(INS%IB(KK,JP),INS%JB(KK,JP),D+1,KSTEP)
    ENDDO

    INS%RES(KK) = INS%VAL(KK) - INS%BAC(KK)

  ELSEIF( INS%PAR(KK).EQ.KKT2M ) THEN

    INS%BAC(KK) = 0._R8

    IF(INS%OTYPE(KK).EQ.KKMOORIN) THEN
     DO JT=1,NFGTSTEPS
      DO JP=1,NPQ
       INS%BAC(KK) = INS%BAC(KK) + INS%PQ(KK,JP)*GRD%T2MB(INS%IB(KK,JP),INS%JB(KK,JP),JT)/REAL(NFGTSTEPS,R8)
      ENDDO
     ENDDO
    ELSE
     DO JP=1,NPQ
       INS%BAC(KK) = INS%BAC(KK) + INS%PQ(KK,JP)*GRD%T2MB(INS%IB(KK,JP),INS%JB(KK,JP),KSTEP) 
     ENDDO
    ENDIF

    INS%RES(KK) = INS%VAL(KK) - INS%BAC(KK)

  ELSEIF( INS%PAR(KK).EQ.KKQ2M ) THEN

    INS%BAC(KK) = 0._R8

    DO JP=1,NPQ
       INS%BAC(KK) = INS%BAC(KK) + INS%PQ(KK,JP)*GRD%Q2MB(INS%IB(KK,JP),INS%JB(KK,JP),KSTEP) 
    ENDDO

    INS%RES(KK) = INS%VAL(KK) - INS%BAC(KK)

  ELSE

     WRITE(IOUNERR,*) 'PARAMETER ',INS%PAR(KK),' NOT SUPPORTED'
     CALL ABOR1('OBS_INS : UNKNOWN PARAMETER')

  ENDIF

 ENDDO CYOBS
#ifdef SHARED_MEMORY
!$OMP END DO
!$OMP END PARALLEL
#endif

 WRITE(IOUNLOG,*)
 WRITE(IOUNLOG,*) ' *** INS OBSERVATIONS MISIFTS COMPUTED'
 CALL FLUSH(IOUNLOG)

CALL MYFRTPROF_WALL('OBS_INS: COMPUTE INSITU BACKGROUND VALUES',1)
 CONTAINS
   REAL(KIND=R8) FUNCTION potemp( ps, pt, pp, ppr )

      REAL(KIND=R8), INTENT(IN) :: ps
      REAL(KIND=R8), INTENT(IN) :: pt
      REAL(KIND=R8), INTENT(IN) :: pp
      REAL(KIND=R8), INTENT(IN) :: ppr

      REAL(KIND=R8) :: zpol
      REAL(KIND=R8), PARAMETER :: a1 =  1.067610e-05
      REAL(KIND=R8), PARAMETER :: a2 = -1.434297e-06
      REAL(KIND=R8), PARAMETER :: a3 = -7.566349e-09
      REAL(KIND=R8), PARAMETER :: a4 = -8.535585e-06
      REAL(KIND=R8), PARAMETER :: a5 =  3.074672e-08
      REAL(KIND=R8), PARAMETER :: a6 =  1.918639e-08
      REAL(KIND=R8), PARAMETER :: a7 =  1.788718e-10

      zpol = a1 + a2 * ps + a3 * ( pp + ppr ) + a4 * pt &
         & + a5 * ps * pt + a6 * pt * pt + a7 * pt * ( pp + ppr )

      potemp = pt + ( pp - ppr ) * zpol

   END FUNCTION potemp

   REAL(KIND=R8) FUNCTION dep_to_p( p_dep, p_lat )
      REAL(KIND=R8), INTENT(IN) :: p_dep    ! Depth in meters
      REAL(KIND=R8), INTENT(IN) :: p_lat    ! Latitude in degrees
      REAL(KIND=R8) :: z_x
      REAL(KIND=R8) :: z_c1
      REAL(KIND=R8) :: z_c2
      REAL(KIND=R8) :: z_d

      z_x = SIN( p_lat / 57.29578 )
      z_x = z_x * z_x
      z_c1 = ( 5.92  + 5.25 * z_x ) * 1e-3
      z_c2 = 2.21e-6
      z_d = ( z_c1 - 1 ) * ( z_c1 - 1  ) - 4 * z_c2 * p_dep
      dep_to_p = (( 1 - z_c1 ) - SQRT( z_d )) / ( 2 * z_c2 )
   END FUNCTION dep_to_p

END SUBROUTINE OBS_INS
