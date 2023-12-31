#undef SHARED_MEMORY
SUBROUTINE THINOUT_SSS(ITER,ZIN1,ZIN2)

!.. THINNING OF SSS OBSERVATIONS
!..
!.. FOR EACH SATELLITE DEFINE 2D THINNING
!.. BOXES, WITHIN WHICH ONLY THE OBS CLOSER TO
!.. ANALYSIS TIME IS RETAINED.
!..
!.. ANDREA STORTO *INGV* 2009-04-6

USE SET_KND
USE OBSDEF
USE OBS_STR, ONLY : SSS
USE GRD_STR
USE IOUNITS, ONLY : IOUNLOG,IOUNOUT
USE MYFRTPROF

IMPLICIT NONE

 INTEGER(I4),INTENT(IN)   ::  ITER
 REAL(R8),   INTENT(IN)   ::  ZIN1, ZIN2
 REAL(R8) :: XRES,YRES,ZAUX,ZM
 INTEGER(I4)   ::  K,JOBS,IIND,JIND,KOKTH,KKOTH,NX,NY,MAXOBS,OIND
 INTEGER(I4)   ::  JX,JY,JSAT,KEVE_HTHN,PKN
 REAL(DP), ALLOCATABLE :: ZDIST(:,:)
 INTEGER(I4), ALLOCATABLE :: KCOUNTTH(:,:),KGOOD(:,:),KOBS(:,:,:)

#include "obs_events.h"

 CALL MYFRTPROF_WALL('THINOUT_SSS: SSS DATA THINNING',0)

  IF(SSS%NO.LE.0) THEN
    CALL MYFRTPROF_WALL('THINOUT_SSS: SSS DATA THINNING',1)
    RETURN
  ENDIF


  WRITE(IOUNOUT,*) ' *** SSS THINNING - ITERATION ',ITER
  WRITE(IOUNLOG,*) ' *** SSS THINNING - ITERATION ',ITER

  KOKTH=0
  KKOTH=0
  XRES=ZIN1
  YRES=ZIN2

  IF(ITER.EQ.1) THEN
     PKN=1
     KEVE_HTHN = KEVE_HTH1
  ELSEIF(ITER.EQ.2) THEN
     PKN=2
     KEVE_HTHN = KEVE_HTH2
  ELSE
     CALL ABOR1('THINOUT: UNSUPPORTED NUMBER OF ITERATION')
  ENDIF

  NX=INT( GRD%IM / XRES ) + 1
  NY=INT( GRD%JM / YRES ) + 1

  MAXOBS = 1000*SSS%NO/(NX*NY)

  ALLOCATE(KCOUNTTH(NX,NY),  &
         & ZDIST(NX,NY),  &
         & KGOOD(NX,NY) )


 SSSSAT : DO JSAT=1,NOSSSSATS

  IF(SSS%ISSSOBS(JSAT) .LE. 0 ) CYCLE SSSSAT

  MAXOBS = COUNT(SSS%FLC(1:SSS%NO).EQ.1 .AND. SSS%KSAT(1:SSS%NO) .EQ. JSAT)
  WRITE(IOUNLOG,*) '    SAT ', JSAT, ' MAXOBS = ',MAXOBS
  ALLOCATE( KOBS(NX,NY,MAXOBS) )

  WRITE(IOUNOUT,*) 'THINNING FOR SATELLITE '//TRIM(CSSSSATID(JSAT))
  WRITE(IOUNLOG,*) 'THINNING FOR SATELLITE '//TRIM(CSSSSATID(JSAT)),&
  & ', FOUND',SSS%ISSSOBS(JSAT),' OBS'

  KCOUNTTH=0
  ZDIST = 1000000._DP
  KGOOD=0

#ifdef SHARED_MEMORY
!$OMP PARALLEL DEFAULT(SHARED), PRIVATE(K,IIND,JIND,....)
!$OMP DO SCHEDULE(STATIC,1)
#endif
  CYOBS5: DO K=1,SSS%NO

       IF(SSS%FLC(K).EQ.0 .OR. SSS%KSAT(K) .NE. JSAT) CYCLE CYOBS5

       !... COORDINATES ON PROJECTED PLAN
       IIND = INT((REAL(SSS%IB(K,PKN),KIND=R8))/XRES) + 1
       JIND = INT((REAL(SSS%JB(K,PKN),KIND=R8))/YRES) + 1

       KCOUNTTH(IIND,JIND) = KCOUNTTH(IIND,JIND) + 1
       OIND = KCOUNTTH(IIND,JIND)
       IF(OIND.GT. MAXOBS) &
       & CALL ABOR1('THINOUT: MAXOBS MUST BE INCREASED')
       KOBS(IIND,JIND,OIND) = K

       IF(ABS(SSS%TDIST(K)) .LT. ZDIST(IIND,JIND) ) THEN
              ZDIST(IIND,JIND) = ABS(SSS%TDIST(K))
              KGOOD(IIND,JIND) = K
       ENDIF

  ENDDO CYOBS5
#ifdef SHARED_MEMORY
!$OMP END DO
!$OMP END PARALLEL
#endif

  DO JY=1,NY
    CYX : DO JX=1,NX

      IF(KCOUNTTH(JX,JY) .GT. 0 ) THEN
          IF (KGOOD(JX,JY).EQ.0 ) &
          & CALL ABOR1('SSS THINNING, PROBLEMS WITH MINIMUM DIST')

            DO K= 1, KCOUNTTH(JX,JY)
               IF( KGOOD(JX,JY) .NE. KOBS(JX,JY,K) ) THEN
                   JOBS = KOBS(JX,JY,K)
                   SSS%FLC(JOBS) = 0
                   SSS%EVE(JOBS) = KEVE_HTHN
                   KKOTH = KKOTH +1
               ELSE
                   KOKTH = KOKTH + 1
               ENDIF
            ENDDO

      ENDIF

    ENDDO CYX
  ENDDO

  DEALLOCATE( KOBS )

 ENDDO SSSSAT

  WRITE(IOUNOUT,*) ' END SSS THINNING - ITERATION ',ITER
  WRITE(IOUNLOG,*)
  WRITE(IOUNLOG,*) ' *** THINNING REPORT - ITERATION ',ITER
  WRITE(IOUNLOG,*) ' THINNING GRID : ',NX,' X',NY
  WRITE(IOUNLOG,*) ' FILTERED-OUT  : ',KKOTH
  WRITE(IOUNLOG,*) ' RETAINED      : ',KOKTH
  WRITE(IOUNLOG,*) ' TOTAL         : ',KOKTH + KKOTH
  WRITE(IOUNLOG,*) ' FLAGS COUNT   : ',SUM(SSS%FLC(1:SSS%NO))

  SSS%NC = KOKTH

  DEALLOCATE(KCOUNTTH,  &
         & ZDIST,  &
         & KGOOD )

  CALL MYFRTPROF_WALL('THINOUT_SSS: SSS DATA THINNING',1)

END SUBROUTINE THINOUT_SSS
