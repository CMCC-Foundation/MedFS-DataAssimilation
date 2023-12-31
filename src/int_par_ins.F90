SUBROUTINE INT_PAR_INS

!-----------------------------------------------------------------------
!                                                                      !
! GET INTERPOLATION PARAMETERS FOR A GRID                              !
!                                                                      !
! VERSION 1: S.DOBRICIC 2006                                           !
! VERSION 2: A.S. ADAPTED FOR INS OBS IN A GENERIC GRID                     !
!-----------------------------------------------------------------------

 USE SET_KND
 USE GRD_STR
 USE PHINTERP2
#ifndef NECSX
! USE PHINTERP
#endif
 USE OBSDEF
 USE OMP_LIB
 USE OBS_STR
 USE IOUNITS, ONLY : IOUNOUT, IOUNLOG, IOUNERR
 USE RUN
 USE OBSHANDLING, ONLY : ZMAXDEP
 USE MYFRTPROF, ONLY : MYFRTPROF_WALL

 IMPLICIT NONE

  INTEGER(I4)   ::  I, K,II, JJ,IM,JM,MSCOUNT,JLEV,KLEV, JP
  INTEGER(I4)   ::  I1, J1, K1, IDEP,I2,IAUX
  INTEGER(I4)   ::  XIND1,KA,KAG
  REAL(R8)      ::  P1, Q1, R1, ZZSS, PQT(4)
  REAL(R8)      ::  MSK4, DIV_X, DIV_Y, NEWLON
  INTEGER(I4)   :: MYT,NNT,KSTART,KEND,NCHNK(999),NSIZE,JTHRD

#include "obs_events.h"

CALL MYFRTPROF_WALL('INT_PAR_INS: SET INTERPOLATION FOR INSITU',0)

  IF(INS%NO.LE.0) THEN
        WRITE(IOUNLOG,*) '*** WARNING: ROUTINE INT_PAR_INS: NO OBS!'
        CALL MYFRTPROF_WALL('INT_PAR_INS: SET INTERPOLATION FOR INSITU',1)
        RETURN
  ENDIF

  CALL PREPINTERP_INIT

  MSCOUNT = 0
  INS%IB = 0
  INS%JB = 0
  IM=GRD%IM
  JM=GRD%JM
  INS%FLC=0
  INS%TDIST(1:INS%NO)=INS%TIM(1:INS%NO) - ZANJUL1950

#ifdef SHARED_MEMORY
!$OMP PARALLEL DEFAULT(SHARED), PRIVATE(K,KSTART,KEND,NNT,MYT,NSIZE,&
!$OMP & NCHNK,JTHRD )
  NNT = OMP_GET_NUM_THREADS()
  MYT = OMP_GET_THREAD_NUM() + 1
  NSIZE = INS%NO / NNT
  NCHNK(:) = NSIZE
  DO JTHRD = 1, MOD( INS%NO, NNT )
    NCHNK(JTHRD) = NCHNK(JTHRD) + 1
  ENDDO
  IF( MYT .NE. 1 ) THEN
     KSTART = SUM(NCHNK(1:MYT-1)) + 1
     KEND   = SUM(NCHNK(1:MYT))
  ELSE
     KSTART = 1
     KEND = NCHNK(1 )
  ENDIF
  WRITE(*,*) ' THREAD ', MYT, ' LIMITS ', KSTART,KEND
#else
  KSTART=1
  KEND=INS%NO
#endif
  OBS_LOOP1 : DO K = KSTART,KEND

  ! CHECK IF ALREADY COMPUTED
      IF(K.GT.KSTART) THEN
            IF(INS%PROF(K).EQ. INS%PROF(K-1) .AND. INS%FLC(K-1) .EQ. 1 ) THEN
               INS%IB(K,:) = INS%IB(K-1,:)
               INS%JB(K,:) = INS%JB(K-1,:)
               INS%PQ(K,1:NPQ) = INS%PQ(K-1,1:NPQ)
               CYCLE OBS_LOOP1
            ENDIF
      ENDIF

      IF( INS%EVE(K) .EQ. KEVE_PDOM ) CYCLE OBS_LOOP1

#ifndef NECSX
!      IF( NN_INTERPM .EQ. 1 ) THEN
!         CALL PREPINTERP(INS%LON(K),INS%LAT(K),INS%IB(K,:),INS%JB(K,:),INS%PQ(K,1:NPQ))
!      ELSEIF ( NN_INTERPM .EQ. 2 ) THEN
         CALL PREPINTERP2(INS%LON(K),INS%LAT(K),INS%IB(K,:),INS%JB(K,:),INS%PQ(K,1:NPQ))
!      ENDIF
#else
      CALL PREPINTERP2(INS%LON(K),INS%LAT(K),INS%IB(K,:),INS%JB(K,:),INS%PQ(K,1:NPQ))
#endif

#ifdef INTERP_DEBUG
      IF( ANY( INS%IB(K,:) .EQ. 0 ) ) THEN
          WRITE(2001,*) K, INS%LON(K),INS%LAT(K)
      ENDIF
#endif

ENDDO OBS_LOOP1
#ifdef SHARED_MEMORY
!$OMP END PARALLEL
!$OMP PARALLEL DEFAULT(SHARED), PRIVATE(K,KLEV,JLEV,ZZSS,R1,JP)
!$OMP DO SCHEDULE(DYNAMIC)
#endif
OBS_LOOP2 : DO K = 1,INS%NO

      IF( ALL(INS%IB(K,:).GE.1 ) .AND. ALL(INS%JB(K,:).GE.1 ) .AND. &
         ALL(INS%IB(K,:).LE.GRD%IM ) .AND. ALL(INS%JB(K,:).LE.GRD%JM ) ) THEN
         INS%FLC(K) = 1
      ELSE
         INS%FLC(K) = 0
         INS%EVE(K) = KEVE_INTE
#ifdef INTERP_DEBUG
         WRITE(9998,*) INS%LON(K),INS%LAT(K),INS%IB(K,:), INS%JB(K,:)
#endif
         CYCLE OBS_LOOP2
      ENDIF

      IF(INS%FLC(K) .EQ. 0) CYCLE OBS_LOOP2
      SELECT CASE(INS%OTYPE(K))
            CASE(100:112,201,401)   ! XBT/MBT
               INS%KTY(K)=KSXBT
            CASE(741)   ! TESAC
               INS%KTY(K)=KSTESAC
            CASE(831)   ! ARGO
               INS%KTY(K)=KSARGO
            CASE(820)   ! BUOYS
               INS%KTY(K)=KSBUOY
            CASE(-109,-114:-111,-124:-119)   ! BUOYS
               INS%KTY(K)=KSAIR2
            CASE DEFAULT
               CALL ABOR1('INTPAR3D : FOUND UNSUPPORTED OBSTYPE')
       ENDSELECT

       IF( INS%DPT(K) .GT. ZMAXDEP(INS%KTY(K)) &
       & .AND. INS%DSRC(K)(1:3) .NE. 'ICO' &
       & .AND. INS%DSRC(K)(1:3) .NE. 'DRF' &
       & .AND. INS%DSRC(K)(1:3) .NE. 'BUF' ) THEN
            INS%FLC(K) = 0
            INS%EVE(K) = KEVE_TOOD
            CYCLE OBS_LOOP2
       ENDIF

  ! FIND VERTICAL LEVELS
       IF( INS%DSRC(K)(1:3) .EQ. 'ICO' .OR. INS%DSRC(K)(1:3) .EQ. 'DRF') THEN
        IF( LL_KEEP_SURFACE ) THEN
         INS%KB(K) = 1
         INS%RB(K) = 0._R8
        ELSE
         INS%FLC(K) = 0
         INS%EVE(K) = KEVE_DEPL
         CYCLE OBS_LOOP2
        ENDIF
       ELSEIF ( INS%DSRC(K)(1:3) .EQ. 'BUF' ) THEN
         INS%KB(K) = 1
         INS%RB(K) = 0._R8
       ELSE 
        KLEV=0
        LEVS : DO JLEV = 1,GRD%KM-1
         IF(INS%DPT(K).GE.GRD%DEP(JLEV).AND.&
         & INS%DPT(K).LT.GRD%DEP(JLEV+1)) THEN
            KLEV = JLEV
            EXIT LEVS
         ENDIF
        ENDDO LEVS

        IF(KLEV.EQ. 0 ) THEN
          IF(INS%DPT(K).LT.0._R8 .OR. .NOT. LL_KEEP_SURFACE) THEN
            INS%FLC(K) = 0
            INS%EVE(K) = KEVE_DEPL
            CYCLE OBS_LOOP2
          ELSE
            INS%KB(K) = 1
            INS%RB(K) = 0._R8
          ENDIF
        ELSE
          INS%KB(K) = KLEV
          INS%RB(K) = MAX(1.E-12_R8,ABS((INS%DPT(K) - GRD%DEP(KLEV)))) / &
                 & MAX(1.E-12_R8,ABS((GRD%DEP(KLEV+1) - GRD%DEP(KLEV))))
        ENDIF
       ENDIF

       R1=INS%RB(K)

       PQT(1:4) = INS%PQ(K,1:4)

       INS%PQ(K,1) = (1._R8-R1) * PQT(1) * GRD%MSK(INS%IB(K,1),INS%JB(K,1),INS%KB(K))
       INS%PQ(K,2) = (1._R8-R1) * PQT(2) * GRD%MSK(INS%IB(K,2),INS%JB(K,2),INS%KB(K))
       INS%PQ(K,3) = (1._R8-R1) * PQT(3) * GRD%MSK(INS%IB(K,3),INS%JB(K,3),INS%KB(K))
       INS%PQ(K,4) = (1._R8-R1) * PQT(4) * GRD%MSK(INS%IB(K,4),INS%JB(K,4),INS%KB(K))
       INS%PQ(K,5) =        R1  * PQT(1) * GRD%MSK(INS%IB(K,1),INS%JB(K,1),INS%KB(K)+1)
       INS%PQ(K,6) =        R1  * PQT(2) * GRD%MSK(INS%IB(K,2),INS%JB(K,2),INS%KB(K)+1)
       INS%PQ(K,7) =        R1  * PQT(3) * GRD%MSK(INS%IB(K,3),INS%JB(K,3),INS%KB(K)+1)
       INS%PQ(K,8) =        R1  * PQT(4) * GRD%MSK(INS%IB(K,4),INS%JB(K,4),INS%KB(K)+1)

       ZZSS=0._R8
       DO JP=1,NPQ
         ZZSS = ZZSS + GRD%MSK(INS%IB(K,JP),INS%JB(K,JP),INS%KB(K))
         ZZSS = ZZSS + GRD%MSK(INS%IB(K,JP),INS%JB(K,JP),INS%KB(K)+1)
       ENDDO

       IF( ZZSS .LT. 1._R8 ) THEN
           INS%FLC(K) = 0
           INS%EVE(K) = KEVE_MASK
       ELSE
           INS%PQ(K,:) = INS%PQ(K,:) / SUM(INS%PQ(K,:))
       ENDIF

ENDDO OBS_LOOP2
#ifdef SHARED_MEMORY
!$OMP END DO
!$OMP END PARALLEL
#endif

MSCOUNT=COUNT(INS%EVE(1:INS%NO).EQ.KEVE_DEPL)
WRITE(IOUNLOG,*) &
& ' INS OBS FILTERED OUT FOR UNVALID DEPTH: ',MSCOUNT

MSCOUNT=COUNT(INS%EVE(1:INS%NO).EQ.KEVE_INTE)
WRITE(IOUNLOG,*) &
& ' INS OBS FILTERED OUT FOR INTERP PROBLEMS: ',MSCOUNT

MSCOUNT=COUNT(INS%EVE(1:INS%NO).EQ.KEVE_TOOD)
WRITE(IOUNLOG,*) &
& ' INS OBS FILTERED OUT FOR MAX DEPTH EXCEEDED: ',MSCOUNT

MSCOUNT=COUNT(INS%EVE(1:INS%NO).EQ.KEVE_MASK)
WRITE(IOUNLOG,*) &
& ' INS OBS FILTERED OUT FOR MASKING INCONSISTENCIES: ',MSCOUNT

! ---
! COUNT GOOD OBSERVATIONS
  INS%NC = COUNT(INS%FLC(1:INS%NO) .EQ. 1)

  IF(OBS%AIR .EQ. 1) THEN
    KA=0
    KAG=0
    DO K = 1,INS%NO
       IF(INS%OTYPE(K).LT.0) THEN
         KA=KA+1
          IF(INS%FLC(K).EQ.1) KAG=KAG+1
       ENDIF
    ENDDO
    WRITE(IOUNLOG,*)
    WRITE(IOUNLOG,*) ' AIR OBS REPORT -- PART 1'
    WRITE(IOUNLOG,*) ' TOTAL : ', KA
    WRITE(IOUNLOG,*) ' GOOD  : ', KAG
    WRITE(IOUNLOG,*) 

    KA=0
    KAG=0
    DO K = 1,INS%NO
       IF(INS%KTY(K).EQ.KSAIR2) THEN
         KA=KA+1
          IF(INS%FLC(K).EQ.1) KAG=KAG+1
       ENDIF
    ENDDO
    WRITE(IOUNLOG,*)
    WRITE(IOUNLOG,*) ' AIR OBS REPORT -- PART 1'
    WRITE(IOUNLOG,*) ' TOTAL : ', KA
    WRITE(IOUNLOG,*) ' GOOD  : ', KAG
    WRITE(IOUNLOG,*)
  ENDIF

  WRITE(IOUNLOG,*)
  WRITE(IOUNLOG,*) ' *** INS OBS AFTER INTERPOLATION SET-UP'
  WRITE(IOUNLOG,*) ' TOTAL NUMBER OF INS OBS           :',INS%NO
  WRITE(IOUNLOG,*) ' TOTAL NUMBER OF RETAINED INS  OBS :',INS%NC

CALL MYFRTPROF_WALL('INT_PAR_INS: SET INTERPOLATION FOR INSITU',1)
END SUBROUTINE INT_PAR_INS
