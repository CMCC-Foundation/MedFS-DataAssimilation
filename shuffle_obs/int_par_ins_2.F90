#undef SHARED_MEMORY
SUBROUTINE INT_PAR_INS_2

!-----------------------------------------------------------------------
!                                                                      !
! GET INTERPOLATION PARAMETERS FOR A GRID                              !
!                                                                      !
! VERSION 1: S.DOBRICIC 2006                                           !
! VERSION 2: A.S. ADAPTED FOR INS OBS IN A GENERIC GRID                     !
!-----------------------------------------------------------------------

 USE SET_KND, ONLY : R8, I4
 USE PHINTERP
 USE OBSDEF
 USE OMP_LIB
 USE OBS_STR

 IMPLICIT NONE

  INTEGER(I4)   ::  NX, NY, NZ
  INTEGER(I4)   ::  I, K,II, JJ,MSCOUNT,JLEV,KLEV, JP
  INTEGER(I4)   ::  I1, J1, K1, IDEP,I2,IAUX
  INTEGER(I4)   ::  XIND1
  REAL(R8)      ::  P1, Q1, R1, ZZSS
  REAL(R8)      ::  MSK4, DIV_X, DIV_Y, NEWLON
  REAL(R8), ALLOCATABLE   ::  GRD_MSK(:,:,:)
  INTEGER(I4)   :: MYT,NNT,KSTART,KEND,NCHNK(999),NSIZE,JTHRD

  IF(INS%NO.LE.0) THEN
        WRITE(6,*) '*** WARNING: ROUTINE INT_PAR_INS: NO OBS!'
        RETURN
  ENDIF

  CALL GETNCDIM('GRID_CR.nc','x',NX)
  CALL GETNCDIM('GRID_CR.nc','y',NY)
  CALL GETNCDIM('GRID_CR.nc','z',NZ)
  ALLOCATE( GRD_MSK(NX,NY,NZ) )
  CALL GETNCVAR('GRID_CR.nc','tmsk',NX,NY,NZ,GRD_MSK)

  INS%IB = 0
  INS%JB = 0

#ifdef SHARED_MEMORY
!$OMP PARALLEL DEFAULT(PRIVATE), SHARED(INS)
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

      CALL PREPINTERP(INS%LON(K),INS%LAT(K),INS%IB(K,:),INS%JB(K,:),INS%PQ(K,1:NPQ))

ENDDO OBS_LOOP1
#ifdef SHARED_MEMORY
!$OMP END PARALLEL
!$OMP PARALLEL DEFAULT(SHARED), PRIVATE(K,KLEV,JLEV,ZZSS,R1)
!$OMP DO SCHEDULE(DYNAMIC)
#endif
OBS_LOOP2 : DO K = 1,INS%NO

      INS%FLC(K) = 1

  ! FIND VERTICAL LEVELS

       R1=INS%RB(K)

       INS%PQ(K,1) = (1._R8-R1) * INS%PQ(K,1) * GRD_MSK(INS%IB(K,1),INS%JB(K,1),INS%KB(K))
       INS%PQ(K,2) = (1._R8-R1) * INS%PQ(K,2) * GRD_MSK(INS%IB(K,2),INS%JB(K,2),INS%KB(K))
       INS%PQ(K,3) = (1._R8-R1) * INS%PQ(K,3) * GRD_MSK(INS%IB(K,3),INS%JB(K,3),INS%KB(K))
       INS%PQ(K,4) = (1._R8-R1) * INS%PQ(K,4) * GRD_MSK(INS%IB(K,4),INS%JB(K,4),INS%KB(K))
       INS%PQ(K,5) =     R1  * INS%PQ(K,1) * GRD_MSK(INS%IB(K,1),INS%JB(K,1),INS%KB(K)+1)
       INS%PQ(K,6) =     R1  * INS%PQ(K,2) * GRD_MSK(INS%IB(K,2),INS%JB(K,2),INS%KB(K)+1)
       INS%PQ(K,7) =     R1  * INS%PQ(K,3) * GRD_MSK(INS%IB(K,3),INS%JB(K,3),INS%KB(K)+1)
       INS%PQ(K,8) =     R1  * INS%PQ(K,4) * GRD_MSK(INS%IB(K,4),INS%JB(K,4),INS%KB(K)+1)

       INS%PQ(K,:) = INS%PQ(K,:) / SUM(INS%PQ(K,:))

ENDDO OBS_LOOP2
#ifdef SHARED_MEMORY
!$OMP END DO
!$OMP END PARALLEL
#endif

  WRITE(6,*)
  WRITE(6,*) ' *** INS OBS AFTER INTERPOLATION SET-UP'
  WRITE(6,*) ' TOTAL NUMBER OF INS OBS           :',INS%NO

END SUBROUTINE INT_PAR_INS_2
