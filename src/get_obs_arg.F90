SUBROUTINE GET_OBS_ARG

!-----------------------------------------------------------------------
!                                                                      !
! LOAD ARGO OBSERVATIONS                                                !
!                                                                      !
! VERSION 1: S.DOBRICIC 2006                                           !
!-----------------------------------------------------------------------

 USE SET_KND
 USE GRD_STR
 USE OBS_STR
 USE IOUNITS, ONLY : IOUNOUT, IOUNLOG, IOUNERR

 IMPLICIT NONE

  INTEGER(I4)   ::  K
  INTEGER(I4)   ::  I1, KK, I
!  LOGICAL       ::  INS

 WRITE(IOUNLOG,*)
 WRITE(IOUNLOG,*)
 WRITE(IOUNLOG,*) ' ******** GET ARGO OBSERVATIONS'

!  INS(I,I1) = I.GE.1 .AND. I.LT.I1

    ARG%NO = 0
    ARG%NC = 0

  IF(OBS%ARG.EQ.0)RETURN

  OPEN(511,FILE='ARG_MIS.DAT',FORM='UNFORMATTED',STATUS='OLD')

! ---
! ALLOCATE MEMORY FOR OBSERVATIONS

   READ(511) ARG%NO

   WRITE(IOUNOUT,*) ' TOTAL NO OF ARGO OBS: ',ARG%NO

   IF(ARG%NO.EQ.0)RETURN

   ALLOCATE ( ARG%INO(ARG%NO), ARG%FLG(ARG%NO), ARG%FLC(ARG%NO), ARG%PAR(ARG%NO))
   ALLOCATE ( ARG%LON(ARG%NO), ARG%LAT(ARG%NO), ARG%DPT(ARG%NO), ARG%TIM(ARG%NO))
   ALLOCATE ( ARG%VAL(ARG%NO), ARG%BAC(ARG%NO), ARG%INC(ARG%NO))
   ALLOCATE ( ARG%BIA(ARG%NO), ARG%ERR(ARG%NO))
   ALLOCATE ( ARG%RES(ARG%NO), ARG%B_A(ARG%NO))
   ALLOCATE ( ARG%IB(ARG%NO), ARG%JB(ARG%NO), ARG%KB(ARG%NO))
   ALLOCATE ( ARG%PB(ARG%NO), ARG%QB(ARG%NO), ARG%RB(ARG%NO))
   ALLOCATE ( ARG%PQ1(ARG%NO), ARG%PQ2(ARG%NO), ARG%PQ3(ARG%NO), ARG%PQ4(ARG%NO))
   ALLOCATE ( ARG%PQ5(ARG%NO), ARG%PQ6(ARG%NO), ARG%PQ7(ARG%NO), ARG%PQ8(ARG%NO))

! ---
! INITIALISE QUALITY FLAG
   ARG%FLG(:) = 1
   ARG%FLC(:) = 1

       READ (511)                                               &
        ARG%INO(1:ARG%NO), ARG%FLG(1:ARG%NO), ARG%PAR(1:ARG%NO) &
       ,ARG%LON(1:ARG%NO), ARG%LAT(1:ARG%NO)                    &
       ,ARG%DPT(1:ARG%NO), ARG%TIM(1:ARG%NO)                    &
       ,ARG%VAL(1:ARG%NO), ARG%BAC(1:ARG%NO)                    &
       ,ARG%ERR(1:ARG%NO), ARG%RES(1:ARG%NO)                    &
       ,ARG%IB(1:ARG%NO), ARG%JB(1:ARG%NO), ARG%KB(1:ARG%NO)    &
       ,ARG%PB(1:ARG%NO), ARG%QB(1:ARG%NO), ARG%RB(1:ARG%NO)
    CLOSE(511)

!      ARG%ERR(1:ARG%NO) = 10000.

! ---
! VERTICAL INTERPOLATION PARAMETERS
    DO K = 1,ARG%NO
     WRITE(749,*) ARG%LON(K),ARG%LAT(K),ARG%TIM(K),ARG%VAL(K)
     IF(ARG%FLG(K).EQ.1)THEN
       ARG%KB(K) = GRD%KM-1
     DO KK = 1,GRD%KM-1
      IF( ARG%DPT(K).GE.GRD%DEP(KK) .AND. ARG%DPT(K).LT.GRD%DEP(KK+1) ) THEN
       ARG%KB(K) = KK
       ARG%RB(K) = (ARG%DPT(K) - GRD%DEP(KK)) / (GRD%DEP(KK+1) - GRD%DEP(KK))
      ENDIF
     ENDDO
     ENDIF
    ENDDO


! RESIDUAL CHECK
!  DO K=1,ARG%NO
!   IF(ARG%PAR(K).EQ.1 .AND. ABS(ARG%RES(K)).GT.5.0) ARG%FLG(K) = 0
!   IF(ARG%PAR(K).EQ.2 .AND. ABS(ARG%RES(K)).GT.2.0) ARG%FLG(K) = 0
!  ENDDO


! ---
! COUNT GOOD OBSERVATIONS
    ARG%NC = 0
  DO K=1,ARG%NO
   WRITE(156,*) K, ARG%FLG(K), ARG%RES(K)
   IF(ARG%FLG(K).EQ.1)THEN
    ARG%NC = ARG%NC + 1
   ELSE
    ARG%FLC(K) = 0
    ARG%BIA(K) = 0.
    ARG%RES(K) = 0.
    ARG%INC(K) = 0.
    ARG%B_A(K) = 0.
    ARG%PQ1(K) = 0.
    ARG%PQ2(K) = 0.
    ARG%PQ3(K) = 0.
    ARG%PQ4(K) = 0.
    ARG%PQ5(K) = 0.
    ARG%PQ6(K) = 0.
    ARG%PQ7(K) = 0.
    ARG%PQ8(K) = 0.
   ENDIF
  ENDDO

  ARG%FLC(:) = ARG%FLG(:)

  WRITE(IOUNLOG,*) ' TOTAL NUMBER OF ARGO OBS          :',ARG%NO
  WRITE(IOUNLOG,*) ' TOTAL NUMBER OF RETAINED ARGO OBS :',ARG%NC

END SUBROUTINE GET_OBS_ARG



SUBROUTINE INT_PAR_ARG

!-----------------------------------------------------------------------
!                                                                      !
! GET INTERPOLATION PARAMETERS FOR A GRID                              !
!                                                                      !
! VERSION 1: S.DOBRICIC 2006                                           !
!-----------------------------------------------------------------------

 USE SET_KND
 USE GRD_STR
 USE OBS_STR
 USE IOUNITS, ONLY : IOUNOUT, IOUNLOG, IOUNERR

 IMPLICIT NONE

  INTEGER(I4)   ::  I, K,II, JJ,IM,JM,MSCOUNT
  INTEGER(I4)   ::  I1, J1, K1, IDEP,PRECI1,PRECJ1
  REAL(R8)      ::  P1, Q1, R1,PRECLAT,PRECLON,PRECP1,PRECQ1
  REAL(R8)      ::  MSK4, DIV_X, DIV_Y

  IF(ARG%NO.LE.0) THEN
        WRITE(IOUNLOG,*) '*** WARNING: ROUTINE INT_PAR_ARG: NO ARGO OBS!'
        RETURN
  ENDIF

  MSCOUNT = 0
  ARG%FLC(:) = ARG%FLG(:)

  IM=GRD%IM
  JM=GRD%JM
  PRECLON=999.
  PRECLAT=999.
  PRECI1=0
  PRECJ1=0
  PRECP1=0
  PRECQ1=0

  OBS_LOOP1 : DO K = 1,ARG%NO

  ! CHECK FLAG
      IF (ARG%FLC(K).EQ.0) CYCLE OBS_LOOP1
  ! CHECK CORRECTNESS OF LATLON VALUES
        IF ( ARG%LON(K).LT.-180._R8 .OR. ARG%LON(K).GT.180._R8  .OR. &
           & ARG%LAT(K).LT.-90._R8  .AND. ARG%LAT(K).LT. 90._R8 ) THEN
             WRITE (IOUNERR,*) 'NOT VALID COORDINATE FOR ARGO OBS'
             WRITE (IOUNERR,*) 'K, LON, LAT', K, ARG%LON(K), ARG%LAT(K)
             CALL ABOR1('WRONG COORDINATES IN OBS FILES')
        ENDIF
  ! CHECK IF ALREADY COMPUTED
        IF(K.GT.1) THEN
            IF(ARG%LON(K).EQ. ARG%LON(K-1) .AND. &
             & ARG%LAT(K).EQ. ARG%LAT(K-1) ) THEN

               ARG%FLC(K) = 1
               ARG%IB(K) = ARG%IB(K-1)
               ARG%JB(K) = ARG%JB(K-1)
               ARG%PB(K) = ARG%PB(K-1)
               ARG%QB(K) = ARG%QB(K-1)
               CYCLE OBS_LOOP1
            ENDIF
          ENDIF
          OUTER_LOOP : DO JJ=1,JM-1
           DO II=1,IM-1
             IF ( ARG%LON(K).GE.GRD%LON(II,JJ)      .AND. &
                & ARG%LON(K).LT.GRD%LON(II+1,JJ)    .AND. &
                & ARG%LAT(K).GE.GRD%LAT(II,JJ)      .AND. &
                & ARG%LAT(K).LT.GRD%LAT(II,JJ+1)    ) THEN

                  ARG%FLC(K) = 1

                  ARG%JB(K) = JJ
                  ARG%IB(K) = II
                  ARG%QB(K) = (ARG%LAT(K)-GRD%LAT(II,JJ))/&
                    & (GRD%LAT(II,JJ+1)-GRD%LAT(II,JJ))
                  ARG%PB(K) = (ARG%LON(K)-GRD%LON(II,JJ))/&
                    & (GRD%LON(II+1,JJ)-GRD%LON(II,JJ))

                  ! EXIT THE LOOP
#ifndef CRITICAL_LOOP_DEPTH
                  CYCLE OBS_LOOP1
#else
                  EXIT OUTER_LOOP
#endif
             END IF
          END DO
        END DO OUTER_LOOP

    ENDDO OBS_LOOP1

! ---
! UNDEFINE MASKED FOR MULTIGRID

    OBS_LOOP2 :  DO K = 1,ARG%NO

      IF(ARG%FLC(K).EQ.0) CYCLE OBS_LOOP2

      I1 = ARG%IB(K)
      J1 = ARG%JB(K)
      IDEP = ARG%KB(K)+1
      MSK4 = GRD%MSK(I1,J1,IDEP) + GRD%MSK(I1+1,J1,IDEP) + &
           & GRD%MSK(I1,J1+1,IDEP) + GRD%MSK(I1+1,J1+1,IDEP)

      IF(MSK4.LT.1.) THEN
        WRITE(135,*) K,I1,J1,IDEP,GRD%MSK(I1,J1,IDEP),GRD%MSK(I1+1,J1,IDEP),&
                   & GRD%MSK(I1,J1+1,IDEP), GRD%MSK(I1+1,J1+1,IDEP)
        WRITE(136,*) K,ARG%LON(K),ARG%LAT(K),&
        & GRD%LON(ARG%IB(K),ARG%JB(K)),GRD%LAT(ARG%IB(K),ARG%JB(K))
        WRITE(137,*) ARG%LON(K),ARG%LAT(K),1
        ARG%FLC(K) = 0
        MSCOUNT = MSCOUNT + 1
      ENDIF

    ENDDO OBS_LOOP2

    WRITE(IOUNLOG,*) &
    & ' ARGO OBS FILTERED OUT FOR MASKING INCONSISTENCIES: ',MSCOUNT

! ---
! HORIZONTAL INTERPOLATION PARAMETERS FOR EACH MASKED GRID

     OBS_LOOP3 :  DO K = 1,ARG%NO

         IF(ARG%FLC(K) .EQ. 0) CYCLE OBS_LOOP3

         I1=ARG%IB(K)
         P1=ARG%PB(K)
         J1=ARG%JB(K)
         Q1=ARG%QB(K)

         K1=ARG%KB(K)
         DIV_Y =  (1.-Q1) * MAX(GRD%MSK(I1,J1  ,K1),GRD%MSK(I1+1,J1  ,K1))     &
                 +    Q1  * MAX(GRD%MSK(I1,J1+1,K1),GRD%MSK(I1+1,J1+1,K1))
         DIV_X =  (1.-P1) * GRD%MSK(I1  ,J1,K1) + P1 * GRD%MSK(I1+1,J1,K1)
          ARG%PQ1(K) = GRD%MSK(I1,J1,K1)                                       &
                      * MAX(GRD%MSK(I1,J1,K1),GRD%MSK(I1+1,J1,K1))             &
                       * (1.-P1) * (1.-Q1)                                     &
                      /( DIV_X * DIV_Y + 1.E-16 )
          ARG%PQ2(K) = GRD%MSK(I1+1,J1,K1)                                     &
                      * MAX(GRD%MSK(I1,J1,K1),GRD%MSK(I1+1,J1,K1))             &
                      *     P1  * (1.-Q1)                                      &
                      /( DIV_X * DIV_Y + 1.E-16 )
         DIV_X =  (1.-P1) * GRD%MSK(I1  ,J1+1,K1) + P1 * GRD%MSK(I1+1,J1+1,K1)
          ARG%PQ3(K) = GRD%MSK(I1,J1+1,K1)                                     &
                      * MAX(GRD%MSK(I1,J1+1,K1),GRD%MSK(I1+1,J1+1,K1))         &
                      * (1.-P1) *     Q1                                       &
                      /( DIV_X * DIV_Y + 1.E-16 )
          ARG%PQ4(K) = GRD%MSK(I1+1,J1+1,K1)                                   &
                      * MAX(GRD%MSK(I1,J1+1,K1),GRD%MSK(I1+1,J1+1,K1))         &
                      *     P1  *     Q1                                       &
                      /( DIV_X * DIV_Y + 1.E-16 )

         K1=ARG%KB(K) + 1
         DIV_Y =  (1.-Q1) * MAX(GRD%MSK(I1,J1  ,K1),GRD%MSK(I1+1,J1  ,K1))     &
                 +    Q1  * MAX(GRD%MSK(I1,J1+1,K1),GRD%MSK(I1+1,J1+1,K1))
         DIV_X =  (1.-P1) * GRD%MSK(I1  ,J1,K1) + P1 * GRD%MSK(I1+1,J1,K1)
          ARG%PQ5(K) = GRD%MSK(I1,J1,K1)                                       &
                      * MAX(GRD%MSK(I1,J1,K1),GRD%MSK(I1+1,J1,K1))             &
                       * (1.-P1) * (1.-Q1)                                     &
                      /( DIV_X * DIV_Y + 1.E-16 )
          ARG%PQ6(K) = GRD%MSK(I1+1,J1,K1)                                     &
                      * MAX(GRD%MSK(I1,J1,K1),GRD%MSK(I1+1,J1,K1))             &
                      *     P1  * (1.-Q1)                                      &
                      /( DIV_X * DIV_Y + 1.E-16 )
         DIV_X =  (1.-P1) * GRD%MSK(I1  ,J1+1,K1) + P1 * GRD%MSK(I1+1,J1+1,K1)
          ARG%PQ7(K) = GRD%MSK(I1,J1+1,K1)                                     &
                      * MAX(GRD%MSK(I1,J1+1,K1),GRD%MSK(I1+1,J1+1,K1))         &
                      * (1.-P1) *     Q1                                       &
                      /( DIV_X * DIV_Y + 1.E-16 )
          ARG%PQ8(K) = GRD%MSK(I1+1,J1+1,K1)                                   &
                      * MAX(GRD%MSK(I1,J1+1,K1),GRD%MSK(I1+1,J1+1,K1))         &
                      *     P1  *     Q1                                       &
                      /( DIV_X * DIV_Y + 1.E-16 )

         R1=ARG%RB(K)

         ARG%PQ1(K) = (1.-R1) * ARG%PQ1(K)
         ARG%PQ2(K) = (1.-R1) * ARG%PQ2(K)
         ARG%PQ3(K) = (1.-R1) * ARG%PQ3(K)
         ARG%PQ4(K) = (1.-R1) * ARG%PQ4(K)
         ARG%PQ5(K) =     R1  * ARG%PQ5(K)
         ARG%PQ6(K) =     R1  * ARG%PQ6(K)
         ARG%PQ7(K) =     R1  * ARG%PQ7(K)
         ARG%PQ8(K) =     R1  * ARG%PQ8(K)

       ENDDO OBS_LOOP3

! ---
! COUNT GOOD OBSERVATIONS
    ARG%NC = 0
  DO K=1,ARG%NO
   IF(ARG%FLC(K).EQ.1)THEN
    ARG%NC = ARG%NC + 1
   ENDIF
  ENDDO

  WRITE(IOUNLOG,*)
  WRITE(IOUNLOG,*) ' *** ARGO OBS AFTER INTERPOLATION SET-UP'
  WRITE(IOUNLOG,*) ' TOTAL NUMBER OF ARGO OBS          :',ARG%NO
  WRITE(IOUNLOG,*) ' TOTAL NUMBER OF RETAINED ARGO OBS :',ARG%NC
  WRITE(IOUNLOG,*) ' TOTAL NUMBER OF RETAINED ARGO OBS :',SUM(ARG%FLC)

END SUBROUTINE INT_PAR_ARG
