SUBROUTINE GET_OBS_TRJ

!-----------------------------------------------------------------------
!                                                                      !
! LOAD SLA OBSERVATIONS                                                !
!                                                                      !
! VERSION 1: S.DOBRICIC 2006                                           !
!-----------------------------------------------------------------------

 USE SET_KND
 USE GRD_STR
 USE OBS_STR

 IMPLICIT NONE

  INTEGER(I4)   ::  K
  INTEGER(I4)   ::  I1, J1, KK, I, J, ITER
  REAL(R8)      ::  P1, Q1, P, Q, SUMT, SUMI, TIMP, DUM1, DUM2

   TRJ%NO = 0
   TRJ%NC = 0

  IF(OBS%TRJ.EQ.0)RETURN

! ---
! ALLOCATE MEMORY FOR OBSERVATIONS


    TRJ%IMB = 871
    TRJ%JMB = 253

    TRJ%JPT = 20
    TRJ%JPN = 4

    TRJ%NO = 2*TRJ%JPN

!   READ(511) TRJ%NO
   PRINT*,'NUMBER OF TRJOCITY OBSERVATIONS: ',  TRJ%NO


   IF(TRJ%NO.EQ.0)RETURN

   ALLOCATE ( TRJ%INO(TRJ%NO), TRJ%FLG(TRJ%NO), TRJ%FLC(TRJ%NO), TRJ%PAR(TRJ%NO))
   ALLOCATE ( TRJ%LON(TRJ%NO), TRJ%LAT(TRJ%NO), TRJ%DPH(TRJ%NO), TRJ%TIM(TRJ%NO))
   ALLOCATE ( TRJ%VAL(TRJ%NO), TRJ%BAC(TRJ%NO), TRJ%INC(TRJ%NO))
   ALLOCATE ( TRJ%BIA(TRJ%NO), TRJ%ERR(TRJ%NO))
   ALLOCATE ( TRJ%RES(TRJ%NO), TRJ%B_A(TRJ%NO))
   ALLOCATE ( TRJ%IB(TRJ%NO), TRJ%JB(TRJ%NO))
   ALLOCATE ( TRJ%PB(TRJ%NO), TRJ%QB(TRJ%NO))
   ALLOCATE ( TRJ%PQ1(TRJ%NO), TRJ%PQ2(TRJ%NO), TRJ%PQ3(TRJ%NO), TRJ%PQ4(TRJ%NO))
   ALLOCATE ( TRJ%DPT(TRJ%NO))
   ALLOCATE ( TRJ%UMN(TRJ%IMB,TRJ%JMB), TRJ%VMN(TRJ%IMB,TRJ%JMB) )

   ALLOCATE ( TRJ%XMN(TRJ%JPT+1,TRJ%JPN), TRJ%YMN(TRJ%JPT+1,TRJ%JPN) )
   ALLOCATE ( TRJ%TIM(TRJ%JPN) )
   ALLOCATE ( TRJ%XOB(TRJ%JPN), TRJ%YOB(TRJ%JPN) )
   ALLOCATE ( TRJ%XTL(TRJ%JPN), TRJ%YTL(TRJ%JPN) )
   ALLOCATE ( TRJ%XTL_AD(TRJ%JPN), TRJ%YTL_AD(TRJ%JPN) )


! ---
! INITIALISE QUALITY FLAG
   TRJ%FLG(:) = 1
   TRJ%FLC(:) = 1

! ---
! LEVEL CORRESPONDING TO THE MINIMUM DEPTH
   TRJ%KDP=GRD%KM
  DO K=GRD%KM, 1, -1
   IF(GRD%DEP(K).GE.TRJ%DEP) TRJ%KDP = K
  ENDDO

      TRJ%KDP =  33

   PRINT*,' TRJ%DEP: ',TRJ%KDP

    TRJ%UMN(:,:) = 0.0
    TRJ%VMN(:,:) = 0.0

!    DO J=137,137+92-1
!    DO I=295,295+170-1
!     READ(511,*) TRJ%UMN(I,J), TRJ%VMN(I,J), DUM1, DUM2
!    ENDDO
!    ENDDO
    OPEN(511,FILE='VEL_MEAN.DAT',FORM='FORMATTED',STATUS='OLD')
     READ(511) TRJ%UMN, TRJ%VMN
    CLOSE(511)
    TRJ%UMN(:,:) = TRJ%UMN(:,:)
    TRJ%VMN(:,:) = TRJ%VMN(:,:)

    OPEN(511,FILE='TRJ_BCK.DAT',FORM='FORMATTED',STATUS='OLD')
    DO I=1,TRJ%JPT+1
    DO J=1,TRJ%JPN
     READ(511,*) TRJ%XMN(I,J), TRJ%YMN(I,J)
     TRJ%XMN(I,J) = TRJ%XMN(I,J) !+ 294.
     TRJ%YMN(I,J) = TRJ%YMN(I,J) !+ 136.
    ENDDO
    ENDDO
    CLOSE(511)

    OPEN(511,FILE='OBSERVATIONS.DAT',FORM='FORMATTED',STATUS='OLD')
    DO J=1,TRJ%JPN
     READ(511,*) DUM1, DUM2, TRJ%XOB(J), TRJ%YOB(J), TRJ%TIM(J)
     TRJ%XOB(J) = (TRJ%XOB(J)+18.125)/0.0625 + 1.
     TRJ%YOB(J) = (TRJ%YOB(J)-30.250)/0.0625 + 1.
    ENDDO
    CLOSE(511)


    TRJ%ERR(:) = 0.1


    DO J=1,TRJ%JPN
     TRJ%RES(        J) = TRJ%XOB(J) - TRJ%XMN(TRJ%JPT+1,J)
     TRJ%PAR(        J) = 1
     TRJ%RES(TRJ%JPN+J) = TRJ%YOB(J) - TRJ%YMN(TRJ%JPT+1,J)
     TRJ%PAR(TRJ%JPN+J) = 2
    ENDDO
     PRINT*,'-------------------------------'
       PRINT*,TRJ%RES(:)
     PRINT*,'-------------------------------'
!     TRJ%RES(2:4) = 0.0
!     TRJ%RES(6:8) = 0.0


! RESIDUAL CHECK
  DO K=1,TRJ%NO
!   IF(ABS(TRJ%RES(K)).GT.3.) TRJ%FLG(K) = 0
  ENDDO

! ---
! COUNT GOOD OBSERVATIONS
    TRJ%NC = 0
  DO K=1,TRJ%NO
   IF(TRJ%FLG(K).EQ.1)THEN
    TRJ%NC = TRJ%NC + 1
   ELSE
    TRJ%FLC(K) = 0
    TRJ%BIA(K) = 0.
    TRJ%RES(K) = 0.
    TRJ%INC(K) = 0.
    TRJ%B_A(K) = 0.
    TRJ%PQ1(K) = 0.
    TRJ%PQ2(K) = 0.
    TRJ%PQ3(K) = 0.
    TRJ%PQ4(K) = 0.
   ENDIF
  ENDDO

  TRJ%FLC(:) = TRJ%FLG(:)

  PRINT*,'QUALITY FLAGS: ',TRJ%FLC(:)
  PRINT*,'RESIDUALS: ',TRJ%RES(:)

END SUBROUTINE GET_OBS_TRJ

SUBROUTINE INT_PAR_TRJ

!-----------------------------------------------------------------------
!                                                                      !
! GET INTERPOLATION PARAMETERS FOR A GRID                              !
!                                                                      !
! VERSION 1: S.DOBRICIC 2006                                           !
!-----------------------------------------------------------------------

 USE SET_KND
 USE GRD_STR
 USE OBS_STR

 IMPLICIT NONE

  INTEGER(I4)   ::  K, N_TOPEX, N_ERS
  INTEGER(I4)   ::  I1, J1, KK, I, J, ITER
  REAL(R8)      ::  P1, Q1, P, Q, SUMT, SUMI, TIMP
  REAL(R8)      ::  MSK4, DIV_X, DIV_Y
  LOGICAL       ::  INSL

  INSL(I,I1) = I.GE.1 .AND. I.LT.I1

  IF(TRJ%NO.GT.0) THEN

       TRJ%FLC(:) = TRJ%FLG(:)

! ---
! INTERPOLATION PARAMETERS
    DO KK = 1,TRJ%NO
     Q1 = (TRJ%LAT(KK) - GRD%LAT(1,1)) / GRD%DLT + 1.0
     J1 = INT(Q1)
     P1 = (TRJ%LON(KK) - GRD%LON(1,1)) / GRD%DLN + 1.0
     I1 = INT(P1)
     PRINT*,'---> ',I1,J1,KK,TRJ%LAT(KK),TRJ%LON(KK)
     IF(INSL(J1,GRD%JM) .AND. INSL(I1,GRD%IM)) THEN
       MSK4 = GRD%MSK(I1,J1,TRJ%KDP) + GRD%MSK(I1+1,J1,TRJ%KDP) + GRD%MSK(I1,J1+1,TRJ%KDP) + GRD%MSK(I1+1,J1+1,TRJ%KDP)
      IF(MSK4.GE.1.0)THEN
       TRJ%IB(KK) = I1
       TRJ%JB(KK) = J1
       TRJ%PB(KK) = (P1-I1)
       TRJ%QB(KK) = (Q1-J1)
      ELSE
       TRJ%FLC(KK) = 0
      ENDIF
     ELSE
       TRJ%FLC(KK) = 0
     ENDIF
    ENDDO


! ---
! HORIZONTAL INTERPOLATION PARAMETERS FOR EACH MASKED GRID
       DO K = 1,TRJ%NO
          PRINT*,TRJ%FLC(K)
        IF(TRJ%FLC(K) .EQ. 1) THEN

         I1=TRJ%IB(K)
         P1=TRJ%PB(K)
         J1=TRJ%JB(K)
         Q1=TRJ%QB(K)

         DIV_Y =  (1.-Q1) * MAX(GRD%MSK(I1,J1  ,TRJ%KDP),GRD%MSK(I1+1,J1  ,TRJ%KDP))     &
                 +    Q1  * MAX(GRD%MSK(I1,J1+1,TRJ%KDP),GRD%MSK(I1+1,J1+1,TRJ%KDP))
         DIV_X =  (1.-P1) * GRD%MSK(I1  ,J1,TRJ%KDP) + P1 * GRD%MSK(I1+1,J1,TRJ%KDP)
          TRJ%PQ1(K) = GRD%MSK(I1,J1,TRJ%KDP)                                      &
                      * MAX(GRD%MSK(I1,J1,TRJ%KDP),GRD%MSK(I1+1,J1,TRJ%KDP))             &
                       * (1.-P1) * (1.-Q1)                                   &
                      /( DIV_X * DIV_Y + 1.E-16 )
          TRJ%PQ2(K) = GRD%MSK(I1+1,J1,TRJ%KDP)                                    &
                      * MAX(GRD%MSK(I1,J1,TRJ%KDP),GRD%MSK(I1+1,J1,TRJ%KDP))             &
                      *     P1  * (1.-Q1)                                    &
                      /( DIV_X * DIV_Y + 1.E-16 )
         DIV_X =  (1.-P1) * GRD%MSK(I1  ,J1+1,TRJ%KDP) + P1 * GRD%MSK(I1+1,J1+1,TRJ%KDP)
          TRJ%PQ3(K) = GRD%MSK(I1,J1+1,TRJ%KDP)                                    &
                      * MAX(GRD%MSK(I1,J1+1,TRJ%KDP),GRD%MSK(I1+1,J1+1,TRJ%KDP))         &
                      * (1.-P1) *     Q1                                     &
                      /( DIV_X * DIV_Y + 1.E-16 )
          TRJ%PQ4(K) = GRD%MSK(I1+1,J1+1,TRJ%KDP)                                  &
                      * MAX(GRD%MSK(I1,J1+1,TRJ%KDP),GRD%MSK(I1+1,J1+1,TRJ%KDP))         &
                      *     P1  *     Q1                                     &
                      /( DIV_X * DIV_Y + 1.E-16 )

        ENDIF
       ENDDO

! ---
! COUNT GOOD OBSERVATIONS
    TRJ%NC = 0
  DO K=1,TRJ%NO
   IF(TRJ%FLC(K).EQ.1)THEN
    TRJ%NC = TRJ%NC + 1
   ENDIF
  ENDDO


 ENDIF




END SUBROUTINE INT_PAR_TRJ
