SUBROUTINE GET_OBS_GLD

!-----------------------------------------------------------------------
!                                                                      !
! LOAD GLIDER OBSERVATIONS                                                !
!                                                                      !
! VERSION 1: S.DOBRICIC 2007                                           !
!-----------------------------------------------------------------------

 USE SET_KND
 USE GRD_STR
 USE OBS_STR

 IMPLICIT NONE

  INTEGER(I4)   ::  K
  INTEGER(I4)   ::  I1, KK, I


    GLD%NO = 0
    GLD%NC = 0

  IF(OBS%GLD.EQ.0)RETURN

  OPEN(511,FILE='GLD_MIS.DAT',FORM='UNFORMATTED',STATUS='OLD')

! ---
! ALLOCATE MEMORY FOR OBSERVATIONS

   READ(511) GLD%NO

   PRINT*,'NUMBER OF GLIDER OBSERVATIONS: ',GLD%NO

   IF(GLD%NO.EQ.0)RETURN

   ALLOCATE ( GLD%INO(GLD%NO), GLD%FLG(GLD%NO), GLD%FLC(GLD%NO), GLD%PAR(GLD%NO))
   ALLOCATE ( GLD%LON(GLD%NO), GLD%LAT(GLD%NO), GLD%DPT(GLD%NO), GLD%TIM(GLD%NO))
   ALLOCATE ( GLD%VAL(GLD%NO), GLD%BAC(GLD%NO), GLD%INC(GLD%NO))
   ALLOCATE ( GLD%BIA(GLD%NO), GLD%ERR(GLD%NO))
   ALLOCATE ( GLD%RES(GLD%NO), GLD%B_A(GLD%NO))
   ALLOCATE ( GLD%IB(GLD%NO), GLD%JB(GLD%NO), GLD%KB(GLD%NO))
   ALLOCATE ( GLD%PB(GLD%NO), GLD%QB(GLD%NO), GLD%RB(GLD%NO))
   ALLOCATE ( GLD%PQ1(GLD%NO), GLD%PQ2(GLD%NO), GLD%PQ3(GLD%NO), GLD%PQ4(GLD%NO))
   ALLOCATE ( GLD%PQ5(GLD%NO), GLD%PQ6(GLD%NO), GLD%PQ7(GLD%NO), GLD%PQ8(GLD%NO))

! ---
! INITIALISE QUALITY FLAG
   GLD%FLG(:) = 1
   GLD%FLC(:) = 1

       READ (511)                                               &
        GLD%INO(1:GLD%NO), GLD%FLG(1:GLD%NO), GLD%PAR(1:GLD%NO) &
       ,GLD%LON(1:GLD%NO), GLD%LAT(1:GLD%NO)                    &
       ,GLD%DPT(1:GLD%NO), GLD%TIM(1:GLD%NO)                    &
       ,GLD%VAL(1:GLD%NO), GLD%BAC(1:GLD%NO)                    &
       ,GLD%ERR(1:GLD%NO), GLD%RES(1:GLD%NO)                    &
       ,GLD%IB(1:GLD%NO), GLD%JB(1:GLD%NO), GLD%KB(1:GLD%NO)    &
       ,GLD%PB(1:GLD%NO), GLD%QB(1:GLD%NO), GLD%RB(1:GLD%NO)
    CLOSE(511)

!      GLD%ERR(1:GLD%NO) = 10000.

! ---
! VERTICAL INTERPOLATION PARAMETERS
    DO K = 1,GLD%NO
     IF(GLD%FLG(K).EQ.1)THEN
       GLD%KB(K) = GRD%KM-1
     DO KK = 1,GRD%KM-1
      IF( GLD%DPT(K).GE.GRD%DEP(KK) .AND. GLD%DPT(K).LT.GRD%DEP(KK+1) ) THEN
       GLD%KB(K) = KK
       GLD%RB(K) = (GLD%DPT(K) - GRD%DEP(KK)) / (GRD%DEP(KK+1) - GRD%DEP(KK))
      ENDIF
     ENDDO
     ENDIF
    ENDDO


! RESIDUAL CHECK
  DO K=1,GLD%NO
   IF(GLD%PAR(K).EQ.1 .AND. ABS(GLD%RES(K)).GT.5.0) GLD%FLG(K) = 0
   IF(GLD%PAR(K).EQ.2 .AND. ABS(GLD%RES(K)).GT.2.0) GLD%FLG(K) = 0
  ENDDO


! ---
! COUNT GOOD OBSERVATIONS
    GLD%NC = 0
  DO K=1,GLD%NO
   IF(GLD%FLG(K).EQ.1)THEN
    GLD%NC = GLD%NC + 1
   ELSE
    GLD%FLC(K) = 0
    GLD%BIA(K) = 0.
    GLD%RES(K) = 0.
    GLD%INC(K) = 0.
    GLD%B_A(K) = 0.
    GLD%PQ1(K) = 0.
    GLD%PQ2(K) = 0.
    GLD%PQ3(K) = 0.
    GLD%PQ4(K) = 0.
    GLD%PQ5(K) = 0.
    GLD%PQ6(K) = 0.
    GLD%PQ7(K) = 0.
    GLD%PQ8(K) = 0.
   ENDIF
  ENDDO

  GLD%FLC(:) = GLD%FLG(:)

   PRINT*,'NUMBER OF GOOD GLIDER OBSERVATIONS: ',GLD%NC

END SUBROUTINE GET_OBS_GLD



SUBROUTINE INT_PAR_GLD

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

  INTEGER(I4)   ::  I, K
  INTEGER(I4)   ::  I1, J1, K1, IDEP
  REAL(R8)      ::  P1, Q1, R1
  REAL(R8)      ::  MSK4, DIV_X, DIV_Y
  LOGICAL       ::  INSL

  INSL(I,I1) = I.GE.1 .AND. I.LT.I1

 IF(GLD%NO.GT.0) THEN

   GLD%FLC(:) = GLD%FLG(:)

! ---
! HORIZONTAL INTERPOLATION PARAMETERS
    DO K = 1,GLD%NO
     Q1 = (GLD%LAT(K) - GRD%LAT(1,1)) / GRD%DLT + 1.0
     J1 = INT(Q1)
     P1 = (GLD%LON(K) - GRD%LON(1,1)) / GRD%DLN + 1.0
     I1 = INT(P1)
     IF(INSL(J1,GRD%JM) .AND. INSL(I1,GRD%IM)) THEN
       GLD%IB(K) = I1
       GLD%JB(K) = J1
       GLD%PB(K) = (P1-I1)
       GLD%QB(K) = (Q1-J1)
     ELSE
       GLD%FLC(K) = 0
     ENDIF
    ENDDO

! ---
! UNDEFINE MASKED FOR MULTIGRID
    DO K = 1,GLD%NO
     IF(GLD%FLC(K).EQ.1)THEN
      I1 = GLD%IB(K)
      J1 = GLD%JB(K)
      IDEP = GLD%KB(K)+1
      MSK4 = GRD%MSK(I1,J1,IDEP) + GRD%MSK(I1+1,J1,IDEP) + GRD%MSK(I1,J1+1,IDEP) + GRD%MSK(I1+1,J1+1,IDEP)
      IF(MSK4.LT.1.) GLD%FLC(K) = 0
     ENDIF
    ENDDO

! ---
! HORIZONTAL INTERPOLATION PARAMETERS FOR EACH MASKED GRID
       DO K = 1,GLD%NO
        IF(GLD%FLC(K) .EQ. 1) THEN

         I1=GLD%IB(K)
         P1=GLD%PB(K)
         J1=GLD%JB(K)
         Q1=GLD%QB(K)


         K1=GLD%KB(K)
         DIV_Y =  (1.-Q1) * MAX(GRD%MSK(I1,J1  ,K1),GRD%MSK(I1+1,J1  ,K1))     &
                 +    Q1  * MAX(GRD%MSK(I1,J1+1,K1),GRD%MSK(I1+1,J1+1,K1))
         DIV_X =  (1.-P1) * GRD%MSK(I1  ,J1,K1) + P1 * GRD%MSK(I1+1,J1,K1)
          GLD%PQ1(K) = GRD%MSK(I1,J1,K1)                                       &
                      * MAX(GRD%MSK(I1,J1,K1),GRD%MSK(I1+1,J1,K1))             &
                       * (1.-P1) * (1.-Q1)                                     &
                      /( DIV_X * DIV_Y + 1.E-16 )
          GLD%PQ2(K) = GRD%MSK(I1+1,J1,K1)                                     &
                      * MAX(GRD%MSK(I1,J1,K1),GRD%MSK(I1+1,J1,K1))             &
                      *     P1  * (1.-Q1)                                      &
                      /( DIV_X * DIV_Y + 1.E-16 )
         DIV_X =  (1.-P1) * GRD%MSK(I1  ,J1+1,K1) + P1 * GRD%MSK(I1+1,J1+1,K1)
          GLD%PQ3(K) = GRD%MSK(I1,J1+1,K1)                                     &
                      * MAX(GRD%MSK(I1,J1+1,K1),GRD%MSK(I1+1,J1+1,K1))         &
                      * (1.-P1) *     Q1                                       &
                      /( DIV_X * DIV_Y + 1.E-16 )
          GLD%PQ4(K) = GRD%MSK(I1+1,J1+1,K1)                                   &
                      * MAX(GRD%MSK(I1,J1+1,K1),GRD%MSK(I1+1,J1+1,K1))         &
                      *     P1  *     Q1                                       &
                      /( DIV_X * DIV_Y + 1.E-16 )

         K1=GLD%KB(K) + 1
         DIV_Y =  (1.-Q1) * MAX(GRD%MSK(I1,J1  ,K1),GRD%MSK(I1+1,J1  ,K1))     &
                 +    Q1  * MAX(GRD%MSK(I1,J1+1,K1),GRD%MSK(I1+1,J1+1,K1))
         DIV_X =  (1.-P1) * GRD%MSK(I1  ,J1,K1) + P1 * GRD%MSK(I1+1,J1,K1)
          GLD%PQ5(K) = GRD%MSK(I1,J1,K1)                                       &
                      * MAX(GRD%MSK(I1,J1,K1),GRD%MSK(I1+1,J1,K1))             &
                       * (1.-P1) * (1.-Q1)                                     &
                      /( DIV_X * DIV_Y + 1.E-16 )
          GLD%PQ6(K) = GRD%MSK(I1+1,J1,K1)                                     &
                      * MAX(GRD%MSK(I1,J1,K1),GRD%MSK(I1+1,J1,K1))             &
                      *     P1  * (1.-Q1)                                      &
                      /( DIV_X * DIV_Y + 1.E-16 )
         DIV_X =  (1.-P1) * GRD%MSK(I1  ,J1+1,K1) + P1 * GRD%MSK(I1+1,J1+1,K1)
          GLD%PQ7(K) = GRD%MSK(I1,J1+1,K1)                                     &
                      * MAX(GRD%MSK(I1,J1+1,K1),GRD%MSK(I1+1,J1+1,K1))         &
                      * (1.-P1) *     Q1                                       &
                      /( DIV_X * DIV_Y + 1.E-16 )
          GLD%PQ8(K) = GRD%MSK(I1+1,J1+1,K1)                                   &
                      * MAX(GRD%MSK(I1,J1+1,K1),GRD%MSK(I1+1,J1+1,K1))         &
                      *     P1  *     Q1                                       &
                      /( DIV_X * DIV_Y + 1.E-16 )

         R1=GLD%RB(K)
          GLD%PQ1(K) = (1.-R1) * GLD%PQ1(K)
          GLD%PQ2(K) = (1.-R1) * GLD%PQ2(K)
          GLD%PQ3(K) = (1.-R1) * GLD%PQ3(K)
          GLD%PQ4(K) = (1.-R1) * GLD%PQ4(K)
          GLD%PQ5(K) =     R1  * GLD%PQ5(K)
          GLD%PQ6(K) =     R1  * GLD%PQ6(K)
          GLD%PQ7(K) =     R1  * GLD%PQ7(K)
          GLD%PQ8(K) =     R1  * GLD%PQ8(K)

        ENDIF
       ENDDO


! ---
! COUNT GOOD OBSERVATIONS
    GLD%NC = 0
  DO K=1,GLD%NO
   IF(GLD%FLC(K).EQ.1)THEN
    GLD%NC = GLD%NC + 1
   ENDIF
  ENDDO

 ENDIF

END SUBROUTINE INT_PAR_GLD
