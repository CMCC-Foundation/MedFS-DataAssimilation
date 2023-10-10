SUBROUTINE GET_OBS_VEL

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
  REAL(R8)      ::  P1, Q1, P, Q, SUMT, SUMI, TIMP
  REAL(R8), DIMENSION(871,253)      ::  UTMP, VTMP

   VEL%NO = 0
   VEL%NC = 0

  IF(OBS%VEL.EQ.0)RETURN

  OPEN(511,FILE='VEL_MIS.DAT',FORM='UNFORMATTED',STATUS='OLD')

! ---
! ALLOCATE MEMORY FOR OBSERVATIONS

    READ(511) UTMP,VTMP
    CLOSE(511)

   VEL%NO = 0
    DO J=1,253
    DO I=1,871
     IF(UTMP(I,J).NE.0.0)THEN
      VEL%NO = VEL%NO + 1
     ENDIF
     IF(VTMP(I,J).NE.0.0)THEN
      VEL%NO = VEL%NO + 1
     ENDIF
    ENDDO
    ENDDO
!   READ(511) VEL%NO
   VEL%NO = 0
   PRINT*,'NUMBER OF VELOCITY OBSERVATIONS: ',  VEL%NO


   IF(VEL%NO.EQ.0)RETURN

   ALLOCATE ( VEL%INO(VEL%NO), VEL%FLG(VEL%NO), VEL%FLC(VEL%NO), VEL%PAR(VEL%NO))
   ALLOCATE ( VEL%LON(VEL%NO), VEL%LAT(VEL%NO), VEL%DPH(VEL%NO), VEL%TIM(VEL%NO))
   ALLOCATE ( VEL%VAL(VEL%NO), VEL%BAC(VEL%NO), VEL%INC(VEL%NO))
   ALLOCATE ( VEL%BIA(VEL%NO), VEL%ERR(VEL%NO))
   ALLOCATE ( VEL%RES(VEL%NO), VEL%B_A(VEL%NO))
   ALLOCATE ( VEL%IB(VEL%NO), VEL%JB(VEL%NO))
   ALLOCATE ( VEL%PB(VEL%NO), VEL%QB(VEL%NO))
   ALLOCATE ( VEL%PQ1(VEL%NO), VEL%PQ2(VEL%NO), VEL%PQ3(VEL%NO), VEL%PQ4(VEL%NO))
   ALLOCATE ( VEL%DPT(VEL%NO))

   K = 0
    DO J=1,253
    DO I=1,871
     IF(UTMP(I,J).NE.0.0)THEN
      K = K + 1
      VEL%RES(K) = UTMP(I,J)
      VEL%PAR(K) = 1
      VEL%ERR(K) = 0.005
      VEL%LON(K) = -18.125 + (I-1)*0.0625
      VEL%LAT(K) =  30.250 + (J-1)*0.0625
      VEL%DPT(K) = 350.00
      VEL%FLG(K) = 1
     ENDIF
     IF(VTMP(I,J).NE.0.0)THEN
      K = K + 1
      VEL%RES(K) = VTMP(I,J)
      VEL%PAR(K) = 2
      VEL%ERR(K) = 0.005
      VEL%LON(K) = -18.125 + (I-1.)*0.0625
      VEL%LAT(K) =  30.250 + (J-1.)*0.0625
      VEL%DPT(K) = 350.00
      VEL%FLG(K) = 1
     ENDIF
    ENDDO
    ENDDO

! ---
! INITIALISE QUALITY FLAG
   VEL%FLG(:) = 1
   VEL%FLC(:) = 1

! ---
! LEVEL CORRESPONDING TO THE MINIMUM DEPTH
   VEL%KDP=GRD%KM
  DO K=GRD%KM, 1, -1
   IF(GRD%DEP(K).GE.VEL%DEP) VEL%KDP = K
  ENDDO
      VEL%KDP =  32

   PRINT*,' VEL%DEP: ',VEL%KDP

!       READ (511)                                              &
!        VEL%INO(1:VEL%NO), VEL%FLG(1:VEL%NO)                   &
!       ,VEL%LON(1:VEL%NO), VEL%LAT(1:VEL%NO), VEL%TIM(1:VEL%NO)&
!       ,VEL%VAL(1:VEL%NO), VEL%BAC(1:VEL%NO)                   &
!       ,VEL%ERR(1:VEL%NO), VEL%RES(1:VEL%NO)                   &
!       ,VEL%IB(1:VEL%NO), VEL%JB(1:VEL%NO)                     &
!       ,VEL%PB(1:VEL%NO), VEL%QB(1:VEL%NO)
!    CLOSE(511)

! RESIDUAL CHECK
  DO K=1,VEL%NO
   IF(ABS(VEL%RES(K)).GT.3.) VEL%FLG(K) = 0
  ENDDO

! ---
! COUNT GOOD OBSERVATIONS
    VEL%NC = 0
  DO K=1,VEL%NO
   IF(VEL%FLG(K).EQ.1)THEN
    VEL%NC = VEL%NC + 1
   ELSE
    VEL%FLC(K) = 0
    VEL%BIA(K) = 0.
    VEL%RES(K) = 0.
    VEL%INC(K) = 0.
    VEL%B_A(K) = 0.
    VEL%PQ1(K) = 0.
    VEL%PQ2(K) = 0.
    VEL%PQ3(K) = 0.
    VEL%PQ4(K) = 0.
   ENDIF
  ENDDO

  VEL%FLC(:) = VEL%FLG(:)

END SUBROUTINE GET_OBS_VEL

SUBROUTINE INT_PAR_VEL

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

 IF(VEL%NO.GT.0) THEN

       VEL%FLC(:) = VEL%FLG(:)

! ---
! INTERPOLATION PARAMETERS
    DO KK = 1,VEL%NO
     Q1 = (VEL%LAT(KK) - GRD%LAT(1,1)) / GRD%DLT + 1.0
     J1 = INT(Q1)
     P1 = (VEL%LON(KK) - GRD%LON(1,1)) / GRD%DLN + 1.0
     I1 = INT(P1)
     PRINT*,'---> ',I1,J1,KK,VEL%LAT(KK),VEL%LON(KK)
     IF(INSL(J1,GRD%JM) .AND. INSL(I1,GRD%IM)) THEN
       MSK4 = GRD%MSK(I1,J1,VEL%KDP) + GRD%MSK(I1+1,J1,VEL%KDP) + GRD%MSK(I1,J1+1,VEL%KDP) + GRD%MSK(I1+1,J1+1,VEL%KDP)
      IF(MSK4.GE.1.0)THEN
       VEL%IB(KK) = I1
       VEL%JB(KK) = J1
       VEL%PB(KK) = (P1-I1)
       VEL%QB(KK) = (Q1-J1)
      ELSE
       VEL%FLC(KK) = 0
      ENDIF
     ELSE
       VEL%FLC(KK) = 0
     ENDIF
    ENDDO


! ---
! HORIZONTAL INTERPOLATION PARAMETERS FOR EACH MASKED GRID
       DO K = 1,VEL%NO
          PRINT*,VEL%FLC(K)
        IF(VEL%FLC(K) .EQ. 1) THEN

         I1=VEL%IB(K)
         P1=VEL%PB(K)
         J1=VEL%JB(K)
         Q1=VEL%QB(K)

         DIV_Y =  (1.-Q1) * MAX(GRD%MSK(I1,J1  ,VEL%KDP),GRD%MSK(I1+1,J1  ,VEL%KDP))     &
                 +    Q1  * MAX(GRD%MSK(I1,J1+1,VEL%KDP),GRD%MSK(I1+1,J1+1,VEL%KDP))
         DIV_X =  (1.-P1) * GRD%MSK(I1  ,J1,VEL%KDP) + P1 * GRD%MSK(I1+1,J1,VEL%KDP)
          VEL%PQ1(K) = GRD%MSK(I1,J1,VEL%KDP)                                      &
                      * MAX(GRD%MSK(I1,J1,VEL%KDP),GRD%MSK(I1+1,J1,VEL%KDP))             &
                       * (1.-P1) * (1.-Q1)                                   &
                      /( DIV_X * DIV_Y + 1.E-16 )
          VEL%PQ2(K) = GRD%MSK(I1+1,J1,VEL%KDP)                                    &
                      * MAX(GRD%MSK(I1,J1,VEL%KDP),GRD%MSK(I1+1,J1,VEL%KDP))             &
                      *     P1  * (1.-Q1)                                    &
                      /( DIV_X * DIV_Y + 1.E-16 )
         DIV_X =  (1.-P1) * GRD%MSK(I1  ,J1+1,VEL%KDP) + P1 * GRD%MSK(I1+1,J1+1,VEL%KDP)
          VEL%PQ3(K) = GRD%MSK(I1,J1+1,VEL%KDP)                                    &
                      * MAX(GRD%MSK(I1,J1+1,VEL%KDP),GRD%MSK(I1+1,J1+1,VEL%KDP))         &
                      * (1.-P1) *     Q1                                     &
                      /( DIV_X * DIV_Y + 1.E-16 )
          VEL%PQ4(K) = GRD%MSK(I1+1,J1+1,VEL%KDP)                                  &
                      * MAX(GRD%MSK(I1,J1+1,VEL%KDP),GRD%MSK(I1+1,J1+1,VEL%KDP))         &
                      *     P1  *     Q1                                     &
                      /( DIV_X * DIV_Y + 1.E-16 )
        PRINT*,VEL%PQ1(K),VEL%PQ2(K),VEL%PQ3(K),VEL%PQ4(K)

        ENDIF
       ENDDO

!     VEL%FLC(2:VEL%NO) = 0

! ---
! COUNT GOOD OBSERVATIONS
    VEL%NC = 0
  DO K=1,VEL%NO
   IF(VEL%FLC(K).EQ.1)THEN
    VEL%NC = VEL%NC + 1
   ENDIF
  ENDDO

  PRINT*,'____________________________________________________________________'
  PRINT*,'NUMBER OF USEFULL:',VEL%NC
  PRINT*,VEL%RES(:)
  PRINT*,VEL%FLC(:)
  PRINT*,'____________________________________________________________________'


 ENDIF




END SUBROUTINE INT_PAR_VEL
