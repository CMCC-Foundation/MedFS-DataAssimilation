SUBROUTINE GET_OBS_XBT

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
  REAL(R8)      ::  P1, Q1, P, Q, R, SUMT, SUMI, TIMP, CNTI


   XBT%NO = 0
   XBT%NC = 0

  IF(OBS%XBT.EQ.0)RETURN

  OPEN(511,FILE='XBT_MIS.DAT',FORM='UNFORMATTED',STATUS='OLD')

! ---
! ALLOCATE MEMORY FOR OBSERVATIONS

  READ(511) XBT%NO

   PRINT*,'NO XBT: ',XBT%NO

   IF(XBT%NO.EQ.0)RETURN

   ALLOCATE ( XBT%INO(XBT%NO), XBT%FLG(XBT%NO), XBT%FLC(XBT%NO), XBT%PAR(XBT%NO))
   ALLOCATE ( XBT%LON(XBT%NO), XBT%LAT(XBT%NO), XBT%DPT(XBT%NO), XBT%TIM(XBT%NO))
   ALLOCATE ( XBT%VAL(XBT%NO), XBT%BAC(XBT%NO), XBT%INC(XBT%NO))
   ALLOCATE ( XBT%BIA(XBT%NO), XBT%ERR(XBT%NO))
   ALLOCATE ( XBT%RES(XBT%NO), XBT%B_A(XBT%NO))
   ALLOCATE ( XBT%IB(XBT%NO), XBT%JB(XBT%NO), XBT%KB(XBT%NO))
   ALLOCATE ( XBT%PB(XBT%NO), XBT%QB(XBT%NO), XBT%RB(XBT%NO))
   ALLOCATE ( XBT%PQ1(XBT%NO), XBT%PQ2(XBT%NO), XBT%PQ3(XBT%NO), XBT%PQ4(XBT%NO))
   ALLOCATE ( XBT%PQ5(XBT%NO), XBT%PQ6(XBT%NO), XBT%PQ7(XBT%NO), XBT%PQ8(XBT%NO))

! ---
! INITIALISE QUALITY FLAG
   XBT%FLG(:) = 1
   XBT%FLC(:) = 1

       READ (511)                                               &
        XBT%INO(1:XBT%NO), XBT%FLG(1:XBT%NO), XBT%PAR(1:XBT%NO) &
       ,XBT%LON(1:XBT%NO), XBT%LAT(1:XBT%NO)                    &
       ,XBT%DPT(1:XBT%NO), XBT%TIM(1:XBT%NO)                    &
       ,XBT%VAL(1:XBT%NO), XBT%BAC(1:XBT%NO)                    &
       ,XBT%ERR(1:XBT%NO), XBT%RES(1:XBT%NO)                    &
       ,XBT%IB(1:XBT%NO), XBT%JB(1:XBT%NO), XBT%KB(1:XBT%NO)    &
       ,XBT%PB(1:XBT%NO), XBT%QB(1:XBT%NO), XBT%RB(1:XBT%NO)
      CLOSE (511)

!	WRITE(*,*)'TEST SULLE FLAGS ',XBT%FLG(1:XBT%NO)
!	WRITE(*,*)'TEST SULLE LAT ',XBT%LAT(1:XBT%NO)
!	WRITE(*,*)'TEST SUI VAL ',XBT%VAL(1:XBT%NO)
!	WRITE(*,*)'TEST SULLE IB',XBT%IB(1:XBT%NO)

!!!    DO K = 1,XBT%NO
!!!     IF(XBT%LON(K).NE. 18.39400 .AND. XBT%LAT(K).NE. 39.73530 ) THEN
!!!        XBT%FLG(K) = 0
!!!     ELSE
!        XBT%LON(K) = 25.
!        XBT%LAT(K) = 37.35
!!!        XBT%RES(K) = 1.
!!!        XBT%ERR(K) = 0.01
!!!     ENDIF
!!!    ENDDO

! ---
! VERTICAL INTERPOLATION PARAMETERS
    DO K = 1,XBT%NO
     IF(XBT%FLG(K).EQ.1)THEN
       XBT%KB(K) = GRD%KM-1
!     WRITE(*,*) 'PRE VERT INT: ',K,XBT%KB(K),XBT%RB(K)
     DO KK = 1,GRD%KM-1
!     WRITE(*,*) 'VERT. INT: ',K,KK,XBT%DPT(K),GRD%DEP(KK),GRD%DEP(KK+1)
      IF( XBT%DPT(K).GE.GRD%DEP(KK) .AND. XBT%DPT(K).LT.GRD%DEP(KK+1) ) THEN
       XBT%KB(K) = KK
       XBT%RB(K) = (XBT%DPT(K) - GRD%DEP(KK)) / (GRD%DEP(KK+1) - GRD%DEP(KK))
      ENDIF
     ENDDO
!     WRITE(*,*) 'POST VERT INT: ',K,XBT%KB(K),XBT%RB(K)
     ENDIF
    ENDDO

! ---
! ELIMINIAMO IL THINNING DELLE OSSERVAZIONI, IN ATTESA DI VEDERE IN CHE MODO SERVIRSENE
!! THIN OBSERVATIONS
!    DO K = 1,XBT%NO-1
!     IF(XBT%FLG(K).EQ.1)THEN
!     PRINT*,'THIN CONTROL: XBT%KB(K), XBT%FLG(K),K ',XBT%KB(K),XBT%FLG(K),K
!         KK = K + 1
!         CNTI = 1.
!       DO WHILE(KK.LE.XBT%NO .AND. XBT%KB(K).EQ.XBT%KB(MIN(KK,XBT%NO)) .AND. XBT%FLG(MIN(KK,XBT%NO)).EQ.1)
!         XBT%VAL(K) = XBT%VAL(K) + XBT%VAL(KK)
!         XBT%BAC(K) = XBT%BAC(K) + XBT%BAC(KK)
!         XBT%RES(K) = XBT%RES(K) + XBT%RES(KK)
!         XBT%DPT(K) = XBT%DPT(K) + XBT%DPT(KK)
!         XBT%FLG(KK) = 0
!         CNTI = CNTI + 1.
!         KK = KK + 1
!       ENDDO
!         XBT%VAL(K) = XBT%VAL(K)/CNTI
!         XBT%BAC(K) = XBT%BAC(K)/CNTI
!         XBT%RES(K) = XBT%RES(K)/CNTI
!         XBT%DPT(K) = XBT%DPT(K)/CNTI
!         KK = XBT%KB(K)
!         XBT%RB(K) = (XBT%DPT(K) - GRD%DEP(KK)) / (GRD%DEP(KK+1) - GRD%DEP(KK))
!     ENDIF
!    ENDDO


! RESIDUAL CHECK
!!  DO K=1,XBT%NO
!!   IF(XBT%PAR(K).EQ.1 .AND. ABS(XBT%RES(K)).GT.5.0) THEN
!!	XBT%FLG(K) = 0
!!	WRITE(*,*)' IL RESIDUAL CHECK HA ELIMINATO IL DATO',K,XBT%RES(K)
!!   END IF
!!  ENDDO


! ---
! COUNT GOOD OBSERVATIONS
    XBT%NC = 0
  DO K=1,XBT%NO
   IF(XBT%FLG(K).EQ.1)THEN
    XBT%NC = XBT%NC + 1
   ELSE
    XBT%FLC(K) = 0
    XBT%BIA(K) = 0.
    XBT%RES(K) = 0.
    XBT%INC(K) = 0.
    XBT%B_A(K) = 0.
    XBT%PQ1(K) = 0.
    XBT%PQ2(K) = 0.
    XBT%PQ3(K) = 0.
    XBT%PQ4(K) = 0.
    XBT%PQ5(K) = 0.
    XBT%PQ6(K) = 0.
    XBT%PQ7(K) = 0.
    XBT%PQ8(K) = 0.
   ENDIF
  ENDDO


   XBT%FLC(:) = XBT%FLG(:)

    PRINT*,' HANNO PASSATO IL THINNING CHECK', XBT%NC, ' OSSERVAZIONI'
!    PRINT*,' STAMPA DELLE OSSERVAZIONI'
!    DO K=1,XBT%NO
!      IF(XBT%FLG(K).EQ.1)THEN
!    PRINT*,'K: ',K,XBT%FLC(K),XBT%BIA(K),XBT%RES(K),XBT%INC(K),XBT%B_A(K),XBT%PQ1(K), &
!    XBT%PQ2(K),XBT%PQ3(K),XBT%PQ4(K),XBT%PQ5(K),XBT%PQ6(K),XBT%PQ7(K),XBT%PQ8(K)
!      END IF
!    END DO

END SUBROUTINE GET_OBS_XBT



SUBROUTINE INT_PAR_XBT

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

  INTEGER(I4)   ::  K, K1,II, JJ,IM,JM,PIER1,PIER2,PIER3,PIER4
  INTEGER(I4)   ::  I1, J1, KK, I, J, ITER, IDEP,PRECI1,PRECJ1
  REAL(R8)      ::  P1, Q1, R1, P, Q, R, SUMT, SUMI, TIMP, CNTI
  REAL(R8)      ::  MSK4, DIV_X, DIV_Y,PRECLAT,PRECLON,PRECP1,PRECQ1


 IF(XBT%NO.GT.0) THEN

    XBT%FLC(:) = XBT%FLG(:)
	PIER1=0
	PIER2=0
	DO K=1,XBT%NO
	 IF (XBT%FLG(K).EQ.1) PIER1=PIER1+1
	 IF (XBT%FLG(K).EQ.0) PIER2=PIER2+1
	 IF (XBT%FLC(K).EQ.1) PIER3=PIER3+1
	 IF (XBT%FLC(K).EQ.0) PIER4=PIER4+1
	END DO
	WRITE(*,*) 'ABBIAMO IN INGRESSO ', XBT%NO, ' PROFILI'
	WRITE(*,*) 'RISULTANO ', PIER1, 'ACCETTATI IN XBT%FLG'
	WRITE(*,*) 'RISULTANO ', PIER2, 'RIFIUTATI IN XBT%FLG'
	WRITE(*,*) 'RISULTANO ', PIER3, 'ACCETTATI IN XBT%FLC'
	WRITE(*,*) 'RISULTANO ', PIER4, 'RIFIUTATI IN XBT%FLC'
	

! ---
! HORIZONTAL INTERPOLATION PARAMETERS
! IN QUESTO PUNTO VA MODIFICATA LA PROCEDURA PER FARE L'INTERPOLAZIONE SU ORCA2
! PER ORA NE USO UNA PIUTTOSTO INEFFICIENTE CHE LAVORA SULL'INTERO DOMINIO
! VA AGGIUNTO IL CASO DELLE LONGITUDINI CHE CADONO NEL CAMBIO DI DATA
! LE VARIABILI PRECLON E PRECLAT IMMAGAZZINANO L'E ULTIME LAT,LON ACCETTATE COME
! VALIDE, IN MODO DA NON DOVER FARE LA CERCA PER TUTTI I DATI APPARTENENTI AL PROFILO
! LE VARIABILI PRECI1,PRECJ1,PRECP1,PRECQ1 SERVONO AD ASSEGNARE DIRETTAMENTE I VALORI
! CALCOLATI PER IL PRIMO PUNTO DEL PROFILO, SENZA PASSARE PER L'ALGORITMO DI RICERCA

       WRITE(*,*) 'PIERLUIGI: CHECK DI IMPORTAZIONE XBT IN 3DVAR'

!	AGGIUNGO TEMPORANEAMENTE LE DIMENSIONI DELLA GRIGLIA ORCA2 NELLE VARIABILI IM E JM
!	QUESTO ALGORITMO E' PIUTTOSTO INEFFICIENTE E NON RECUPERA I DATI TRA 178E E 180W
! BISOGNA INSERIRE UN ALGORITMO MIGLIORE DI RICERCA BINARIA, UTILIZZANDO UNA GRIGLIA COARSE
	IM=182
	JM=149
        PRECLON=999.
        PRECLAT=999.
        PRECI1=0
        PRECJ1=0
        PRECP1=0
        PRECQ1=0
    DO K = 1,XBT%NO
	IF (XBT%FLC(K).EQ.1) THEN
        IF ((XBT%LON(K).NE.PRECLON).AND.(XBT%LAT(K).NE.PRECLAT)) THEN
	XBT%FLC(K) = 0
        DO II=1,IM-1
          DO JJ=1,JM-1
            IF ((XBT%LON(K).GE.GRD%LON(II,JJ).AND.XBT%LON(K).LT.GRD%LON(II+1,JJ))  &
        .AND.(XBT%LAT(K).GE.GRD%LAT(II,JJ).AND.XBT%LAT(K).LT.GRD%LAT(II,JJ+1))) THEN
!       WRITE(*,*) '3DVAR LAT: ',XBT%LAT(K),GRD%LAT(II,JJ),GRD%LAT(II,JJ+1),II,JJ,K
!       WRITE(*,*) '3DVAR LON: ',XBT%LON(K),GRD%LON(II,JJ),GRD%LON(II+1,JJ),II,JJ,K
	     XBT%FLC(K) = 1
             J1=JJ
             I1=II
             P1=(XBT%LAT(K)-GRD%LAT(II,JJ))/(GRD%LAT(II,JJ+1)-GRD%LAT(II,JJ))
             Q1=(XBT%LON(K)-GRD%LON(II,JJ))/(GRD%LON(II+1,JJ)-GRD%LON(II,JJ))
	     XBT%IB(K) = I1
             XBT%JB(K) = J1
             XBT%PB(K) = P1
             XBT%QB(K) = Q1
            END IF
          END DO
        END DO
        PRECLON=XBT%LON(K)
        PRECLAT=XBT%LAT(K)
        PRECI1=I1
        PRECJ1=J1
        PRECP1=P1
        PRECQ1=Q1
	ELSE
            XBT%IB(K) = PRECI1
            XBT%JB(K) = PRECJ1
            XBT%PB(K) = PRECP1
            XBT%QB(K) = PRECQ1
	END IF
	END IF
     END DO

! ---
! UNDEFINE MASKED FOR MULTIGRID
    DO K = 1,XBT%NO
     IF(XBT%FLC(K).EQ.1)THEN
      I1 = XBT%IB(K)
      J1 = XBT%JB(K)
      IDEP = XBT%KB(K)+1
       MSK4 = GRD%MSK(I1,J1,IDEP) + GRD%MSK(I1+1,J1,IDEP) + GRD%MSK(I1,J1+1,IDEP) + GRD%MSK(I1+1,J1+1,IDEP)
!	WRITE(*,*)'COME STIAMO TRASFERENDO I1 E J1? ', I1,J1,IDEP,XBT%IB(K),XBT%JB(K),XBT%KB(K),K
!	WRITE(*,*) 'GRD%MSK(I1,J1,IDEP),GRD%MSK(I1+1,J1,IDEP),GRD%MSK(I1,J1+1,IDEP),GRD%MSK(I1+1,J1+1,IDEP)'
!	WRITE(*,*) GRD%MSK(I1,J1,IDEP),GRD%MSK(I1+1,J1,IDEP),GRD%MSK(I1,J1+1,IDEP),GRD%MSK(I1+1,J1+1,IDEP)
      IF(MSK4.LT.1.0)XBT%FLC(K) = 0
     ENDIF
    ENDDO

! ---
! HORIZONTAL INTERPOLATION PARAMETERS FOR EACH MASKED GRID
       DO K = 1,XBT%NO
        IF(XBT%FLC(K) .EQ. 1) THEN

         I1=XBT%IB(K)
         P1=XBT%PB(K)
         J1=XBT%JB(K)
         Q1=XBT%QB(K)
         K1=XBT%KB(K)

!	WRITE(*,*) 'I1,P1,J1,Q1,K1,K'
!	WRITE(*,*) I1,P1,J1,Q1,K1,K
!	WRITE(*,*) 'GRD%MSK(I1,J1  ,K1),GRD%MSK(I1+1,J1  ,K1)'
!	WRITE(*,*) GRD%MSK(I1,J1  ,K1),GRD%MSK(I1+1,J1  ,K1)
!	WRITE(*,*)'GRD%MSK(I1,J1+1,K1),GRD%MSK(I1+1,J1+1,K1)'
!	WRITE(*,*)GRD%MSK(I1,J1+1,K1),GRD%MSK(I1+1,J1+1,K1)


         DIV_Y =  (1.-Q1) * MAX(GRD%MSK(I1,J1  ,K1),GRD%MSK(I1+1,J1  ,K1))     &
                 +    Q1  * MAX(GRD%MSK(I1,J1+1,K1),GRD%MSK(I1+1,J1+1,K1))
         DIV_X =  (1.-P1) * GRD%MSK(I1  ,J1,K1) + P1 * GRD%MSK(I1+1,J1,K1)

!	 WRITE(*,*)'FIRST: DIV_Y,DIV_X ',DIV_Y,DIV_X

          XBT%PQ1(K) = GRD%MSK(I1,J1,K1)                                       &
                      * MAX(GRD%MSK(I1,J1,K1),GRD%MSK(I1+1,J1,K1))             &
                       * (1.-P1) * (1.-Q1)                                     &
                      /( DIV_X * DIV_Y + 1.E-16 )
          XBT%PQ2(K) = GRD%MSK(I1+1,J1,K1)                                     &
                      * MAX(GRD%MSK(I1,J1,K1),GRD%MSK(I1+1,J1,K1))             &
                      *     P1  * (1.-Q1)                                      &
                      /( DIV_X * DIV_Y + 1.E-16 )

! 	 WRITE(*,*)'XBT%PQ1(K),XBT%PQ2(K)',XBT%PQ1(K),XBT%PQ2(K)

         DIV_X =  (1.-P1) * GRD%MSK(I1  ,J1+1,K1) + P1 * GRD%MSK(I1+1,J1+1,K1)

!	 WRITE(*,*)'SECND: DIV_Y,DIV_X ',DIV_Y,DIV_X

          XBT%PQ3(K) = GRD%MSK(I1,J1+1,K1)                                     &
                      * MAX(GRD%MSK(I1,J1+1,K1),GRD%MSK(I1+1,J1+1,K1))         &
                      * (1.-P1) *     Q1                                       &
                      /( DIV_X * DIV_Y + 1.E-16 )
          XBT%PQ4(K) = GRD%MSK(I1+1,J1+1,K1)                                   &
                      * MAX(GRD%MSK(I1,J1+1,K1),GRD%MSK(I1+1,J1+1,K1))         &
                      *     P1  *     Q1                                       &
                      /( DIV_X * DIV_Y + 1.E-16 )

!	 WRITE(*,*)'XBT%PQ3(K),XBT%PQ4(K)',XBT%PQ3(K),XBT%PQ4(K)

         K1=XBT%KB(K) + 1
!	 WRITE(*,*)'K1=XBT%KB(K) + 1 ', K1
	
         DIV_Y =  (1.-Q1) * MAX(GRD%MSK(I1,J1  ,K1),GRD%MSK(I1+1,J1  ,K1))     &
                 +    Q1  * MAX(GRD%MSK(I1,J1+1,K1),GRD%MSK(I1+1,J1+1,K1))
         DIV_X =  (1.-P1) * GRD%MSK(I1  ,J1,K1) + P1 * GRD%MSK(I1+1,J1,K1)

!	  WRITE(*,*)' THIRD: DIV_Y,DIV_X ',DIV_Y,DIV_X

          XBT%PQ5(K) = GRD%MSK(I1,J1,K1)                                       &
                      * MAX(GRD%MSK(I1,J1,K1),GRD%MSK(I1+1,J1,K1))             &
                       * (1.-P1) * (1.-Q1)                                     &
                      /( DIV_X * DIV_Y + 1.E-16 )
          XBT%PQ6(K) = GRD%MSK(I1+1,J1,K1)                                     &
                      * MAX(GRD%MSK(I1,J1,K1),GRD%MSK(I1+1,J1,K1))             &
                      *     P1  * (1.-Q1)                                      &
                      /( DIV_X * DIV_Y + 1.E-16 )

!	 WRITE(*,*)'XBT%PQ5(K),XBT%PQ6(K) ',XBT%PQ5(K),XBT%PQ6(K)

         DIV_X =  (1.-P1) * GRD%MSK(I1  ,J1+1,K1) + P1 * GRD%MSK(I1+1,J1+1,K1)

!	 WRITE(*,*)' FOURTH: DIV_Y,DIV_X ',DIV_Y,DIV_X

          XBT%PQ7(K) = GRD%MSK(I1,J1+1,K1)                                     &
                      * MAX(GRD%MSK(I1,J1+1,K1),GRD%MSK(I1+1,J1+1,K1))         &
                      * (1.-P1) *     Q1                                       &
                      /( DIV_X * DIV_Y + 1.E-16 )
          XBT%PQ8(K) = GRD%MSK(I1+1,J1+1,K1)                                   &
                      * MAX(GRD%MSK(I1,J1+1,K1),GRD%MSK(I1+1,J1+1,K1))         &
                      *     P1  *     Q1                                       &
                      /( DIV_X * DIV_Y + 1.E-16 )

!	 WRITE(*,*)'XBT%PQ7(K),XBT%PQ8(K) ',XBT%PQ7(K),XBT%PQ8(K)

         R1=XBT%RB(K)
          XBT%PQ1(K) = (1.-R1) * XBT%PQ1(K)
          XBT%PQ2(K) = (1.-R1) * XBT%PQ2(K)
          XBT%PQ3(K) = (1.-R1) * XBT%PQ3(K)
          XBT%PQ4(K) = (1.-R1) * XBT%PQ4(K)
          XBT%PQ5(K) =     R1  * XBT%PQ5(K)
          XBT%PQ6(K) =     R1  * XBT%PQ6(K)
          XBT%PQ7(K) =     R1  * XBT%PQ7(K)
          XBT%PQ8(K) =     R1  * XBT%PQ8(K)
	
!	 WRITE(*,*)'R1=XBT%RB(K) ',R1
!	 WRITE(*,*)'XBT%PQ1(K), ... , XBT%PQ8(K) '
!	 WRITE(*,*) XBT%PQ1(K),XBT%PQ2(K),XBT%PQ3(K),XBT%PQ4(K)
!	 WRITE(*,*) XBT%PQ5(K),XBT%PQ6(K),XBT%PQ7(K),XBT%PQ8(K)

        ENDIF
       ENDDO


! ---
! COUNT GOOD OBSERVATIONS
    XBT%NC = 0
  DO K=1,XBT%NO
   IF(XBT%FLC(K).EQ.1)THEN
    XBT%NC = XBT%NC + 1
   ENDIF
  ENDDO


 ENDIF


END SUBROUTINE INT_PAR_XBT
