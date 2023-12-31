SUBROUTINE REARRALPHA_LV

USE RECFILTER

!-- RECOMPUTE ALPHAS ACCORDING TO EXTENDED GRID
!   THIS ROUTINE IS SUITABLE TO BE RECALLED WITHIN RF ROUTINES

USE GRD_STR
USE MYNETCDF
USE MYFRTPROF
USE MPIREL

IMPLICIT  NONE

INTEGER(KIND=I4) :: NX,NY,NZ
INTEGER(KIND=I4) :: I,J,K,KK,K2

CALL MYFRTPROF_WALL('REARRALPHA_LV: REARRANGE ALPHA',0)

WRITE(IOUNLOG,*)
WRITE(IOUNLOG,*) ' /// REARRANGE ALPHA'
WRITE(IOUNLOG,*)

#ifndef USE_POINTERS
IF( .NOT.ALLOCATED(GRD%DX).OR..NOT.ALLOCATED(GRD%DY).OR.&
  & .NOT.ALLOCATED(GRD%MSR)) THEN
   WRITE(IOUNERR,*) 'REARRALPHA CALLED BEFORE GRID INITIALIZATION'
   CALL ABOR1('ERROR IN REARRALPHA, CANNOT PROCEED')
ENDIF
#endif

IF(.NOT.LLRFINIT) &
& CALL ABOR1('REARRALPHA CALLED BEFORE RF INITIALIZATION')

NX=GRD%IM
NY=GRD%JM
NZ=GRD%KM

WRITE(IOUNLOG,*) ' ABOUT TO ALLOCATE FCT:',NX,NY,NTOTV

ALLOCATE( FCT(NX,NY,NTOTV) )
FCT=0._R8

WRITE(IOUNLOG,*) ' ABOUT TO ALLOCATE RF_X ARRAYS:',NY,IMAX,NTOTV
WRITE(IOUNLOG,*) '               AND RF_Y ARRAYS:',NX,JMAX,NTOTV

ALLOCATE( RF_AEX(NY,IMAX,NTOTV),&
        & RF_BEX(NY,IMAX,NTOTV),&
        & RF_AEY(NX,JMAX,NTOTV),&
        & RF_BEY(NX,JMAX,NTOTV) )

RF_AEX = 0._R8
RF_AEY = 0._R8
RF_BEX = 0._R8
RF_BEY = 0._R8

IF( RCF_TYPE .EQ. 2 ) THEN
  ALLOCATE( RF_GEX(NY,IMAX,NTOTV),&
        & RF_DEX(NY,IMAX,NTOTV),&
        & RF_GEY(NX,JMAX,NTOTV),&
        & RF_DEY(NX,JMAX,NTOTV) )

  RF_GEX = 0._R8
  RF_GEY = 0._R8
  RF_DEX = 0._R8
  RF_DEY = 0._R8
ENDIF


#ifdef SHARED_MEMORY
!$OMP PARALLEL DEFAULT(SHARED),PRIVATE(K,K2,J,KK,I)
!$OMP DO SCHEDULE(DYNAMIC,1)
#endif
DO K = 1, NTOTV
   K2=K

   DO J = 1, NY
      KK = ISTPNS(1,J,K)
      IF(LLMSR(1,J,K2)) THEN
          KK = KK + 1
          RF_AEX(J,1:KK,K) = ALXNS(1,J,K)
          RF_BEX(J,1:KK,K) = BTXNS(1,J,K)
          IF( RCF_TYPE .EQ. 2 ) THEN
            RF_GEX(J,1:KK,K) = GMXNS(1,J,K)
            RF_DEX(J,1:KK,K) = DLXNS(1,J,K)
          ENDIF
      ENDIF
      DO I = 2,NX
         IF(.NOT.LLMSR(I,J,K2) .AND. LLMSR(I-1,J,K2)) THEN
                 RF_AEX(J,KK+1:KK+ISTPNS(I,J,K),K) = ALXNS(I,J,K)
                 RF_BEX(J,KK+1:KK+ISTPNS(I,J,K),K) = BTXNS(I,J,K)
                 IF( RCF_TYPE .EQ. 2 ) THEN
                   RF_GEX(J,KK+1:KK+ISTPNS(I,J,K),K) = GMXNS(I,J,K)
                   RF_DEX(J,KK+1:KK+ISTPNS(I,J,K),K) = DLXNS(I,J,K)
                 ENDIF
                 KK = KK + ISTPNS(I,J,K)
         ELSEIF(LLMSR(I,J,K2) .AND. .NOT.LLMSR(I-1,J,K2)) THEN
                 RF_AEX(J,KK+1:KK+ISTPNS(I,J,K)+1,K) = ALXNS(I,J,K)
                 RF_BEX(J,KK+1:KK+ISTPNS(I,J,K)+1,K) = BTXNS(I,J,K)
                 IF( RCF_TYPE .EQ. 2 ) THEN
                   RF_GEX(J,KK+1:KK+ISTPNS(I,J,K)+1,K) = GMXNS(I,J,K)
                   RF_DEX(J,KK+1:KK+ISTPNS(I,J,K)+1,K) = DLXNS(I,J,K)
                 ENDIF
                 KK = KK + ISTPNS(I,J,K) + 1
         ELSEIF(LLMSR(I,J,K2)) THEN
                 RF_AEX(J,KK+1,K) = ALXNS(I,J,K)
                 RF_BEX(J,KK+1,K) = BTXNS(I,J,K)
                 IF( RCF_TYPE .EQ. 2 ) THEN
                   RF_GEX(J,KK+1,K) = GMXNS(I,J,K)
                   RF_DEX(J,KK+1,K) = DLXNS(I,J,K)
                 ENDIF
                 KK = KK + 1
         ENDIF
      ENDDO
   ENDDO
   IF (LL_PRINT_RECFILT) THEN
      WRITE(IOUNLOG,*) 'RF_AEX ', K, ': ',KK,SUM(RF_AEX(:,1:KK,K))/(NY*KK),&
      & MINVAL(RF_AEX(:,1:KK,K)),MAXVAL(RF_AEX(:,1:KK,K))
      WRITE(IOUNLOG,*) 'RF_BEX ', K, ': ',KK,SUM(RF_BEX(:,1:KK,K))/(NY*KK),&
      & MINVAL(RF_BEX(:,1:KK,K)),MAXVAL(RF_BEX(:,1:KK,K))
   ENDIF
   CALL FLUSH(IOUNLOG)

   DO I = 1, NX
      KK = JSTPNS(I,1,K)
      IF(LLMSR(I,1,K2)) THEN
          KK = KK + 1
          RF_AEY(I,1:KK,K) = ALYNS(I,1,K)
          RF_BEY(I,1:KK,K) = BTYNS(I,1,K)
          IF( RCF_TYPE .EQ. 2 ) THEN
            RF_GEY(I,1:KK,K) = GMYNS(I,1,K)
            RF_DEY(I,1:KK,K) = DLYNS(I,1,K)
          ENDIF
      ENDIF
      DO J = 2, NY
         IF(.NOT.LLMSR(I,J,K2) .AND. LLMSR(I,J-1,K2)) THEN
                 RF_AEY(I,KK+1:KK+JSTPNS(I,J,K),K) = ALYNS(I,J,K)
                 RF_BEY(I,KK+1:KK+JSTPNS(I,J,K),K) = BTYNS(I,J,K)
                 IF( RCF_TYPE .EQ. 2 ) THEN
                   RF_GEY(I,KK+1:KK+JSTPNS(I,J,K),K) = GMYNS(I,J,K)
                   RF_DEY(I,KK+1:KK+JSTPNS(I,J,K),K) = DLYNS(I,J,K)
                 ENDIF
                 KK = KK + JSTPNS(I,J,K)
         ELSEIF(LLMSR(I,J,K2) .AND. .NOT.LLMSR(I,J-1,K2)) THEN
                 RF_AEY(I,KK+1:KK+JSTPNS(I,J,K)+1,K) = ALYNS(I,J,K)
                 RF_BEY(I,KK+1:KK+JSTPNS(I,J,K)+1,K) = BTYNS(I,J,K)
                 IF( RCF_TYPE .EQ. 2 ) THEN
                   RF_GEY(I,KK+1:KK+JSTPNS(I,J,K)+1,K) = GMYNS(I,J,K)
                   RF_DEY(I,KK+1:KK+JSTPNS(I,J,K)+1,K) = DLYNS(I,J,K)
                 ENDIF
                 KK = KK + JSTPNS(I,J,K) + 1
         ELSEIF(LLMSR(I,J,K2)) THEN
                 RF_AEY(I,KK+1,K) = ALYNS(I,J,K)
                 RF_BEY(I,KK+1,K) = BTYNS(I,J,K)
                 IF( RCF_TYPE .EQ. 2 ) THEN
                   RF_GEY(I,KK+1,K) = GMYNS(I,J,K)
                   RF_DEY(I,KK+1,K) = DLYNS(I,J,K)
                 ENDIF
                 KK = KK + 1
         ENDIF
      ENDDO
   ENDDO
   IF (LL_PRINT_RECFILT) THEN
      WRITE(IOUNLOG,*) 'RF_AEY ', K, ': ',KK,SUM(RF_AEY(:,1:KK,K))/(NX*KK),&
      & MINVAL(RF_AEY(:,1:KK,K)),MAXVAL(RF_AEY(:,1:KK,K))
      WRITE(IOUNLOG,*) 'RF_BEY ', K, ': ',KK,SUM(RF_BEY(:,1:KK,K))/(NX*KK),&
      & MINVAL(RF_BEY(:,1:KK,K)),MAXVAL(RF_BEY(:,1:KK,K))
   ENDIF
   CALL FLUSH(IOUNLOG)

   WHERE( LLMSR(:,:,K2) ) FCT(:,:,K) = 1._R8

ENDDO
#ifdef SHARED_MEMORY
!$OMP END DO
!$OMP END PARALLEL
#endif

IF(LLPRINTRCFLC) THEN
  CALL WRITE_VARNCDF('RF_AEX.NC'//TRIM(CMPIDOM),GRD%JM,IMAX,NTOTV,&
  & RF_AEX,'rf_aex',' ','RF_AEX RF coeff')
  CALL WRITE_VARNCDF('RF_BEX.NC'//TRIM(CMPIDOM),GRD%JM,IMAX,NTOTV,&
  & RF_BEX,'rf_bex',' ','RF_BEX RF coeff')
  CALL WRITE_VARNCDF('RF_AEY.NC'//TRIM(CMPIDOM),GRD%IM,JMAX,NTOTV,&
  & RF_AEY,'rf_aey',' ','RF_AEY RF coeff')
  CALL WRITE_VARNCDF('RF_BEY.NC'//TRIM(CMPIDOM),GRD%IM,JMAX,NTOTV,&
  & RF_BEY,'rf_bey',' ','RF_BEY RF coeff')
  IF( RCF_TYPE .EQ. 2 ) THEN
        CALL WRITE_VARNCDF('RF_GEX.NC'//TRIM(CMPIDOM),GRD%JM,IMAX,NTOTV,&
        & RF_GEX,'rf_gex',' ','RF coeff')
        CALL WRITE_VARNCDF('RF_DEX.NC'//TRIM(CMPIDOM),GRD%JM,IMAX,NTOTV,&
        & RF_DEX,'rf_dex',' ','RF coeff')
        CALL WRITE_VARNCDF('RF_GEY.NC'//TRIM(CMPIDOM),GRD%IM,JMAX,NTOTV,&
        & RF_GEY,'rf_gey',' ','RF coeff')
        CALL WRITE_VARNCDF('RF_DEY.NC'//TRIM(CMPIDOM),GRD%IM,JMAX,NTOTV,&
        & RF_DEY,'rf_dey',' ','RF coeff')   
  ENDIF
  CALL WRITE_VARNCDF('RF_INX.NC'//TRIM(CMPIDOM),GRD%IM,GRD%JM,NTOTV,&
  & INX,'inx',' ','INX RF coeff')
  CALL WRITE_VARNCDF('RF_JNX.NC'//TRIM(CMPIDOM),GRD%IM,GRD%JM,NTOTV,&
  & JNX,'jnx',' ','JNX RF coeff')
  CALL WRITE_VARNCDF('RF_IMX.NC'//TRIM(CMPIDOM),NTOTV,&
  & IMX,'imx',' ','IMX RF coeff')
  CALL WRITE_VARNCDF('RF_JMX.NC'//TRIM(CMPIDOM),NTOTV,&
  & JMX,'jmx',' ','JMX RF coeff')
ENDIF

DEALLOCATE( ALXNS, BTXNS, ALYNS, BTYNS)
IF( RCF_TYPE .EQ. 2 ) DEALLOCATE( GMXNS, DLXNS, GMYNS, DLYNS)

CALL MYFRTPROF_WALL('REARRALPHA_LV: REARRANGE ALPHA',1)

END SUBROUTINE REARRALPHA_LV
