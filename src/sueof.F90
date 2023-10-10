SUBROUTINE SUEOF

!-----------------------------------------------------------------------
!                                                                      !
! DEFINE FILTER CONSTANTS, EOFS, ETC.                                  !
! 2010 SUMMER REVISION, NETCDF SUPPORT                                 !
!                                                                      !
! VERSION 1: S.DOBRICIC 2006                                           !
!-----------------------------------------------------------------------

  USE SET_KND
  USE GRD_STR
  USE ICE
  USE EOF_STR
  USE MYNETCDF
  USE IOUNITS
  USE MYFRTPROF,ONLY : MYFRTPROF_WALL
  USE RUN
  USE STATS, ONLY : MEAN

  IMPLICIT NONE

  INTEGER(I4)                 :: NRG, NEC, K, NSPL, I, J, KP, KP2, KP3
  INTEGER(I4)                 :: KSTP, NSLC,REGNUM, IERR
  INTEGER(I4) :: JPAR,JLEV,KK,KPARS,KRANK,KDIMS(3),KLEV

  REAL(R8) :: ZL(ROS%NEOF,3)
  REAL(R8) :: ZC(ROS%NEOF,ROS%KMT,3)
  REAL(R8), ALLOCATABLE :: MLD(:,:), BGFAC2(:,:), BGFAC1(:)
  INTEGER(I4), ALLOCATABLE :: CNT(:)
  REAL(R8) :: ZMLD, ZL1, ZL2, ZDL, RPI
  LOGICAL  :: LL_MLD_VLOC

CALL MYFRTPROF_WALL('SUEOF: SET-UP EOF',0)

! ---
! EOFS

     WRITE(IOUNLOG,*)
     WRITE(IOUNLOG,*) '  /// EOFS SETUP'

     IF(ROS%EOGNX .LE. 0) ROS%EOGNX=GRD%IM
     IF(ROS%EOGNY .LE. 0) ROS%EOGNY=GRD%JM
     IF(ROS%EOGNZ .LE. 0) ROS%EOGNZ=GRD%KM
     IF(ROS%EOGNZ .NE. GRD%KM) THEN
          WRITE(IOUNERR,*) 'ROS%EOGNZ != GRD%KM, ',ROS%EOGNZ,GRD%KM
          WRITE(IOUNERR,*) '   THIS IS NOT SUPPORTED FOR THE TIME BEING'
          CALL ABOR1('SUEOF: VERTICAL DIMENSION OF EOFS AND FIRST GUESS MISMATCH')
     ENDIF
     IF(ROS%KMT .NE. GRD%KM*PSV3D + PSV2D) THEN
          CALL ABOR1('SUEOF : KMT != STATE VECTOR')
     ENDIF

     IF(ROS%NREG .EQ. 0 ) ROS%NREG = GRD%IM*GRD%JM
     WRITE(IOUNLOG,*) ' EOF MODES    : ',ROS%NEOF
     WRITE(IOUNLOG,*) ' EOF REGIONS  : ',ROS%NREG
     WRITE(IOUNLOG,*) ' EOF LEVELS   : ',ROS%KMT
     CALL FLUSH(IOUNLOG)

     ALLOCATE ( GRD%RO( GRD%IM, GRD%JM, ROS%NEOF))
     ALLOCATE ( GRD%RO_AD( GRD%IM, GRD%JM, ROS%NEOF))
     ALLOCATE ( ROS%EVC( ROS%NREG, ROS%KMT, ROS%NEOF), ROS%EVA( ROS%NREG, ROS%NEOF) )

     WRITE(IOUNLOG,*) ' EOF ALLOCATED, CALLING READ_NCEOF'
     CALL FLUSH(IOUNLOG)

     CALL READ_NCEOF(TRIM(ROS%EOF_FILE),ROS%NREG,GRD%IM,GRD%JM,&
                   & ROS%KMT,ROS%NEOF,ROS%EVA,ROS%EVC,GRD%REG,LSQUAREDEVA)
   
     WRITE(IOUNLOG,*) ' READ_NCEOF CALLED'
     CALL FLUSH(IOUNLOG)

     IF( LL_ICEBGFILT ) CALL ICEBGFILT

     IF(REDNMC.GT.0.1_R8 .AND. REDNMC.LT.10._R8 ) THEN
       WRITE(IOUNLOG,*)
       WRITE(IOUNLOG,*) ' *** USING REDNMC == ',REDNMC
       WRITE(IOUNLOG,*)
       ROS%EVA=ROS%EVA*REDNMC
     ENDIF

     IF(LL_BGFACT) THEN 
       ALLOCATE( BGFAC2(GRD%IM,GRD%JM) )
       ALLOCATE( BGFAC1(ROS%NREG) )
       ALLOCATE( CNT  (ROS%NREG) )
       CNT = 0
       BGFAC1 = 0._R8
       WRITE(IOUNLOG,*) ' READING GRID-DEPENDENT INFLATION OF B'
       CALL GETNCVAR('RATIO_BGERR.nc','ratio',GRD%IM,GRD%JM,BGFAC2)
       DO J=1,GRD%JM
         DO I=1,GRD%IM
           IF( GRD%MSK(I,J,1) .GT. 0.5_R8 ) THEN
             KK=GRD%REG(I,J)
             BGFAC1(KK) = BGFAC1(KK) + BGFAC2(I,J)
             CNT   (KK) = CNT   (KK) + 1
           ENDIF
         ENDDO
       ENDDO
       WRITE(IOUNLOG,*) ' MIN/MEAN/MAX VAL OF COUNTS :', &
            & MINVAL(CNT), SUM(CNT)/REAL(GRD%IM*GRD%JM), MAXVAL(CNT)
       WRITE(IOUNLOG,*) ' MIN/MEAN/MAX VAL OF COUNTS WITHOUT LAND :', &
            & MINVAL(CNT,MASK=(CNT>0)),&
            & SUM(CNT,MASK=(CNT>0))/REAL(SUM(GRD%MSK(:,:,1))),&
            & MAXVAL(CNT,MASK=(CNT>0))
       DEALLOCATE( BGFAC2 )
       WHERE( CNT .GT. 0 ) BGFAC1 = BGFAC1 / REAL(CNT,R8)
       WHERE( CNT .EQ. 0 ) BGFAC1 = 1._R8
       DO KK=1,ROS%NEOF
           ROS%EVA(:,KK)=ROS%EVA(:,KK)*BGFAC1
       ENDDO
       DEALLOCATE( CNT, BGFAC1 )
     ENDIF

     DO KK=1,ROS%NEOF
      DO J=1,GRD%JM
       DO I=1,GRD%IM
          IF (GRD%REG(I,J) .GT. 0) THEN
             ROS%EVC(GRD%REG(I,J),1:GRD%KM,KK)=ROS%EVC(GRD%REG(I,J),1:GRD%KM,KK)*GRD%MSK(I,J,:)
             ROS%EVC(GRD%REG(I,J),(GRD%KM+1):(2*GRD%KM),KK)=ROS%EVC(GRD%REG(I,J),(GRD%KM+1):(2*GRD%KM),KK)*GRD%MSK(I,J,:)
             IF( LL_TQ2 .AND. GRD%DISTC(I,J) .LT. 100000._R8) ROS%EVC(GRD%REG(I,J),(2*GRD%KM+1):(2*GRD%KM+2),KK) = 0._R8
          ENDIF
       ENDDO
      ENDDO
     ENDDO

     IF( IDOUBLEDOM .EQ. 2 ) THEN

        WRITE(IOUNLOG,*) ' *** INITIALIZING PHYSICAL SPACE CONTROL VECTOR'
        WRITE(IOUNLOG,*)
        WRITE(IOUNLOG,*) ' READING FROM CR ANINCR INCRO, DIMS ARE : ',&
        & GRD%IM,GRD%JM,ROS%NEOF
        WRITE(IOUNLOG,*) ' SIZE OF GRD%RO IS', SIZE(GRD%RO)
        CALL GETNCVAR('ANINCR_HR.NC','INCRO',GRD%IM,GRD%JM,ROS%NEOF,GRD%RO)

        WRITE(IOUNLOG,*) ' *** INITIALIZATION END'
        CALL FLUSH(IOUNLOG)
     ENDIF


     LL_MLD_VLOC = ( NN_MLD_VLOC .NE. 0 )
     ALLOCATE( GRD%MVLOC(GRD%IM,GRD%JM,ROS%KMT) )
     GRD%MVLOC = 1._R8

     IF( LL_MLD_VLOC ) THEN
        ALLOCATE( MLD(GRD%IM,GRD%JM)  )
        CALL GETNCVAR('MLD.nc','somxl010',GRD%IM,GRD%JM,MLD )
        IF( NN_MLD_VLOC .EQ. 1 ) THEN
          DO J=1,GRD%JM 
            DO I=1,GRD%IM 
             IF( GRD%MSK(I,J,1) .GT. 0.5_R8 ) THEN
              DO K=1,GRD%KM
                IF( GRD%DEP(K) .LT. MLD(I,J) ) THEN
                    GRD%MVLOC(I,J,K) = 1._R8
                    KLEV=K
                ELSE
                    GRD%MVLOC(I,J,K) = 0._R8
                ENDIF
              ENDDO
             GRD%MVLOC(I,J,KLEV+1) = 0.5_R8
             ENDIF
            ENDDO
          ENDDO
        ELSE IF ( NN_MLD_VLOC .EQ. 2 ) THEN
         RPI = 2._R8*ASIN(1._R8)
         DO J=1,GRD%JM
          DO I=1,GRD%IM
           ZMLD = MLD(I,J)
           ZL1  = 10._R8
           IF ( ZMLD .LT. 15._R8 ) ZL1 = 5._R8
           ZL2  = ZMLD
           ZDL  = ZL2-ZL1
           DO K=1,GRD%KM
              IF( GRD%DEP(K) .LT. ZL1 ) THEN
                  GRD%MVLOC(I,J,K) = 1._R8
              ELSEIF ( GRD%DEP(K) .GT. ZL2 ) THEN
                  GRD%MVLOC(I,J,K) = 0._R8
              ELSE
                  GRD%MVLOC(I,J,K) = 0.5_R8*(1._R8-COS(RPI*(ZL2-GRD%DEP(K))/ZDL))
              ENDIF
           ENDDO
          ENDDO
         ENDDO
        ELSE
          WRITE(IOUNERR,*) ' NN_MLD_VLOC = ',NN_MLD_VLOC, ' NOT RECOGNIZED'
          CALL ABOR1('SUEOFSE : UNSUPPORTED NN_MLD_VLOC')
        ENDIF
        DEALLOCATE( MLD )
        DO K=2,PSV3D
            GRD%MVLOC(:,:,((K-1)*GRD%KM+1):(K*GRD%KM))=GRD%MVLOC(:,:,1:GRD%KM)
        ENDDO
     ENDIF

     IF (LL_UNMASK_EOFS) THEN
      DO K=2,GRD%KM
       KP=0
       KP2=0
       KP3=0
       DO J=1,GRD%JM 
        DO I=1,GRD%IM 
         KK=GRD%REG(I,J)
         IF( KK .GT. 0 .AND. ALL(ROS%EVC(KK,K,:).EQ.0._R8) .AND. GRD%MSK(I,J,K).GT.0.9_R8 ) THEN
           KP=KP+1
           IF( ALL(ROS%EVC(KK,K-1,:).EQ.0._R8) ) THEN
              KP2=KP2+1
           ELSE
              KP3=KP3+1
              ROS%EVC(KK,K,:) = ROS%EVC(KK,K-1,:)
              ROS%EVC(KK,K+GRD%KM,:) = ROS%EVC(KK,K+GRD%KM-1,:)
           ENDIF
         ENDIF
        ENDDO
       ENDDO
       WRITE(IOUNLOG,*) ' LL_UNMASK_EOFS -- ', K, KP, KP2, KP3
      ENDDO
     ENDIF

CALL MYFRTPROF_WALL('SUEOF: SET-UP EOF',1)
END SUBROUTINE SUEOF
