#undef SHARED_MEMORY
SUBROUTINE MDEOFS

!---------------------------------------------------------------------------
!                                                                          !
!    COPYRIGHT 2006 SRDJAN DOBRICIC, CMCC, BOLOGNA                         !
!    2011 : RECODED FOR ONLINE USE (A.STORTO)                              !  
!                                                                          !
!--------------------------------------------------------------------------- 


  USE SET_KND
  USE GRD_STR
  USE EOF_STR
  USE RUN
  USE MYFRTPROF
  USE IOUNITS

  IMPLICIT NONE

  INTEGER (I4)               :: I, J, K, N
  INTEGER (I4)               :: KI, KJ, KR

  REAL(R8)  :: COVRAD, DST, DAM
  REAL(R8), ALLOCATABLE, DIMENSION (:,:)   :: BV
  REAL(R8), ALLOCATABLE, DIMENSION (:,:)   :: SV
  REAL(R8), ALLOCATABLE, DIMENSION (:)   :: CVR,MLD,KM2
  REAL(R8), ALLOCATABLE, DIMENSION (:,:)   :: DNSR, KM2R, MSKR
  REAL(R8), ALLOCATABLE :: DENS(:,:,:) , MXLD(:,:)

! ARPACK SECTION
  INTEGER(I4)         :: NLAND,NTOTK,MYPOINT(ROS%NREG)
  INTEGER(I4) :: NCVT, NCV2T

  REAL(R8), ALLOCATABLE :: ZFDF(:,:), ZT(:,:,:), ZS(:,:,:)
  INTEGER(I4) :: IAUX, IERR
  INTEGER(I4), PARAMETER :: NUNFD=671
  REAL(R8)    :: ZCF

  CALL MYFRTPROF_WALL('MDEOFS: MODULATE EOFS DEPENDING ON MLD',0)

  NCVT=INT(ROS%NEOF*1.5)
  NCV2T=NCVT*(NCVT+8)

  ALLOCATE ( BV(ROS%KMT,ROS%KMT) )
  ALLOCATE ( SV(GRD%KM,GRD%KM) )
  ALLOCATE ( CVR(ROS%NREG) )
  ALLOCATE ( DNSR(ROS%NREG,GRD%KM) )
  ALLOCATE ( KM2R(ROS%NREG,GRD%KM) )
  ALLOCATE ( MSKR(ROS%NREG,GRD%KM) )
  ALLOCATE (  MLD(ROS%NREG       ) )
  ALLOCATE (  KM2(ROS%NREG       ) )

  ALLOCATE (  DENS(GRD%IM,GRD%JM,GRD%KM) )
  ALLOCATE (  MXLD(GRD%IM,GRD%JM       ) )

  CALL RDDNS( GRD%IM,GRD%JM,GRD%KM,DENS,MXLD )

  IF( LLFLOWDEP ) THEN
      WRITE(IOUNLOG,*)
      WRITE(IOUNLOG,*) ' FLOW-DEPENDENT FACTOR FOR B VARIANCES'
      ALLOCATE ( ZFDF(ROS%NREG,ROS%KMT) )
      ALLOCATE ( ZT(ROS%EOGNX,ROS%EOGNY,GRD%KM) )
      ALLOCATE ( ZS(ROS%EOGNX,ROS%EOGNY,GRD%KM) )
      CALL RDFDF(ROS%EOGNX,ROS%EOGNY,GRD%KM,ZT,ZS)
      K=0
      DO J=1,ROS%EOGNY
        DO I=1,ROS%EOGNX
          K=K+1
          ZFDF(K,1:GRD%KM) = ZT(I,J,:)
          ZFDF(K,GRD%KM+1:2*GRD%KM) = ZS(I,J,:)
        ENDDO
      ENDDO
      DEALLOCATE( ZT, ZS )
      IF( LLFLOWDEP_CALIBR ) THEN
          OPEN(NUNFD, FILE='FD_CALIBR.dat', STATUS='OLD', IOSTAT=IERR)
          IF(IERR.NE.0) CALL ABOR1('MDEOFS: PROBLEM OPENING FD_CALIBR.dat')
          WRITE(IOUNLOG,*)
          WRITE(IOUNLOG,*) ' APPLICATION OF CALIBRATION FACTOR'
          DO J=1,GRD%KM
             READ(NUNFD,*,IOSTAT=IERR) IAUX, ZCF
             IF(IERR.NE.0) CALL ABOR1('MDEOFS: PROBLEM READING FD_CALIBR.dat')
             ZFDF(:,J) = ZFDF(:,J)*ZCF
             WRITE(IOUNLOG,*) ' -- LEVEL ',J,' FACTOR = ',ZCF
             ZFDF(:,GRD%KM+J) = ZFDF(:,GRD%KM+J)*ZCF
          ENDDO
          CLOSE(NUNFD)
          WRITE(IOUNLOG,*)
      ENDIF
  ENDIF

  CVR(:)    = 0._R8
  DNSR(:,:) = 0._R8
  KM2R(:,:) = 0._R8
  KM2 (:  ) = 0._R8
  NLAND     = 0
  MLD(:)    = 0._R8
      

  DO K=1,GRD%KM
   DO J=1,GRD%JM
    DO I=1,GRD%IM

      KR=GRD%REG(I,J)
      IF (KR .GT. 0) THEN
         DAM = GRD%MSK(I,J,K)*GRD%DX(I,J)*GRD%DY(I,J)
         DNSR(KR,K)=DNSR(KR,K) + DENS(I,J,K)*DAM
         KM2R(KR,K)=KM2R(KR,K) + DAM
      ENDIF

    ENDDO
   ENDDO
  ENDDO

   DO J=1,GRD%JM
    DO I=1,GRD%IM
      KR=GRD%REG(I,J)
      IF (KR .GT. 0) THEN
         DAM = GRD%MSK(I,J,1)*GRD%DX(I,J)*GRD%DY(I,J)
         MLD (KR  )=MLD (KR  ) + MXLD(I,J  )*DAM
         KM2 (KR  )=KM2 (KR  ) + DAM
      ENDIF
   ENDDO
  ENDDO

  MSKR = 0._R8
  WHERE( KM2R .GT. 0.001_R8 ) 
     DNSR = DNSR/KM2R
     MSKR = 1._R8
  ENDWHERE
  WHERE( KM2 .GT. 0.001_R8 ) 
     MLD  = MLD /KM2
  ENDWHERE

  DEALLOCATE( MXLD, DENS)
  PRINT*,'FORMED WEIGTHED DENS'

! FIND CORRELATION RADIUS
  DO J=1,ROS%NREG
   
    COVRAD = 0._R8
    IF(MSKR(J,1).GT.0.9) THEN

     DO KI=1,GRD%KM
     
      IF(MSKR(J,KI).EQ.1._R8 .AND. GRD%DEP(KI).LE.ROS%DEP_MAX) THEN
       COVRAD = 0.5_R8 * MAX(1.0_R8,DNSR(J,KI)-DNSR(J,1))
      ENDIF

     ENDDO

    ENDIF

    CVR(J) = COVRAD

  ENDDO

  K=0
  DO J=1,ROS%NREG
    IF(MSKR(J,1).GT.0._R8 .AND. ROS%EVA(J,1) .GT. 0._R8 &
    & .AND. ROS%EVC(J,1,1) .NE. 0._R8)THEN
      K=K+1
      MYPOINT(K)=J
    ENDIF
  ENDDO

  NTOTK=K
  NLAND=ROS%NREG-NTOTK

  WRITE(IOUNLOG,*) ' LAND POINTS ',NLAND
  WRITE(IOUNLOG,*) ' SEA  POINTS ',NTOTK

#ifdef SHARED_MEMORY
!$OMP PARALLEL DEFAULT(SHARED), PRIVATE(K,J,BV,KI,KJ,N,COVRAD,DST,SV)
!$OMP DO SCHEDULE(DYNAMIC)
#endif
  DO K=1,NTOTK
     WRITE(IOUNOUT,*) K, NTOTK
     J=MYPOINT(K)

! FORM BV
      DO KJ=1,ROS%KMT
        DO KI=1,ROS%KMT
          BV(KI,KJ) = SUM( ROS%EVC(J,KI,:) * &
          & ROS%EVA(J,:) * ROS%EVC(J,KJ,:) )
        ENDDO
      ENDDO

! FORM DENSITY COVARIANCES

      COVRAD = 0.5_R8*CVR(J)**2

      DO KJ=1,GRD%KM
        DO KI=1,GRD%KM
          DST = (DNSR(J,KI)-DNSR(J,KJ))**2
          ! IF( ABS(DST) .GT. 100._R8 ) DST = 0._R8
          SV(KI,KJ) = EXP(-DST/COVRAD)
        ENDDO
      ENDDO

! FORM ENTRYWISE PRODUCTS
      DO KJ=1,GRD%KM
        DO KI=1,GRD%KM
          BV(KI         ,KJ         ) = BV(KI         ,KJ         ) * SV(KI,KJ) 
          BV(KI  +GRD%KM,KJ         ) = BV(KI  +GRD%KM,KJ         ) * SV(KI,KJ) 
          BV(KI         ,KJ  +GRD%KM) = BV(KI         ,KJ  +GRD%KM) * SV(KI,KJ) 
          BV(KI  +GRD%KM,KJ  +GRD%KM) = BV(KI  +GRD%KM,KJ  +GRD%KM) * SV(KI,KJ) 
        ENDDO
      ENDDO

      IF( LLFLOWDEP ) THEN
        DO KI=1,ROS%KMT
          BV( KI, KI ) = BV( KI, KI ) * ZFDF(J, KI)  
        ENDDO
      ENDIF

! PARTIAL EIGENVALUE DECOMPOSITION WITH ARPACK

      CALL EIGENV(NCVT, NCV2T,BV,ROS%EVC(J,:,:),ROS%EVA(J,:))

    ENDDO
#ifdef SHARED_MEMORY
!$OMP END DO
!$OMP END PARALLEL
#endif

DEALLOCATE ( CVR )
DEALLOCATE ( BV, SV )
DEALLOCATE ( DNSR, KM2R, MSKR )

IF( ALLOCATED (ZFDF) ) DEALLOCATE(ZFDF)

IF( NCONF .EQ. 105 .OR. LL_WRITE_MDEOFS ) THEN
    WRITE(IOUNOUT,*) ' CALLING WRITE_NCEOF'
    CALL WRITE_NCEOF('EOF_MODULATED.nc',.TRUE.)
ENDIF

! ----------------------------------------------------------------
  CALL MYFRTPROF_WALL('MDEOFS: MODULATE EOFS DEPENDING ON MLD',1)

END SUBROUTINE MDEOFS
