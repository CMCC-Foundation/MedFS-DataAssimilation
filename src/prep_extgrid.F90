SUBROUTINE PREP_EXTGRID

USE SET_KND
USE RECFILTER
USE GRD_STR
USE MYFRTPROF, ONLY : MYFRTPROF_WALL
USE IOUNITS, ONLY : IOUNOUT,IOUNLOG
USE RUN
USE ICE
USE LBCLNK

IMPLICIT NONE

INTEGER(I4) :: ICEOPT, KF, JLEV

CALL MYFRTPROF_WALL('PREP_EXTGRID: PREPARE GRID EXTENSION',0)

  IF( NCONF .NE. 202 ) THEN

    WRITE(IOUNLOG,*) ' ALLOCATING CRX/CRY, DIMENSION: ', &
    & GRD%IM,GRD%JM,PSV3D*GRD%KM+PSV2D

    ALLOCATE( CRX(GRD%IM,GRD%JM,PSV3D*GRD%KM+PSV2D),&
    & CRY(GRD%IM,GRD%JM,PSV3D*GRD%KM+PSV2D) )
     !!!!jenny set RF_L from namelist

    IF( .NOT. LL_RFLV ) THEN
       CRX(:,:,:)=RF_L
       CRY(:,:,:)=RF_L
    ELSE

    KF=0
    IF( LL_TS ) THEN
      IF( .NOT. LL_CORRAD_2D ) THEN
         CALL GETNCVAR('Corrad_onx.nc','votemper',&
         & GRD%IM,GRD%JM,GRD%KM,CRX(:,:,KF+1:KF+GRD%KM) )
         CALL GETNCVAR('Corrad_ony.nc','votemper',&
         & GRD%IM,GRD%JM,GRD%KM,CRY(:,:,KF+1:KF+GRD%KM) )
      ELSE
         CALL GETNCVAR('Corrad_onx.nc','votemper',&
         & GRD%IM,GRD%JM,       CRX(:,:,KF+1) )
         CALL GETNCVAR('Corrad_ony.nc','votemper',&
         & GRD%IM,GRD%JM,       CRY(:,:,KF+1) )
         DO JLEV=2,GRD%KM
            CRX(:,:,KF+JLEV) = CRX(:,:,KF+1)
            CRY(:,:,KF+JLEV) = CRY(:,:,KF+1)
         ENDDO
      ENDIF
      KF=KF+GRD%KM
      IF( .NOT. LL_CORRAD_2D ) THEN
         CALL GETNCVAR('Corrad_onx.nc','vosaline',&
         & GRD%IM,GRD%JM,GRD%KM,CRX(:,:,KF+1:KF+GRD%KM) )
         CALL GETNCVAR('Corrad_ony.nc','vosaline',&
         & GRD%IM,GRD%JM,GRD%KM,CRY(:,:,KF+1:KF+GRD%KM) )
      ELSE
         CALL GETNCVAR('Corrad_onx.nc','vosaline',&
         & GRD%IM,GRD%JM,       CRX(:,:,KF+1) )
         CALL GETNCVAR('Corrad_ony.nc','vosaline',&
         & GRD%IM,GRD%JM,       CRY(:,:,KF+1) )
         DO JLEV=2,GRD%KM
            CRX(:,:,KF+JLEV) = CRX(:,:,KF+1)
            CRY(:,:,KF+JLEV) = CRY(:,:,KF+1)
         ENDDO
      ENDIF
      KF=KF+GRD%KM
    ENDIF

    IF( LL_UV ) THEN
      CALL GETNCVAR('Corrad_onx.nc','vozocrtx',&
      & GRD%IM,GRD%JM,GRD%KM,CRX(:,:,KF+1:KF+GRD%KM) )
      CALL GETNCVAR('Corrad_ony.nc','vozocrtx',&
      & GRD%IM,GRD%JM,GRD%KM,CRY(:,:,KF+1:KF+GRD%KM) )
      KF=KF+GRD%KM
      CALL GETNCVAR('Corrad_onx.nc','vomecrty',&
      & GRD%IM,GRD%JM,GRD%KM,CRX(:,:,KF+1:KF+GRD%KM) )
      CALL GETNCVAR('Corrad_ony.nc','vomecrty',&
      & GRD%IM,GRD%JM,GRD%KM,CRY(:,:,KF+1:KF+GRD%KM) )
      KF=KF+GRD%KM
    ENDIF

    IF( LL_SSH) THEN
      CALL GETNCVAR('Corrad_onx.nc','sossheig',&
      & GRD%IM,GRD%JM,CRX(:,:,KF+1) )
      CALL GETNCVAR('Corrad_ony.nc','sossheig',&
      & GRD%IM,GRD%JM,CRY(:,:,KF+1) )
      KF=KF+1
    ENDIF

    IF( LL_SST) THEN
      CALL GETNCVAR('Corrad_onx.nc','sosstsst',&
      & GRD%IM,GRD%JM,CRX(:,:,KF+1) )
      CALL GETNCVAR('Corrad_ony.nc','sosstsst',&
      & GRD%IM,GRD%JM,CRY(:,:,KF+1) )
      KF=KF+1
    ENDIF

    IF( LL_SSS) THEN
      CALL GETNCVAR('Corrad_onx.nc','sosaline',&
      & GRD%IM,GRD%JM,CRX(:,:,KF+1) )
      CALL GETNCVAR('Corrad_ony.nc','sosaline',&
      & GRD%IM,GRD%JM,CRY(:,:,KF+1) )
      KF=KF+1
    ENDIF

    IF( LL_TQ2) THEN
      CALL GETNCVAR('Corrad_onx.nc','soattaml',&
      & GRD%IM,GRD%JM,CRX(:,:,KF+1) )
      CALL GETNCVAR('Corrad_ony.nc','soattaml',&
      & GRD%IM,GRD%JM,CRY(:,:,KF+1) )
      KF=KF+1
      CALL GETNCVAR('Corrad_onx.nc','soatqaml',&
      & GRD%IM,GRD%JM,CRX(:,:,KF+1) )
      CALL GETNCVAR('Corrad_ony.nc','soatqaml',&
      & GRD%IM,GRD%JM,CRY(:,:,KF+1) )
      KF=KF+1
    ENDIF
    IF( LL_ICE) THEN
      CALL GETNCVAR('Corrad_onx.nc','ileadfra',&
      & GRD%IM,GRD%JM,CRX(:,:,KF+1) )
      CALL GETNCVAR('Corrad_ony.nc','ileadfra',&
      & GRD%IM,GRD%JM,CRY(:,:,KF+1) )
      KF=KF+1
      CALL GETNCVAR('Corrad_onx.nc','iicethic',&
      & GRD%IM,GRD%JM,CRX(:,:,KF+1) )
      CALL GETNCVAR('Corrad_ony.nc','iicethic',&
      & GRD%IM,GRD%JM,CRY(:,:,KF+1) )
      KF=KF+1
      CALL GETNCVAR('Corrad_onx.nc','iicevelu',&
      & GRD%IM,GRD%JM,CRX(:,:,KF+1) )
      CALL GETNCVAR('Corrad_ony.nc','iicevelu',&
      & GRD%IM,GRD%JM,CRY(:,:,KF+1) )
      KF=KF+1
      CALL GETNCVAR('Corrad_onx.nc','iicevelv',&
      & GRD%IM,GRD%JM,CRX(:,:,KF+1) )
      CALL GETNCVAR('Corrad_ony.nc','iicevelv',&
      & GRD%IM,GRD%JM,CRY(:,:,KF+1) )
      KF=KF+1
    ENDIF

    WRITE(IOUNLOG,*) ' LINKING BOUNDARIES'

    CALL LBC_LNK( CRX )
    CALL LBC_LNK( CRY )

  ENDIF ! SET FROM NAMELIST
  ENDIF
  IF( LLVARMDT .OR. NCONF .EQ. 202) THEN
     ALLOCATE( CRMX(GRD%IM,GRD%JM),CRMY(GRD%IM,GRD%JM) )
     IF( ZMDT_CLS .LE. 0._R8 ) THEN
       CALL READ_CRZ(GRD%IM,GRD%JM,CRMX,CRMY)
       CALL LBC_LNK( CRMX )
       CALL LBC_LNK( CRMY )
     ELSE
       CRMX = ZMDT_CLS
       CRMY = ZMDT_CLS
     ENDIF
  ENDIF

  ICEOPT = COUNT( (/ LL_ICEREJ, LL_ICEBGFILT, &
  & LL_ICEANFILT, LL_ICERFMASK, LL_ICEREJ_DCOAST/) )

  IF( .NOT. LL_ICEINIT .AND. ICEOPT.GT.0) CALL ICE_INIT

CALL MYFRTPROF_WALL('PREP_EXTGRID: PREPARE GRID EXTENSION',1)

END SUBROUTINE PREP_EXTGRID
