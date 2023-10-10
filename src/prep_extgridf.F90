SUBROUTINE PREP_EXTGRIDF

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

CALL MYFRTPROF_WALL('PREP_EXTGRIDF: PREPARE GRID EXTENSION',0)

  IF( NCONF .NE. 202 ) THEN

    WRITE(IOUNLOG,*) ' ALLOCATING CRFX/CRFY, DIMENSION: ', &
    & GRD%IM,GRD%JM,PSV3D*GRD%KM+PSV2D

    ALLOCATE( CRFX(GRD%IM,GRD%JM,PSV3D*GRD%KM+PSV2D),&
    & CRFY(GRD%IM,GRD%JM,PSV3D*GRD%KM+PSV2D) )

    KF=0
    IF( LL_TS ) THEN
      IF( .NOT. LL_CORRAD_2D ) THEN
         CALL GETNCVAR('Corrad_fd_onx.nc','votemper',&
         & GRD%IM,GRD%JM,GRD%KM,CRFX(:,:,KF+1:KF+GRD%KM) )
         CALL GETNCVAR('Corrad_fd_ony.nc','votemper',&
         & GRD%IM,GRD%JM,GRD%KM,CRFY(:,:,KF+1:KF+GRD%KM) )
      ELSE
         CALL GETNCVAR('Corrad_fd_onx.nc','votemper',&
         & GRD%IM,GRD%JM,       CRFX(:,:,KF+1) )
         CALL GETNCVAR('Corrad_fd_ony.nc','votemper',&
         & GRD%IM,GRD%JM,       CRFY(:,:,KF+1) )
         DO JLEV=2,GRD%KM
            CRFX(:,:,KF+JLEV) = CRFX(:,:,KF+1)
            CRFY(:,:,KF+JLEV) = CRFY(:,:,KF+1)
         ENDDO
      ENDIF
      KF=KF+GRD%KM
      IF( .NOT. LL_CORRAD_2D ) THEN
         CALL GETNCVAR('Corrad_fd_onx.nc','vosaline',&
         & GRD%IM,GRD%JM,GRD%KM,CRFX(:,:,KF+1:KF+GRD%KM) )
         CALL GETNCVAR('Corrad_fd_ony.nc','vosaline',&
         & GRD%IM,GRD%JM,GRD%KM,CRFY(:,:,KF+1:KF+GRD%KM) )
      ELSE
         CALL GETNCVAR('Corrad_fd_onx.nc','vosaline',&
         & GRD%IM,GRD%JM,       CRFX(:,:,KF+1) )
         CALL GETNCVAR('Corrad_fd_ony.nc','vosaline',&
         & GRD%IM,GRD%JM,       CRFY(:,:,KF+1) )
         DO JLEV=2,GRD%KM
            CRFX(:,:,KF+JLEV) = CRFX(:,:,KF+1)
            CRFY(:,:,KF+JLEV) = CRFY(:,:,KF+1)
         ENDDO
      ENDIF
      KF=KF+GRD%KM
    ENDIF

    IF( LL_UV ) THEN
      CALL GETNCVAR('Corrad_fd_onx.nc','vozocrtx',&
      & GRD%IM,GRD%JM,GRD%KM,CRFX(:,:,KF+1:KF+GRD%KM) )
      CALL GETNCVAR('Corrad_fd_ony.nc','vozocrtx',&
      & GRD%IM,GRD%JM,GRD%KM,CRFY(:,:,KF+1:KF+GRD%KM) )
      KF=KF+GRD%KM
      CALL GETNCVAR('Corrad_fd_onx.nc','vomecrty',&
      & GRD%IM,GRD%JM,GRD%KM,CRFX(:,:,KF+1:KF+GRD%KM) )
      CALL GETNCVAR('Corrad_fd_ony.nc','vomecrty',&
      & GRD%IM,GRD%JM,GRD%KM,CRFY(:,:,KF+1:KF+GRD%KM) )
      KF=KF+GRD%KM
    ENDIF

    IF( LL_SSH) THEN
      CALL GETNCVAR('Corrad_fd_onx.nc','sossheig',&
      & GRD%IM,GRD%JM,CRFX(:,:,KF+1) )
      CALL GETNCVAR('Corrad_fd_ony.nc','sossheig',&
      & GRD%IM,GRD%JM,CRFY(:,:,KF+1) )
      KF=KF+1
    ENDIF

    IF( LL_SST) THEN
      CALL GETNCVAR('Corrad_fd_onx.nc','sosstsst',&
      & GRD%IM,GRD%JM,CRFX(:,:,KF+1) )
      CALL GETNCVAR('Corrad_fd_ony.nc','sosstsst',&
      & GRD%IM,GRD%JM,CRFY(:,:,KF+1) )
      KF=KF+1
    ENDIF

    IF( LL_SSS) THEN
      CALL GETNCVAR('Corrad_fd_onx.nc','sosaline',&
      & GRD%IM,GRD%JM,CRFX(:,:,KF+1) )
      CALL GETNCVAR('Corrad_fd_ony.nc','sosaline',&
      & GRD%IM,GRD%JM,CRFY(:,:,KF+1) )
      KF=KF+1
    ENDIF

    IF( LL_TQ2) THEN
      CALL GETNCVAR('Corrad_fd_onx.nc','soattaml',&
      & GRD%IM,GRD%JM,CRFX(:,:,KF+1) )
      CALL GETNCVAR('Corrad_fd_ony.nc','soattaml',&
      & GRD%IM,GRD%JM,CRFY(:,:,KF+1) )
      KF=KF+1
      CALL GETNCVAR('Corrad_fd_onx.nc','soatqaml',&
      & GRD%IM,GRD%JM,CRFX(:,:,KF+1) )
      CALL GETNCVAR('Corrad_fd_ony.nc','soatqaml',&
      & GRD%IM,GRD%JM,CRFY(:,:,KF+1) )
      KF=KF+1
    ENDIF

    IF( LL_ICE) THEN
      CALL GETNCVAR('Corrad_fd_onx.nc','ileadfra',&
      & GRD%IM,GRD%JM,CRFX(:,:,KF+1) )
      CALL GETNCVAR('Corrad_fd_ony.nc','ileadfra',&
      & GRD%IM,GRD%JM,CRFY(:,:,KF+1) )
      KF=KF+1
      CALL GETNCVAR('Corrad_fd_onx.nc','iicethic',&
      & GRD%IM,GRD%JM,CRFX(:,:,KF+1) )
      CALL GETNCVAR('Corrad_fd_ony.nc','iicethic',&
      & GRD%IM,GRD%JM,CRFY(:,:,KF+1) )
      KF=KF+1
      CALL GETNCVAR('Corrad_fd_onx.nc','iicevelu',&
      & GRD%IM,GRD%JM,CRFX(:,:,KF+1) )
      CALL GETNCVAR('Corrad_fd_ony.nc','iicevelu',&
      & GRD%IM,GRD%JM,CRFY(:,:,KF+1) )
      KF=KF+1
      CALL GETNCVAR('Corrad_fd_onx.nc','iicevelv',&
      & GRD%IM,GRD%JM,CRFX(:,:,KF+1) )
      CALL GETNCVAR('Corrad_fd_ony.nc','iicevelv',&
      & GRD%IM,GRD%JM,CRFY(:,:,KF+1) )
      KF=KF+1
    ENDIF

    WRITE(IOUNLOG,*) ' LINKING BOUNDARIES'

    CALL LBC_LNK( CRFX )
    CALL LBC_LNK( CRFY )

  ENDIF

CALL MYFRTPROF_WALL('PREP_EXTGRIDF: PREPARE GRID EXTENSION',1)

END SUBROUTINE PREP_EXTGRIDF
