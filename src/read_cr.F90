SUBROUTINE READ_CR(NX,NY,NZ,CXTO,CYTO,CXSO,CYSO)

  USE NETCDF
  USE SET_KND
  USE GRD_STR
  USE IOUNITS
  USE RUN, ONLY : LL_CORRAD_2D
  USE DRWN
  USE MYNETCDF

  IMPLICIT NONE

  INTEGER(I4), INTENT(IN) :: NX, NY, NZ

  INTEGER :: NCID,VARID(2)

  REAL(R8), DIMENSION(NX,NY,NZ),INTENT(OUT) :: CXTO,CYTO,CXSO,CYSO

  REAL(R4), DIMENSION(NX,NY,NZ) :: CXT,CYT,CXS,CYS
  INTEGER(I4), DIMENSION(NX,NY) :: IMASK
  INTEGER(I4) :: X, Y, JLEV

  CALL CHECK( NF90_OPEN('Corrad_onx.nc', NF90_NOWRITE, NCID) )
  CALL CHECK( NF90_INQ_VARID(NCID, 'votemper', VARID(1)) )
  CALL CHECK( NF90_INQ_VARID(NCID, 'vosaline', VARID(2)) )
  IF( .NOT. LL_CORRAD_2D ) THEN
    CALL CHECK( NF90_GET_VAR(NCID, VARID(1),CXTO,START=(/GRD%I0,GRD%J0,1/),COUNT=(/GRD%IM,GRD%JM,GRD%KM/) ) )
    CALL CHECK( NF90_GET_VAR(NCID, VARID(2),CXSO,START=(/GRD%I0,GRD%J0,1/),COUNT=(/GRD%IM,GRD%JM,GRD%KM/) ) )
  ELSE
    CALL CHECK( NF90_GET_VAR(NCID, VARID(1),CXTO(:,:,1),START=(/GRD%I0,GRD%J0,1/),COUNT=(/GRD%IM,GRD%JM,1/) ) )
    CALL CHECK( NF90_GET_VAR(NCID, VARID(2),CXSO(:,:,1),START=(/GRD%I0,GRD%J0,1/),COUNT=(/GRD%IM,GRD%JM,1/) ) )
    DO JLEV=2,GRD%KM
      CXTO(:,:,JLEV)=CXTO(:,:,1)
      CXSO(:,:,JLEV)=CXSO(:,:,1)
    ENDDO
  ENDIF
  CALL CHECK( NF90_CLOSE(NCID) )
  WRITE(IOUNLOG,*) ' CORRAD ON X READ'
  CALL FLUSH(IOUNLOG)

  CALL CHECK( NF90_OPEN('Corrad_ony.nc', NF90_NOWRITE, NCID) )
  CALL CHECK( NF90_INQ_VARID(NCID, 'votemper', VARID(1)) )
  CALL CHECK( NF90_INQ_VARID(NCID, 'vosaline', VARID(2)) )
  IF( .NOT. LL_CORRAD_2D ) THEN
    CALL CHECK( NF90_GET_VAR(NCID, VARID(1),CYTO,START=(/GRD%I0,GRD%J0,1/),COUNT=(/GRD%IM,GRD%JM,GRD%KM/) ) )
    CALL CHECK( NF90_GET_VAR(NCID, VARID(2),CYSO,START=(/GRD%I0,GRD%J0,1/),COUNT=(/GRD%IM,GRD%JM,GRD%KM/) ) )
  ELSE
    CALL CHECK( NF90_GET_VAR(NCID, VARID(1),CYTO(:,:,1),START=(/GRD%I0,GRD%J0,1/),COUNT=(/GRD%IM,GRD%JM,1/) ) )
    CALL CHECK( NF90_GET_VAR(NCID, VARID(2),CYSO(:,:,1),START=(/GRD%I0,GRD%J0,1/),COUNT=(/GRD%IM,GRD%JM,1/) ) )
    DO JLEV=2,GRD%KM
      CYTO(:,:,JLEV)=CYTO(:,:,1)
      CYSO(:,:,JLEV)=CYSO(:,:,1)
    ENDDO
  ENDIF
  CALL CHECK( NF90_CLOSE(NCID) )
  WRITE(IOUNLOG,*) ' CORRAD ON Y READ'
  CALL FLUSH(IOUNLOG)

END SUBROUTINE READ_CR
