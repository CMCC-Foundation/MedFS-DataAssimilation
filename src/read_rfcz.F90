SUBROUTINE READ_RFCZ(NX,NY,SCXT,SCYT,GALXT,GALYT,GBTXT,GBTYT)

  USE SET_KND
  USE NETCDF
  USE MYFRTPROF
  USE RECFILTER, ONLY : LLRFCBOUNDED, ZRFMAXC

  IMPLICIT NONE

  INTEGER(I4) :: NX, NY

  INTEGER(I4) :: NCID,VARID(6)

  REAL(R8), DIMENSION(NX,NY) :: SCXT,SCYT,GALXT,GALYT,GBTXT,GBTYT
  INTEGER(I4) :: X, Y

  CALL MYFRTPROF_WALL('READ_RFCZ: READ RF COEFFS',0)

  CALL CHECK( NF90_OPEN('rf_coeffs_SSH.nc', NF90_NOWRITE, NCID) )
  CALL CHECK( NF90_INQ_VARID(NCID, 'scx', VARID(1)) )
  CALL CHECK( NF90_INQ_VARID(NCID, 'scy', VARID(2)) )
  CALL CHECK( NF90_INQ_VARID(NCID, 'galx', VARID(3)) )
  CALL CHECK( NF90_INQ_VARID(NCID, 'galy', VARID(4)) )
  CALL CHECK( NF90_INQ_VARID(NCID, 'gbtx', VARID(5)) )
  CALL CHECK( NF90_INQ_VARID(NCID, 'gbty', VARID(6)) )

  CALL CHECK( NF90_GET_VAR(NCID, VARID(1),SCXT ) )
  CALL CHECK( NF90_GET_VAR(NCID, VARID(2),SCYT ) )
  CALL CHECK( NF90_GET_VAR(NCID, VARID(3),GALXT ) )
  CALL CHECK( NF90_GET_VAR(NCID, VARID(4),GALYT ) )
  CALL CHECK( NF90_GET_VAR(NCID, VARID(5),GBTXT ) )
  CALL CHECK( NF90_GET_VAR(NCID, VARID(6),GBTYT ) )

  CALL CHECK( NF90_CLOSE(NCID) )

  IF(LLRFCBOUNDED) THEN

     WHERE(GALXT .GT. ZRFMAXC ) GALXT = ZRFMAXC
     WHERE(GALYT .GT. ZRFMAXC ) GALYT = ZRFMAXC
     WHERE(GBTXT .GT. ZRFMAXC ) GBTXT = ZRFMAXC
     WHERE(GBTYT .GT. ZRFMAXC ) GBTYT = ZRFMAXC

  ENDIF

  CALL MYFRTPROF_WALL('READ_RFCZ: READ RF COEFFS',1)

END SUBROUTINE READ_RFCZ
