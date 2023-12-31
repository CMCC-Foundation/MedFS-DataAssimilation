SUBROUTINE SET_DIM

  USE SET_KND
  USE NETCDF

  IMPLICIT NONE

    INTEGER(I4) :: NCID, IDVAR

    CALL CHECK( NF90_OPEN('GRID.nc', NF90_NOWRITE, NCID) )
    CALL CHECK( NF90_INQ_DIMID (NCID, 'x', IDVAR) )
    CALL CHECK( NF90_INQUIRE_DIMENSION (NCID, IDVAR, LEN = IM) )
    CALL CHECK( NF90_INQ_DIMID (NCID, 'y', IDVAR) )
    CALL CHECK( NF90_INQUIRE_DIMENSION (NCID, IDVAR, LEN = JM) )
    CALL CHECK( NF90_INQ_DIMID (NCID, 'z', IDVAR) )
    CALL CHECK( NF90_INQUIRE_DIMENSION (NCID, IDVAR, LEN = KM) )
    CALL CHECK( NF90_CLOSE(NCID) )

    WRITE(IOUNLOG,*) ' DIMENSIONS : ',IM,JM,KM

END SUBROUTINE SET_DIM
