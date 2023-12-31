SUBROUTINE READ_ARGLEVS(CFILENAME,NUMLEVS)

USE SET_KND
USE NETCDF

IMPLICIT NONE

CHARACTER(LEN=*) , INTENT(IN) :: CFILENAME
INTEGER(I4),INTENT(OUT) :: NUMLEVS
INTEGER(I4) :: STAT, NCID, IDVAR

    STAT = NF90_OPEN(CFILENAME, NF90_NOWRITE, NCID)
    IF (STAT /= NF90_NOERR) CALL HANDLE_ERR(STAT)
    STAT = NF90_INQ_DIMID (NCID, 'N_LEVELS', IDVAR)
    IF (STAT /= NF90_NOERR) CALL HANDLE_ERR(STAT)
    STAT = NF90_INQUIRE_DIMENSION (NCID, IDVAR, LEN = NUMLEVS)
    IF (STAT /= NF90_NOERR) CALL HANDLE_ERR(STAT)
    STAT = NF90_CLOSE(NCID)
    IF (STAT /= NF90_NOERR) CALL HANDLE_ERR(STAT)

END SUBROUTINE READ_ARGLEVS
