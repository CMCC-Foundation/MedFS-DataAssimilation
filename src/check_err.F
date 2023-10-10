      SUBROUTINE CHECK_ERR(IRET)
      INTEGER IRET
      INCLUDE 'netcdf.inc'
      IF (IRET .NE. NF_NOERR) THEN
      PRINT *, NF_STRERROR(IRET)
      STOP
      ENDIF
      END
