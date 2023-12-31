SUBROUTINE READICOADSFILENC(CFIN,IOUN,NMAX,NOB,KTYP,PAR,DEP,VAL,TIM,LON,LAT,&
& PRF,PLA,INST,BIAS)

! READ OBSERVATION IN REFORMATTED ICOADS FORMAT
!
! A.S. - 26.11.2013

  USE SET_KND
  USE CALENDAR
  USE OBSDEF
  USE NETCDF
  USE RUN , ONLY : ZJULSTART, ZJULEND

  IMPLICIT NONE

  INTEGER(I4), INTENT(IN) :: IOUN,NMAX
  CHARACTER(LEN=*), INTENT(IN) :: CFIN
  INTEGER(I4), INTENT(OUT) :: NOB
  INTEGER(I4), DIMENSION(NMAX),INTENT(OUT) :: KTYP, PAR, PRF, INST
  REAL   (DP), DIMENSION(NMAX),INTENT(OUT) :: TIM
  REAL   (R8), DIMENSION(NMAX),INTENT(OUT) :: DEP,  VAL, LON, LAT, BIAS
  CHARACTER(LEN=8 ), DIMENSION(NMAX),INTENT(OUT) :: PLA

  INTEGER(I4), PARAMETER :: MAXVARS=14
  INTEGER(I4) :: KOBS, NOTW, NFLR, KKTP, KKOP
  INTEGER(I4) :: NCID, ID, JOBS, NOBS
  REAL(DP)    :: ZTIME, ZJULTMP

  REAL(DP), ALLOCATABLE, DIMENSION(:) :: RTIME
  REAL(R8), ALLOCATABLE, DIMENSION(:) :: RLON, RLAT, VALT, &
  & MACROBIAS, MICROBIAS, RERR
  INTEGER(I4), ALLOCATABLE, DIMENSION(:) :: ISTAT, IREPT, ISTATION
  CHARACTER(LEN=8 ), ALLOCATABLE, DIMENSION(:) :: CTMP

  WRITE(IOUN,*) ' ^^^ READICOADSFILENC PROCESSING FILE ',TRIM(CFIN)

  KOBS=0
  KKOP=0
  NOTW=0
  NFLR=0

  CALL YMDS2JU( 1950, 1, 1, 0._DP, ZJULTMP )

  CALL CHECK( NF90_OPEN(TRIM(CFIN), NF90_NOWRITE, NCID ) )
  CALL CHECK( NF90_INQ_DIMID( NCID, 'obs', ID ) )
  CALL CHECK( NF90_INQUIRE_DIMENSION( NCID, ID, LEN = NOBS ) )

  ALLOCATE( RLON(NOBS)  , RLAT(NOBS), ISTAT(NOBS) )
  ALLOCATE( RTIME(NOBS) , VALT(NOBS), MACROBIAS(NOBS) )
  ALLOCATE( RERR(NOBS)  , MICROBIAS(NOBS), CTMP(NOBS), IREPT(NOBS) )
  ALLOCATE( ISTATION(NOBS) )

  CALL CHECK( NF90_INQ_VARID( NCID, 'longitude', ID ) )
  CALL CHECK( NF90_GET_VAR  ( NCID, ID, RLON ) )
  CALL CHECK( NF90_INQ_VARID( NCID, 'latitude', ID ) )
  CALL CHECK( NF90_GET_VAR  ( NCID, ID, RLAT ) )
  CALL CHECK( NF90_INQ_VARID( NCID, 'julian_time', ID ) )
  CALL CHECK( NF90_GET_VAR  ( NCID, ID, RTIME) )
  CALL CHECK( NF90_INQ_VARID( NCID, 'sst_value', ID ) )
  CALL CHECK( NF90_GET_VAR  ( NCID, ID, VALT ) )
  CALL CHECK( NF90_INQ_VARID( NCID, 'macro_bias', ID ) )
  CALL CHECK( NF90_GET_VAR  ( NCID, ID, MACROBIAS ) )
  CALL CHECK( NF90_INQ_VARID( NCID, 'micro_bias', ID ) )
  CALL CHECK( NF90_GET_VAR  ( NCID, ID, MICROBIAS ) )
  CALL CHECK( NF90_INQ_VARID( NCID, 'status', ID ) )
  CALL CHECK( NF90_GET_VAR  ( NCID, ID, ISTAT) )
  CALL CHECK( NF90_INQ_VARID( NCID, 'total_err', ID ) )
  CALL CHECK( NF90_GET_VAR  ( NCID, ID, RERR ) )
  CALL CHECK( NF90_INQ_VARID( NCID, 'platform', ID ) )
  CALL CHECK( NF90_GET_VAR  ( NCID, ID, CTMP ) )
  CALL CHECK( NF90_INQ_VARID( NCID, 'report_type', ID ) )
  CALL CHECK( NF90_GET_VAR  ( NCID, ID, IREPT ) )
  CALL CHECK( NF90_INQ_VARID( NCID, 'station_id', ID ) )
  CALL CHECK( NF90_GET_VAR  ( NCID, ID, ISTATION ) )

  CALL CHECK( NF90_CLOSE ( NCID ) )

  BIAS = 0._R8

  ! Select
  CYOBS : DO JOBS=1, NOBS
     ZTIME = RTIME(JOBS)
     IF( ZTIME .LT. ZJULSTART .OR. ZTIME .GT. ZJULEND) THEN
       NOTW = NOTW + 1
       CYCLE CYOBS
     ENDIF
     IF( ISTAT(JOBS) .NE. 2 ) THEN
       NFLR = NFLR + 1
       CYCLE CYOBS
     ENDIF
     KOBS = KOBS + 1  
     LON(KOBS) = RLON(JOBS)
     LAT(KOBS) = RLAT(JOBS)
     PAR(KOBS) = KKSST
     INST(KOBS) = ISTATION(JOBS)
     SELECT CASE(ISTATION(JOBS))
       CASE(0:5,9,11,12,20) 
            KTYP(KOBS) = KKXBT
       CASE(6:8,13:16,19) 
            KTYP(KOBS) = KKBUOY
       CASE(10,17) 
            KTYP(KOBS) = KKTESAC
       CASE(18,21) 
            KTYP(KOBS) = KKARGO
       CASE DEFAULT 
            KTYP(KOBS) = KKXBT
     ENDSELECT
     DEP(KOBS) = RERR(JOBS)
     VAL(KOBS) = VALT(JOBS)
     IF( MACROBIAS(JOBS) .GT. -99._R8 ) THEN
           VAL(KOBS) = VAL(KOBS) + MACROBIAS(JOBS)
           BIAS(KOBS) = BIAS(KOBS) + MACROBIAS(JOBS)
     ENDIF
     IF( MICROBIAS(JOBS) .GT. -99._R8 ) THEN
           VAL(KOBS) = VAL(KOBS) + MICROBIAS(JOBS)
           BIAS(KOBS) = BIAS(KOBS) + MICROBIAS(JOBS)
     ENDIF
     TIM(KOBS) = ZTIME - ZJULTMP
     PRF(KOBS) = KOBS
     IF(CTMP(JOBS)(5:8) .NE. 'NULL') THEN
       PLA(KOBS)(1:8) = CTMP(JOBS)(1:8)
     ELSE
       PLA(KOBS)(1:8) = '00000000'
     ENDIF
   
  ENDDO CYOBS

  NOB = KOBS

  WRITE(IOUN,*)
  WRITE(IOUN,*) ' ^^^ REPORT '
  WRITE(IOUN,*) ' TOTAL NUMBER OF OBS     :', NOBS
  WRITE(IOUN,*) ' NUMBER OF BAD  OBS COO  :', KKOP
  WRITE(IOUN,*) ' NUMBER OF GOOD OBSERV.  :', KOBS
  WRITE(IOUN,*) ' NUMBER OF P OUT OF TIME :', NOTW
  WRITE(IOUN,*) ' NUMBER OF BAD FLAGS OBS :', NFLR
  WRITE(IOUN,*) 

END SUBROUTINE READICOADSFILENC
