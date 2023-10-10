SUBROUTINE READCORFILE(CFILENAME,NUMLEVS,IOUN,MAXOBS,LLSTRICTPOS,LLSTRICTTIME,&
         & LSTRICTPSU,LSTRICTTEM,LGLOBTQF,LGLOBSQF,LLCHECKTIME, LVTHIN,LVAVG,&
         & NTOTOBS,OBSTYPE, OBSPAR,OBSDEP,OBSVAL, OBSTIME, &
         & OBSLON,OBSLAT,OBSPROF,OBSPLAT,CINST)

!
! -> READ IN OBSERVATIONS FROM CORA
!    FROM NATIVE NETCDF FORMAT
!
!    ANDREA STORTO _2013-05_
!

USE SET_KND
USE OBSDEF
USE RUN
USE NETCDF
USE CALENDAR, ONLY : TIME_DIFF,YMDS2JU

IMPLICIT NONE

CHARACTER(LEN=*) , INTENT(IN) :: CFILENAME
CHARACTER(LEN=2) , INTENT(IN) :: CINST
INTEGER(I4),INTENT(IN) :: NUMLEVS,IOUN, MAXOBS
INTEGER(I4),INTENT(OUT) :: NTOTOBS
!... SWITCHES
LOGICAL, INTENT(IN) :: LLSTRICTPOS,LSTRICTPSU,LLSTRICTTIME,&
         & LSTRICTTEM,LGLOBTQF,LGLOBSQF,LLCHECKTIME, LVTHIN, LVAVG

REAL   (DP) :: DATE_START, DATE_END, ZJULREF
INTEGER(I4) :: IYY,IMM,IDD,IHH,IMI,ISS
REAL   (DP) :: ZSS,ZANDATELOC

REAL   (R8) :: X, G, P, ZDENS, DPS_T, DPS_S, OVAL, ZTC(2)
REAL   (R8), PARAMETER :: DEG2RAD = 0.01745329_R8
INTEGER(I4) :: JI, JJ, NTC(2)

INTEGER(I4) :: KERR(59),KOK(5),KCO(10)
INTEGER(I4) :: STAT, NCID, IDVAR, JOBS
INTEGER(I4) :: NDTSTR,NPROFS,NLEVELS,KDTSTR
INTEGER(I4) :: JLEV,JPROF,ITYPE

!... OBS STORAGE
INTEGER(I4), PARAMETER :: MAXPROFS = 50000
INTEGER(I4), DIMENSION(MAXOBS),INTENT(OUT) :: OBSTYPE, OBSPAR, OBSPROF
REAL   (DP), DIMENSION(MAXOBS),INTENT(OUT) :: OBSTIME
REAL   (R8), DIMENSION(MAXOBS),INTENT(OUT) :: OBSDEP,  OBSVAL, OBSLON,OBSLAT
CHARACTER(LEN=8 ), DIMENSION(MAXOBS),INTENT(OUT) :: OBSPLAT
CHARACTER(LEN=40) :: CREFDT, CAUX

REAL   (R8)      , ALLOCATABLE, DIMENSION(:,:) :: TEMP,SALC,DEPTHC
REAL   (R8)      , ALLOCATABLE, DIMENSION(:,:) :: TEMPA,SALCA,DEPTHCA
REAL   (R8)      , ALLOCATABLE, DIMENSION(:)   :: LATS,LONS
REAL   (DP)      , ALLOCATABLE, DIMENSION(:)   :: TIME
CHARACTER(LEN=4 ), ALLOCATABLE, DIMENSION(:)   :: INST
CHARACTER(LEN=8 ), ALLOCATABLE, DIMENSION(:)   :: PLANO

!!! CHARACTER(LEN=NUMLEVS), DIMENSION(MAXPROFS) :: TQC, SQC
!!! CHARACTER(LEN=MAXPROFS) :: TGQC,SGQC,POSQC,TIMQC

INTEGER(4), DIMENSION(NUMLEVS,MAXPROFS) :: TQC, SQC
INTEGER(4), DIMENSION(MAXPROFS) :: TGQC,SGQC,POSQC,TIMQC

INTEGER(4)  :: IDIFFER1, IDIFFER2,IRETC

INTEGER(4)  :: KPROF

INTEGER(4)       :: planoLength
CHARACTER(LEN=8) :: plano_tmp

LOGICAL :: LPROFT, LPROFS, LOBS, LLPROFOK
LOGICAL :: LXBT

!  INITIALIZE
KCO=0
KERR=0
KOK =0
JOBS=0
LXBT = .FALSE.

WRITE(IOUN,*)
WRITE(IOUN,*) '*** READCORFILE : READING NETCDF FILE '//TRIM(CFILENAME)
WRITE(IOUN,*)

IF( LVAVG .AND. LVTHIN ) CALL ABOR1('READCORFILE : LVAVG .AND. LVTHIN BOTH TRUE')

    STAT = NF90_OPEN(CFILENAME, NF90_NOWRITE, NCID)
    IF (STAT /= NF90_NOERR) CALL HANDLE_ERR(STAT)

    STAT = NF90_INQ_VARID (NCID, 'PSAL', IDVAR)
    IF( STAT .NE.  NF90_NOERR ) LXBT = .TRUE.

    STAT = NF90_INQ_DIMID (NCID, 'DATE_TIME', IDVAR)
    IF (STAT /= NF90_NOERR) THEN
       STAT = NF90_INQ_DIMID (NCID, 'TIME', IDVAR)
    ENDIF
    IF (STAT /= NF90_NOERR) CALL HANDLE_ERR(STAT)
    STAT = NF90_INQUIRE_DIMENSION (NCID, IDVAR, LEN = KDTSTR)
    IF (STAT /= NF90_NOERR) CALL HANDLE_ERR(STAT)
    STAT = NF90_INQ_DIMID (NCID, 'N_PROF', IDVAR)
    IF (STAT /= NF90_NOERR) THEN
       STAT = NF90_INQ_DIMID (NCID, 'POSITION', IDVAR)
    ENDIF
    IF (STAT /= NF90_NOERR) CALL HANDLE_ERR(STAT)
    STAT = NF90_INQUIRE_DIMENSION (NCID, IDVAR, LEN = NPROFS)
    IF (STAT /= NF90_NOERR) CALL HANDLE_ERR(STAT)
    STAT = NF90_INQ_DIMID (NCID, 'N_LEVELS', IDVAR)
    IF (STAT /= NF90_NOERR) THEN
       STAT = NF90_INQ_DIMID (NCID, 'DEPTH', IDVAR)
    ENDIF
    IF (STAT /= NF90_NOERR) CALL HANDLE_ERR(STAT)
    STAT = NF90_INQUIRE_DIMENSION (NCID, IDVAR, LEN = NLEVELS)
    IF (STAT /= NF90_NOERR) CALL HANDLE_ERR(STAT)
    WRITE(IOUN,*) 'NETCDF DIMENSIONS READ'

    IF(NPROFS.GT.MAXPROFS) &
    & CALL ABOR1('READCORFILE : MAXPROFS MUST BE INCREASED!')
    IF(NLEVELS.GT.NUMLEVS) &
    & CALL ABOR1('READCORFILE : NUMLEVS DOES NOT MATCH!')

    STAT = NF90_INQ_VARID (NCID, 'REFERENCE_DATE_TIME', IDVAR)
    IF (STAT /= NF90_NOERR) THEN
       !STAT = NF90_INQ_VARID (NCID, 'TIME', IDVAR)
       !IF (STAT /= NF90_NOERR) CALL HANDLE_ERR(STAT)
       !STAT = NF90_GET_ATT (NCID, IDVAR, 'units', CAUX)
       !IF (STAT /= NF90_NOERR) CALL HANDLE_ERR(STAT)
       !CREFDT(1:19) = '1950-01-01T00:00:00'
       IYY=1950 ; IMM=01 ; IDD=01 ; IHH=0 ; IMI=0 ; ISS=0
    ELSE
       STAT = NF90_GET_VAR (NCID,IDVAR,CREFDT(1:KDTSTR))
       IF (STAT /= NF90_NOERR) CALL HANDLE_ERR(STAT)
       WRITE(IOUN,*) 'REF DATE IS ', CREFDT(1:KDTSTR)
       READ(CREFDT(1:4),*) IYY
       READ(CREFDT(5:6),*) IMM
       READ(CREFDT(7:8),*) IDD
       READ(CREFDT(9:10),*) IHH
       READ(CREFDT(11:12),*) IMI
       READ(CREFDT(13:14),*) ISS
    ENDIF
    WRITE(IOUN,*) 'REFERENCE_DATE_TIME READ'

    IF(LLCHECKTIME) THEN

      ZSS = REAL(IHH,KIND=DP)/REAL(24,KIND=DP)  + &
      & REAL(IMI,KIND=DP)/REAL(24*60,KIND=DP)   + &
      & REAL(ISS,KIND=DP)/REAL(24*60*60,KIND=DP)

      CALL YMDS2JU(IYY,IMM,IDD,ZSS,ZJULREF)

      DATE_START = ZJULSTART - ZJULREF
      DATE_END   = ZJULEND   - ZJULREF

      WRITE(IOUN,*) ' DATESTART : ',DATE_START
      WRITE(IOUN,*) ' DATEEND   : ',DATE_END
      WRITE(IOUN,*) ' ... JULIAN DAYS SINCE '//CREFDT(1:10)

    ENDIF

    ALLOCATE(LATS(NPROFS),LONS(NPROFS),TIME(NPROFS),INST(NPROFS),&
           & PLANO(NPROFS))

    STAT = NF90_INQ_VARID (NCID, 'LATITUDE', IDVAR)
    IF (STAT /= NF90_NOERR) CALL HANDLE_ERR(STAT)
    STAT = NF90_GET_VAR (NCID,IDVAR,LATS)
    IF (STAT /= NF90_NOERR) CALL HANDLE_ERR(STAT)
    WRITE(IOUN,*) 'LATITUDE READ'

    STAT = NF90_INQ_VARID (NCID, 'LONGITUDE', IDVAR)
    IF (STAT /= NF90_NOERR) CALL HANDLE_ERR(STAT)
    STAT = NF90_GET_VAR (NCID,IDVAR,LONS)
    IF (STAT /= NF90_NOERR) CALL HANDLE_ERR(STAT)
    WRITE(IOUN,*) 'LONGITUDE READ'

    STAT = NF90_INQ_VARID (NCID, 'JULD', IDVAR)
    IF (STAT /= NF90_NOERR) THEN
        STAT = NF90_INQ_VARID (NCID, 'TIME', IDVAR)
    ENDIF
    IF (STAT /= NF90_NOERR) CALL HANDLE_ERR(STAT)
    STAT = NF90_GET_VAR (NCID,IDVAR,TIME)
    IF (STAT /= NF90_NOERR) CALL HANDLE_ERR(STAT)
    WRITE(IOUN,*) 'JULD READ'

    STAT = NF90_INQ_VARID (NCID, 'WMO_INST_TYPE', IDVAR)
    IF (STAT /= NF90_NOERR) THEN
       INST(:)(1:4) = ' 999'
    ELSE
       STAT = NF90_GET_VAR (NCID,IDVAR,INST)
       IF (STAT /= NF90_NOERR) CALL HANDLE_ERR(STAT)
       WRITE(IOUN,*) 'WMO_INST_TYPE READ'
    ENDIF

!    STAT = NF90_INQ_VARID (NCID, 'PLATFORM_NUMBER', IDVAR)
!    IF (STAT /= NF90_NOERR) THEN
!        STAT = NF90_INQ_VARID (NCID, 'DC_REFERENCE', IDVAR)
!    ENDIF
!    IF (STAT /= NF90_NOERR) THEN
!       PLANO(:)(1:8) = '00000000'
!    ELSE
!        STAT = NF90_GET_VAR (NCID,IDVAR,PLANO)
!        IF (STAT /= NF90_NOERR) CALL HANDLE_ERR(STAT)
!        WRITE(IOUN,*) 'PLATFORM_NUMBER READ'
!    ENDIF

    STAT = NF90_INQUIRE_ATTRIBUTE (NCID, NF90_GLOBAL, 'wmo_platform_code', len=planoLength)
    IF (STAT /= NF90_NOERR) THEN
       PLANO(:)(1:8) = '00000000'
    ELSE
       
       STAT = NF90_GET_ATT (NCID, NF90_GLOBAL, 'wmo_platform_code', plano_tmp)
       PLANO(:)(1:planoLength) = plano_tmp
       WRITE(IOUN,*) 'PLATFORM_NUMBER READ'
    ENDIF


    IF(LLSTRICTPOS) THEN
          STAT = NF90_INQ_VARID (NCID, 'POSITION_QC', IDVAR)
          IF (STAT /= NF90_NOERR) CALL HANDLE_ERR(STAT)
          STAT = NF90_GET_VAR (NCID,IDVAR,POSQC(1:NPROFS))
          IF (STAT /= NF90_NOERR) CALL HANDLE_ERR(STAT)
          WRITE(IOUN,*) 'POSITION_QC READ'
    ENDIF

    IF(LLSTRICTTIME) THEN
          STAT = NF90_INQ_VARID (NCID, 'JULD_QC', IDVAR)
          IF (STAT /= NF90_NOERR) THEN
              STAT = NF90_INQ_VARID (NCID, 'TIME_QC', IDVAR)
          ENDIF
          IF (STAT /= NF90_NOERR) CALL HANDLE_ERR(STAT)
          STAT = NF90_GET_VAR (NCID,IDVAR,TIMQC(1:NPROFS))
          IF (STAT /= NF90_NOERR) CALL HANDLE_ERR(STAT)
          WRITE(IOUN,*) 'JULD_QC READ'
    ENDIF

    IF(LGLOBTQF) THEN
          STAT = NF90_INQ_VARID (NCID, 'PROFILE_TEMP_QC', IDVAR)
          IF (STAT /= NF90_NOERR) CALL HANDLE_ERR(STAT)
          STAT = NF90_GET_VAR (NCID,IDVAR,TGQC(1:NPROFS))
          IF (STAT /= NF90_NOERR) CALL HANDLE_ERR(STAT)
          WRITE(IOUN,*) 'PROFILE_POTM_QC READ'
    ENDIF

    IF(LGLOBSQF .AND. .NOT. LXBT) THEN
          STAT = NF90_INQ_VARID (NCID, 'PROFILE_PSAL_QC', IDVAR)
          IF (STAT /= NF90_NOERR) CALL HANDLE_ERR(STAT)
          STAT = NF90_GET_VAR (NCID,IDVAR,SGQC(1:NPROFS))
          IF (STAT /= NF90_NOERR) CALL HANDLE_ERR(STAT)
          WRITE(IOUN,*) 'PROFILE_PSAL_QC READ'
    ENDIF

    ALLOCATE(TEMP(NLEVELS,NPROFS),SALC(NLEVELS,NPROFS),DEPTHC(NLEVELS,NPROFS))
    ALLOCATE(TEMPA(NLEVELS,NPROFS),SALCA(NLEVELS,NPROFS),DEPTHCA(NLEVELS,NPROFS))

    STAT = NF90_INQ_VARID (NCID, 'TEMP', IDVAR)
    IF (STAT /= NF90_NOERR) CALL HANDLE_ERR(STAT)
    STAT = NF90_GET_VAR (NCID,IDVAR,TEMP)
    IF (STAT /= NF90_NOERR) CALL HANDLE_ERR(STAT)
    WRITE(IOUN,*) 'TEMP READ'

    STAT = NF90_INQ_VARID (NCID, 'TEMP_ADJUSTED', IDVAR)
    IF (STAT /= NF90_NOERR) THEN
       TEMPA = -999._R8
    ELSE
       STAT = NF90_GET_VAR (NCID,IDVAR,TEMPA)
       IF (STAT /= NF90_NOERR) CALL HANDLE_ERR(STAT)
       WRITE(IOUN,*) 'TEMP_ADJUSTED READ'
    ENDIF

    IF( .NOT. LXBT ) THEN
      STAT = NF90_INQ_VARID (NCID, 'PSAL', IDVAR)
      IF (STAT /= NF90_NOERR) CALL HANDLE_ERR(STAT)
      STAT = NF90_GET_VAR (NCID,IDVAR,SALC)
      IF (STAT /= NF90_NOERR) CALL HANDLE_ERR(STAT)
      WRITE(IOUN,*) 'PSAL READ'
      STAT = NF90_INQ_VARID (NCID, 'PSAL_ADJUSTED', IDVAR)
      IF (STAT /= NF90_NOERR) THEN
          SALCA = -999._R8
      ELSE
          STAT = NF90_GET_VAR (NCID,IDVAR,SALCA)
          IF (STAT /= NF90_NOERR) CALL HANDLE_ERR(STAT)
          WRITE(IOUN,*) 'PSAL_ADJUSTED READ'
      ENDIF
    ENDIF

    STAT = NF90_INQ_VARID (NCID, 'DEPH', IDVAR)
    IF( STAT .EQ. NF90_NOERR ) THEN
      STAT = NF90_GET_VAR (NCID,IDVAR,DEPTHC)
      IF (STAT /= NF90_NOERR) CALL HANDLE_ERR(STAT)
      WRITE(IOUN,*) 'DEPH READ'
      STAT = NF90_INQ_VARID (NCID, 'DEPH_ADJUSTED', IDVAR)
      IF (STAT /= NF90_NOERR) THEN
          DEPTHCA = DEPTHC
      ELSE
          STAT = NF90_GET_VAR (NCID,IDVAR,DEPTHCA)
          IF (STAT /= NF90_NOERR) CALL HANDLE_ERR(STAT)
          WRITE(IOUN,*) 'DEPH_ADJUSTED READ'
      ENDIF
    ELSEIF( NF90_INQ_VARID (NCID, 'DEPTH', IDVAR) .EQ. NF90_NOERR ) THEN
      STAT = NF90_GET_VAR (NCID,IDVAR,DEPTHC)
      IF (STAT /= NF90_NOERR) CALL HANDLE_ERR(STAT)
      WRITE(IOUN,*) 'DEPH READ'
      STAT = NF90_INQ_VARID (NCID, 'DEPTH_ADJUSTED', IDVAR)
      IF (STAT /= NF90_NOERR) THEN
          DEPTHCA = DEPTHC
      ELSE
          STAT = NF90_GET_VAR (NCID,IDVAR,DEPTHCA)
          IF (STAT /= NF90_NOERR) CALL HANDLE_ERR(STAT)
          WRITE(IOUN,*) 'DEPH_ADJUSTED READ'
      ENDIF
    ELSE
      STAT = NF90_INQ_VARID (NCID, 'PRES', IDVAR)
      IF (STAT /= NF90_NOERR) CALL HANDLE_ERR(STAT)
      STAT = NF90_GET_VAR (NCID,IDVAR,DEPTHC)
      IF (STAT /= NF90_NOERR) CALL HANDLE_ERR(STAT)
      WRITE(IOUN,*) 'PRES READ'
      STAT = NF90_INQ_VARID (NCID, 'PRES_ADJUSTED', IDVAR)
      IF (STAT /= NF90_NOERR) THEN
          DEPTHCA = DEPTHC
      ELSE
          STAT = NF90_GET_VAR (NCID,IDVAR,DEPTHCA)
          IF (STAT /= NF90_NOERR) CALL HANDLE_ERR(STAT)
          WRITE(IOUN,*) 'PRES_ADJUSTED READ'
      ENDIF
      ! PRESSURE TO DEPTH CONVERSION
      DO JJ=1,NPROFS
        X = SIN( (LATS(JJ)/57.29578_R8)*DEG2RAD )
        X = X*X
        DO JI=1,NLEVELS
          P = DEPTHC(JI,JJ)
          IF( P .LT. 99990._R8 ) THEN
            G = 9.780318 * ( 1._R8 + ( 5.2788E-3 + 2.36E-5  * X) * X ) + 1.092E-6 * P
            DEPTHC(JI,JJ) = ((((-1.82E-15 * P + 2.279E-10 ) * P - 2.2512E-5 ) * P + 9.72659) * P) / G
          ENDIF
          P = DEPTHCA(JI,JJ)
          IF( P .LT. 99990._R8 ) THEN
            G = 9.780318 * ( 1._R8 + ( 5.2788E-3 + 2.36E-5  * X) * X ) + 1.092E-6 * P
            DEPTHCA(JI,JJ) = ((((-1.82E-15 * P + 2.279E-10 ) * P - 2.2512E-5 ) * P + 9.72659) * P) / G
          ENDIF
        ENDDO
      ENDDO
    ENDIF

    IF(LSTRICTPSU .AND. .NOT. LXBT) THEN
          STAT = NF90_INQ_VARID (NCID, 'PSAL_ADJUSTED_QC', IDVAR)
          IF (STAT /= NF90_NOERR) THEN
              STAT = NF90_INQ_VARID (NCID, 'PSAL_QC', IDVAR)
          ENDIF
          IF (STAT /= NF90_NOERR) CALL HANDLE_ERR(STAT)
          STAT = NF90_GET_VAR (NCID,IDVAR,SQC(NPROFS,1:NUMLEVS))
          IF (STAT /= NF90_NOERR) CALL HANDLE_ERR(STAT)
          WRITE(IOUN,*) 'PSAL_ADJUSTED_QC READ'
    ENDIF

    IF(LSTRICTTEM) THEN
          STAT = NF90_INQ_VARID (NCID, 'TEMP_ADJUSTED_QC', IDVAR)
          IF (STAT /= NF90_NOERR) THEN
              STAT = NF90_INQ_VARID (NCID, 'TEMP_QC', IDVAR)
          ENDIF
          IF (STAT /= NF90_NOERR) CALL HANDLE_ERR(STAT)
          STAT = NF90_GET_VAR (NCID,IDVAR,TQC(NPROFS,1:NUMLEVS))

          IF (STAT /= NF90_NOERR) CALL HANDLE_ERR(STAT)
          WRITE(IOUN,*) 'POTM_ADJUSTED_QC READ'
    ENDIF

    STAT = NF90_CLOSE(NCID)


!========    OBS SELECTION    ========!

    KPROF=0

    CYCLE_PROF : DO JPROF = 1, NPROFS

       LLPROFOK = .FALSE.
       PLANO(JPROF)(1:8) = ADJUSTL(PLANO(JPROF)(1:8))

       !... CHECK TIME
       IF(LLCHECKTIME) THEN
           IF( TIME(JPROF) .LT. DATE_START .OR. &
                TIME(JPROF) .GT. DATE_END   ) THEN
               KERR(15) = KERR(15) + 1
               CYCLE CYCLE_PROF
           ENDIF
       ENDIF

       !... CHECK LATLON
       IF ( LATS(JPROF) .LT. -89._R8 .OR. LATS(JPROF) .GT. 90._R8 .OR. &
          & LONS(JPROF) .LT.-180._R8 .OR. LONS(JPROF) .GT.180._R8 ) THEN
           KERR(1) = KERR(1) + 1
           CYCLE CYCLE_PROF
       ENDIF
       IF ( LLSTRICTPOS ) THEN
           ! IF( ISNOGOOD( POSQC(JPROF:JPROF) ) ) THEN
           IF( ISNOGOOD( POSQC(JPROF) ) ) THEN
               KERR(4) = KERR(4) + 1
               CYCLE CYCLE_PROF
           ENDIF
       ENDIF

       !... CHECK DATE
       IF (TIME(JPROF) .LE. 0._DP  ) THEN
           KERR(2) = KERR(2) + 1
           CYCLE CYCLE_PROF
       ENDIF
       IF ( LLSTRICTTIME ) THEN
           ! IF( ISNOGOOD( TIMQC(JPROF:JPROF) ) ) THEN
           IF( ISNOGOOD( TIMQC(JPROF) ) ) THEN
               KERR(16) = KERR(16) + 1
               CYCLE CYCLE_PROF
           ENDIF
       ENDIF


!... A.S.: OBS TYPE DEFINITION IS IN ACCORDANCE WITH
!    UKMO SIMPLIFIED VERSION OF WMO 1770 TABLE

       READ(INST(JPROF)(1:4),*) ITYPE
!
! XBT file BA     , 996
! ARGO file PF    , 850
! GLIDER file GL  , 999
!
       SELECT CASE(ITYPE)
            CASE(1:510,800,810,825,900,995,996)   ! XBT/MBT
               ITYPE=KKXBT
               KCO(1) = KCO(1) + 1
            CASE(700:781,830)   ! TESAC
               ITYPE=KKTESAC
               KCO(2) = KCO(2) + 1
            CASE(831:861)   ! ARGO
               ITYPE=KKARGO
               KCO(3) = KCO(3) + 1
            CASE(820)   ! BUOYS
               ITYPE=KKBUOY
               KCO(4) = KCO(4) + 1
            CASE(999)   ! BUOYS
               SELECT CASE(CINST)
                 CASE('XB')
                  ITYPE=KKXBT
                  KCO(1) = KCO(1) + 1
                 CASE('BA')
                  ITYPE=KKXBT
                  KCO(1) = KCO(1) + 1
                 CASE('TE')
                  ITYPE=KKTESAC
                  KCO(2) = KCO(2) + 1
                 CASE('OC')
                  ITYPE=KKXBT
                  KCO(1) = KCO(1) + 1
                 CASE('PF')
                  ITYPE=KKARGO
                  KCO(3) = KCO(3) + 1
                 CASE('HF')
                  ITYPE=KKBUOY
                  KCO(4) = KCO(4) + 1
                 CASE('CT')
                  ITYPE=KKTESAC
                  KCO(2) = KCO(2) + 1
                 CASE('MO')
                  ITYPE=KKBUOY
                  KCO(4) = KCO(4) + 1
                 CASE('GL')
                  ITYPE=KKBUOY
                  KCO(4) = KCO(4) + 1
                 CASE DEFAULT
                  WRITE(IOUN,*) 'UNRECOGNIZED INSTRUMENTAL FILE'
                  CALL ABOR1('READCORFILE : UNRECOGNIZED INSTRUMENTAL FILE')
               ENDSELECT
            CASE DEFAULT
               WRITE(IOUN,*) &
               & 'READCORFILE : FOUND OBSTYPE = ',TRIM(INST(JPROF)),&
               & ' WHICH IS UNKNOWN!'
               KERR(3) = KERR(3) + 1
               CYCLE CYCLE_PROF
       ENDSELECT


       LPROFT=.TRUE.
       IF(LGLOBTQF) THEN
            ! IF( ISNOGOOD( TGQC(JPROF:JPROF) ) ) THEN
            IF( ISNOGOOD( TGQC(JPROF) ) ) THEN
               KERR(7) = KERR(7) +1
               LPROFT=.FALSE.
            ENDIF
       ENDIF

       LPROFS=.TRUE.
       IF(LGLOBSQF .AND. .NOT. LXBT) THEN
            ! IF( ISNOGOOD( SGQC(JPROF:JPROF) ) ) THEN
            IF( ISNOGOOD( SGQC(JPROF) ) ) THEN
               KERR(8) = KERR(8) +1
               LPROFS=.FALSE.
            ENDIF
       ENDIF
       IF( LXBT ) LPROFS=.FALSE.

       IF(.NOT.LPROFT .AND. .NOT.LPROFS) THEN
            KERR(5) = KERR(5) + 1
            CYCLE CYCLE_PROF
       ENDIF

       KOK(4) = KOK(4) + 1
       DPS_T = -999._R8
       DPS_S = -999._R8
       NTC   = 0
       ZTC   = 0._R8

       CYCLE_LEV : DO JLEV = 1,NLEVELS

          IF(DEPTHC(JLEV,JPROF) .GT. 99990._R8 ) THEN
                 KERR(6) = KERR(6) +1
                 CYCLE CYCLE_LEV
          ELSE
             IF(DEPTHCA(JLEV,JPROF) .LT. 99990._R8 ) DEPTHC(JLEV,JPROF) = DEPTHCA(JLEV,JPROF)
          ENDIF

          ! BETTER TO DO HERE
          IF( LL_DENSLEVS ) THEN
              ZDENS = SIGMA0(TEMP(JLEV,JPROF),SALC(JLEV,JPROF))
              IF( ZDENS .GT. 7._R8 .AND. ZDENS .LT. 18.5_R8 ) THEN
                 DEPTHC(JLEV,JPROF) = ZDENS
              ELSE
                 KERR(17) = KERR(17) +1
                 CYCLE CYCLE_LEV
              ENDIF
          ENDIF

        !... TEMPERATURE
          IF(LPROFT) THEN
            IF (TEMP(JLEV,JPROF).LE. 40._R8 .AND. TEMP(JLEV,JPROF) .GE. -2._R8 ) THEN
                  LOBS=.TRUE.
                  IF (LOBS .AND. TEMPA(JLEV,JPROF).LE. 40._R8 .AND. TEMPA(JLEV,JPROF) .GE. -2._R8 ) &
                  & TEMP(JLEV,JPROF) = TEMPA(JLEV,JPROF)
                  IF(LSTRICTTEM) THEN
                    ! IF( ISNOGOOD( TQC(JPROF)(JLEV:JLEV) ) ) LOBS=.FALSE.
                    IF( ISNOGOOD( TQC(JLEV,JPROF) ) ) LOBS=.FALSE.
                  ENDIF
                  IF( DPS_T .GT. 0._R8 .AND. DEPTHC(JLEV,JPROF)-DPS_T .LT. &
                  VTHIN_INT(DEPTHC(JLEV,JPROF) ) .AND. LVTHIN ) THEN
                     KERR(21) = KERR(21) +1
                     LOBS=.FALSE.
                  ENDIF
                  IF( LVAVG .AND. DPS_T .GT. 0._R8 ) THEN
                   IF( DEPTHC(JLEV,JPROF)-DPS_T .LT. VTHIN_INT(DEPTHC(JLEV,JPROF) ) ) THEN
                     NTC(2) = NTC(2)+1
                     ZTC(2) = ZTC(2)+TEMP(JLEV,JPROF)
                     LOBS=.FALSE.
                   ELSE
                     NTC(2) = NTC(2)+1
                     ZTC(2) = ZTC(2)+TEMP(JLEV,JPROF)
                     OVAL   = ZTC(2)/REAL(NTC(2),R8)
                     NTC(2) = 0
                     ZTC(2) = 0._R8
                   ENDIF
                  ELSE
                     OVAL = TEMP(JLEV,JPROF)
                  ENDIF
                  IF(LOBS) THEN
                    IF( .NOT. LLPROFOK ) THEN
                       LLPROFOK=.TRUE.
                       KPROF=KPROF+1
                    ENDIF
                    JOBS=JOBS+1
                    IF(JOBS.GT.MAXOBS) &
                    & CALL ABOR1('READCORFILE : MAXOBS MUST BE INCREASED!')
                    OBSTYPE(JOBS) = ITYPE
                    OBSDEP(JOBS) = DEPTHC(JLEV,JPROF)
                    OBSPAR(JOBS) = KKTEMP
                    OBSVAL(JOBS) = OVAL
                    OBSTIME(JOBS) = TIME(JPROF)
                    OBSLON(JOBS) = LONS(JPROF)
                    OBSLAT(JOBS) = LATS(JPROF)
                    OBSPROF(JOBS) = KPROF
                    OBSPLAT(JOBS) = PLANO(JPROF)
                    DPS_T = OBSDEP(JOBS)
                    IF( SALC(JLEV,JPROF) .LE. 40._R8 .AND. SALC(JLEV,JPROF) .GT.  0._R8 ) THEN
                      IF( SALCA(JLEV,JPROF) .LE. 40._R8 .AND. SALCA(JLEV,JPROF) .GT.  0._R8 ) THEN
                        OBSVAL(JOBS) = POTEMP( SALCA(JLEV,JPROF), OBSVAL(JOBS),&
                        & DEP_TO_P(DEPTHC(JLEV,JPROF),LATS(JPROF)), 0._R8 )
                      ELSE
                        OBSVAL(JOBS) = POTEMP( SALC(JLEV,JPROF), OBSVAL(JOBS),&
                        & DEP_TO_P(DEPTHC(JLEV,JPROF),LATS(JPROF)), 0._R8 )
                      ENDIF
                    ELSE
                      OBSVAL(JOBS) = OBSVAL(JOBS) + 200._R8
                    ENDIF
                    KOK(1) = KOK(1) + 1
                  ELSE
                    KERR(9) = KERR(9) +1
                  ENDIF
            ELSE
                  KERR(10) = KERR(10) +1
            ENDIF
          ENDIF

        !... SALINITY
          IF(LPROFS) THEN
            IF (SALC(JLEV,JPROF).LE. 42._R8 .AND. SALC(JLEV,JPROF) .GE. 0._R8 ) THEN
                  LOBS=.TRUE.
                  IF (LOBS .AND. SALCA(JLEV,JPROF).LE. 42._R8 .AND. SALCA(JLEV,JPROF) .GE. 0._R8 ) &
                  & SALC(JLEV,JPROF) = SALCA(JLEV,JPROF)
                  IF(LSTRICTPSU) THEN
                    ! IF( ISNOGOOD( SQC(JPROF)(JLEV:JLEV) ) ) LOBS=.FALSE.
                    IF( ISNOGOOD( SQC(JLEV,JPROF) ) ) LOBS=.FALSE.
                  ENDIF
                  IF( DPS_S .GE. 0._R8 .AND. DEPTHC(JLEV,JPROF)-DPS_S .LT. &
                  VTHIN_INT(DEPTHC(JLEV,JPROF) ) .AND. LVTHIN ) THEN
                     KERR(22) = KERR(22) +1
                     LOBS=.FALSE.
                  ENDIF
                  IF( LVAVG .AND. DPS_S .GT. 0._R8 ) THEN
                   IF( DEPTHC(JLEV,JPROF)-DPS_S .LT. VTHIN_INT(DEPTHC(JLEV,JPROF) ) ) THEN
                     NTC(1) = NTC(1)+1
                     ZTC(1) = ZTC(1)+SALC(JLEV,JPROF)
                     LOBS=.FALSE.
                   ELSE
                     NTC(1) = NTC(1)+1
                     ZTC(1) = ZTC(1)+SALC(JLEV,JPROF)
                     OVAL   = ZTC(1)/REAL(NTC(1),R8)
                     NTC(1) = 0
                     ZTC(1) = 0._R8
                   ENDIF
                  ELSE
                     OVAL = SALC(JLEV,JPROF)
                  ENDIF
                   IF(LOBS) THEN
                    IF( .NOT. LLPROFOK ) THEN
                       LLPROFOK=.TRUE.
                       KPROF=KPROF+1
                    ENDIF
                    JOBS=JOBS+1
                    IF(JOBS.GT.MAXOBS) &
                    & CALL ABOR1('READCORFILE : MAXOBS MUST BE INCREASED!')
                    OBSTYPE(JOBS) = ITYPE
                    OBSDEP(JOBS) = DEPTHC(JLEV,JPROF)
                    OBSPAR(JOBS) = KKSAL
                    OBSVAL(JOBS) = OVAL
                    OBSTIME(JOBS) = TIME(JPROF)
                    OBSLON(JOBS) = LONS(JPROF)
                    OBSLAT(JOBS) = LATS(JPROF)
                    OBSPROF(JOBS) = KPROF
                    OBSPLAT(JOBS) = PLANO(JPROF)
                    DPS_S = OBSDEP(JOBS)
                    KOK(2) = KOK(2) + 1
                  ELSE
                    KERR(11) = KERR(11) +1
                  ENDIF
            ELSE
                  KERR(12) = KERR(12) +1
            ENDIF
          ENDIF

       ENDDO CYCLE_LEV

    ENDDO CYCLE_PROF

  NTOTOBS=JOBS

  WRITE(IOUN,*)
  WRITE(IOUN,*) '*** REPORT OF READCORFILE : ',&
  & 'READING OF ENSEMBLES INS OBSERVATIONS'
  WRITE(IOUN,*)
  WRITE(IOUN,*) ' // OPTIONS'
  WRITE(IOUN,*) ' IS TIME OF PROFILE CHECKED ?        : ',LLCHECKTIME
  WRITE(IOUN,*) ' IS QUALITY OF POSITION CHECKED ?    : ',LLSTRICTPOS
  WRITE(IOUN,*) ' IS QUALITY OF TIME     CHECKED ?    : ',LLSTRICTTIME
  WRITE(IOUN,*) ' IS QUALITY OF SALINITY CHECKED ?    : ',LSTRICTPSU
  WRITE(IOUN,*) ' IS QUALITY OF TEMPERAT. USED FOR T ?: ',LSTRICTTEM
  WRITE(IOUN,*) ' IS GLOBAL QUALITY OF T PROF. USED ? : ',LGLOBTQF
  WRITE(IOUN,*) ' IS GLOBAL QUALITY OF S PROF. USED ? : ',LGLOBSQF
  WRITE(IOUN,*) ' IS VERTICAL THINNING FOR T/S USED ? : ',LVTHIN
  WRITE(IOUN,*)
  WRITE(IOUN,*) ' // DIMENSION'
  WRITE(IOUN,*) ' LEVELS      : ',NLEVELS
  WRITE(IOUN,*) ' PROFILES    : ',NPROFS
  WRITE(IOUN,*)
  WRITE(IOUN,*) ' // OBS PROFILES (WITHOUT UNVALID TIME/POSITION PROFILES)'
  WRITE(IOUN,*) ' XBT/MBT     : ',KCO(1)
  WRITE(IOUN,*) ' TESAC       : ',KCO(2)
  WRITE(IOUN,*) ' ARGO        : ',KCO(3)
  WRITE(IOUN,*) ' BUOYS       : ',KCO(4)
  WRITE(IOUN,*)
  WRITE(IOUN,*) ' // REPORT'
  WRITE(IOUN,*) ' TOTAL OBS STORED                :', NTOTOBS
  WRITE(IOUN,*) ' PROFILES WITH TIME OUT OF WINDOW:', KERR(15)
  WRITE(IOUN,*) ' PROFILES WITH POSITION PROBLEMS :', KERR(1)
  WRITE(IOUN,*) ' PROFILES WITH DATETIME PROBLEMS :', KERR(2)
  WRITE(IOUN,*) ' PROFILES WITH UNKNOWN INSTRUMENT:', KERR(3)
  WRITE(IOUN,*) ' POSITION QUALITY CHECK NEGATIVE :', KERR(4)
  WRITE(IOUN,*) ' TIME     QUALITY CHECK NEGATIVE :', KERR(16)
  WRITE(IOUN,*) ' PROFILE QC NEGATIVE FOR TEMPER. :', KERR(7)
  WRITE(IOUN,*) ' PROFILE QC NEGATIVE FOR SALINITY:', KERR(8)
  WRITE(IOUN,*) ' PROFILE QC NEGATIVE FOR T AND S :', KERR(5)
  WRITE(IOUN,*) ' NOT VALID DEPTHS                :', KERR(6)

  WRITE(IOUN,*) ' OBSERVATION QC NEGATIVE FOR T   :', KERR(9)
  WRITE(IOUN,*) ' NOT VALID VALUE FOR T           :', KERR(10)
  WRITE(IOUN,*) ' OBSERVATION QC NEGATIVE FOR S   :', KERR(11)
  WRITE(IOUN,*) ' NOT VALID VALUE FOR S           :', KERR(12)
  WRITE(IOUN,*) ' OBSERVATION QC NEGATIVE FOR PT  :', KERR(13)
  WRITE(IOUN,*) ' NOT VALID VALUE FOR PT          :', KERR(14)
  WRITE(IOUN,*) ' DENSITY LEVEL NOT AVAILABLE     :', KERR(17)
  WRITE(IOUN,*) ' VERTICAL THINNING REJECTION - T :', KERR(21)
  WRITE(IOUN,*) ' VERTICAL THINNING REJECTION - S :', KERR(22)

  WRITE(IOUN,*) ' VALID PROFILES                  :', KOK (4)
  WRITE(IOUN,*) ' EXTRACTED OBS OF SALINITY       :', KOK (2)
  WRITE(IOUN,*) ' EXTRACTED OBS OF TEMPERATURE    :', KOK (1)
  WRITE(IOUN,*) ' EXTRACTED OBS OF POT. TEMPER.   :', KOK (3)
  WRITE(IOUN,*)
  WRITE(IOUN,*) '*** END OF READCORFILE'
  WRITE(IOUN,*)

    DEALLOCATE(LATS,LONS,TIME,INST,PLANO)
    DEALLOCATE(TEMP,SALC,DEPTHC)
    DEALLOCATE(TEMPA,SALCA,DEPTHCA)

CONTAINS

LOGICAL FUNCTION ISNOGOOD( CHIN )
IMPLICIT NONE
!!! CHARACTER(LEN=1) :: CHIN
INTEGER(4), INTENT(IN) :: CHIN
ISNOGOOD = .FALSE.
IF( CHIN .EQ. 3 .OR. CHIN .EQ. 4 ) ISNOGOOD = .TRUE.
END FUNCTION ISNOGOOD

   REAL(KIND=R8) FUNCTION dep_to_p( p_dep, p_lat )
    IMPLICIT NONE
      REAL(KIND=R8), INTENT(IN) :: p_dep    ! Depth in meters
      REAL(KIND=R8), INTENT(IN) :: p_lat    ! Latitude in degrees
      REAL(KIND=R8) :: z_x
      REAL(KIND=R8) :: z_c1
      REAL(KIND=R8) :: z_c2
      REAL(KIND=R8) :: z_d
 
      z_x = SIN( p_lat / 57.29578_r8 )
      z_x = z_x * z_x
      z_c1 = ( 5.92_r8  + 5.25_r8 * z_x ) * 1e-3_r8
      z_c2 = 2.21e-6_r8
      z_d = ( z_c1 - 1._r8 ) * ( z_c1 - 1._r8  ) - 4._r8 * z_c2 * p_dep
      dep_to_p = (( 1._r8 - z_c1 ) - SQRT( z_d )) / ( 2._r8 * z_c2 )

   END FUNCTION dep_to_p

   REAL(KIND=R8) FUNCTION potemp( ps, pt, pp, ppr )

    IMPLICIT NONE
      REAL(KIND=R8), INTENT(IN) :: ps
      REAL(KIND=R8), INTENT(IN) :: pt
      REAL(KIND=R8), INTENT(IN) :: pp
      REAL(KIND=R8), INTENT(IN) :: ppr

      REAL(KIND=R8) :: zpol
      REAL(KIND=R8), PARAMETER :: a1 =  1.067610e-05_r8
      REAL(KIND=R8), PARAMETER :: a2 = -1.434297e-06_r8
      REAL(KIND=R8), PARAMETER :: a3 = -7.566349e-09_r8
      REAL(KIND=R8), PARAMETER :: a4 = -8.535585e-06_r8
      REAL(KIND=R8), PARAMETER :: a5 =  3.074672e-08_r8
      REAL(KIND=R8), PARAMETER :: a6 =  1.918639e-08_r8
      REAL(KIND=R8), PARAMETER :: a7 =  1.788718e-10_r8

        zpol = a1 + a2 * ps + a3 * ( pp + ppr ) + a4 * pt &
         & + a5 * ps * pt + a6 * pt * pt + a7 * pt * ( pp + ppr )

      potemp = pt + ( pp - ppr ) * zpol
   END FUNCTION potemp

  FUNCTION sigma0 ( ptem, psal)
    IMPLICIT NONE
    !!----------------------------------------------------------------------
    REAL(KIND=R8), INTENT(in) :: ptem, psal 
    REAL(KIND=R8)  :: sigma0      ! returned value

    REAL(KIND=R8)                      :: zt, zs, zsr, zrau0=1000._R8, zws
    REAL(KIND=R8)                      :: zr1, zr2, zr3, zr4
    !!----------------------------------------------------------------------
    zws = 0._R8
    sigma0 = 0._R8
    zws = SQRT( ABS(psal) )
    zt  = ptem
    zs  = psal
    zsr = zws
    ! compute volumic mass pure water at atm pressure
    zr1 = ( ( ( ( 6.536332e-9*zt-1.120083e-6 )*zt+1.001685e-4)*zt   &
           -9.095290e-3 )*zt+6.793952e-2 )*zt+999.842594
    ! seawater volumic mass atm pressure
    zr2= ( ( ( 5.3875e-9*zt-8.2467e-7 )*zt+7.6438e-5 ) *zt   &
           -4.0899e-3 ) *zt+0.824493
    zr3= ( -1.6546e-6*zt+1.0227e-4 ) *zt-5.72466e-3
    zr4= 4.8314e-4
    ! potential volumic mass (reference to the surface)
    sigma0 = ( zr4*zs + zr3*zsr + zr2 ) *zs + zr1 - zrau0
  END FUNCTION sigma0
  FUNCTION VTHIN_INT(DEP) 
    USE GRD_STR
    IMPLICIT NONE
    REAL(KIND=R8), INTENT(in) :: DEP
    REAL(KIND=R8)  :: VTHIN_INT
    INTEGER(I4) :: KKZ(1)
  IF( NN_VERTMESH .EQ. 1 ) THEN
    IF( DEP .LE. 100._R8 ) VTHIN_INT = 1._R8
    IF( DEP .LE. 1500._R8 .AND. DEP .GT. 100._R8 ) VTHIN_INT = 10._R8
    IF( DEP .LE. 6000._R8 .AND. DEP .GT. 1500._R8) VTHIN_INT = 50._R8
  ELSEIF ( NN_VERTMESH .EQ. 2 ) THEN
    KKZ = MINLOC( ABS(DEP-GRD%DEP) )
    VTHIN_INT = GRD%DZ(KKZ(1))
  ELSE
    CALL ABOR1('READCORFILE : UNRECOGNIZED OPTION FOR NN_VERTMESH')
  ENDIF
  END FUNCTION VTHIN_INT

END SUBROUTINE READCORFILE
