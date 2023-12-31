SUBROUTINE READSSS

!! READ SPACE-BORNE MW OR IR SSS FROM FILE

USE IOUNITS, ONLY : IOUNERR,IOUNOUT,IOUNLOG
USE RUN, ONLY  : CDATETIME,IANTW,IANDATE,ZANJUL1950,LL_TWCENTERED
USE SET_KND
USE OBS_STR
USE GRD_STR
USE ALLOCOBS
USE CALENDAR, ONLY : JU2YMDS

IMPLICIT NONE

INTEGER(I4) :: NOBS,IERR,JSAT,KFILE,KSAT,ISTART,&
             & K,IJ,JOBS,IDIFFER1,IRET,ISAT,KOS
LOGICAL :: LLEXIST,LEX
INTEGER(I4),PARAMETER :: MAX_SSSOBS = 10368000
CHARACTER(LEN=10) :: CSSSDATE,CTEMPDATE
CHARACTER(LEN=60) :: CFILEIN
REAL(R8),DIMENSION(MAX_SSSOBS) :: ZLON,ZLAT,ZDAT
REAL(R8),DIMENSION(MAX_SSSOBS,3) :: ZBIA
REAL(DP),DIMENSION(MAX_SSSOBS) :: ZTIM
INTEGER(I4) :: YYYY,MM,DD
REAL(DP) :: SS, RDAY
REAL(DP)          :: DTSTART,DTEND
CHARACTER(LEN=10) :: CDTSTART,CDTEND

WRITE(IOUNLOG,*)
WRITE(IOUNLOG,*) '... ENTERING READSSS ...'

DTSTART = ZJULSTART - ZJUL1950
DTEND   = ZJULEND   - ZJUL1950

WRITE(IOUNLOG,*) ' STARTDATE (JULIAN DAYS SINCE 1950) : ', DTSTART
WRITE(IOUNLOG,*) ' END  DATE (JULIAN DAYS SINCE 1950) : ', DTEND
WRITE(IOUNLOG,*)

SSS%NO=0

! TO DO : BIAS PREDICTORS

CALL ALLOCSSS(MAX_SSSOBS,.TRUE.)

   IF ( SSS%LLSACD ) THEN

       KOS=0
       ISAT = 3

       CSSSDATE=CDTSTART
       CALL JU2YMDS (ZJULSTART,YYYY,MM,DD,SS)
       WRITE(CSSSDATE,'(I4.4,I2.2,I2.2)') YYYY,MM,DD
       CALL JU2YMDS (ZJULEND,YYYY,MM,DD,SS)
       WRITE(CDTEND,'(I4.4,I2.2,I2.2)') YYYY,MM,DD
       RDAY=0._R8

       DO WHILE( CSSSDATE .NE. CDTEND )

          WRITE(CFILEIN,'(A,A,A)') 'SACD_',CSSSDATE(1:8),'.NC'
          INQUIRE(FILE=CFILEIN,EXIST=LEX)

          IF(LEX) THEN

           CALL READ_SACD(CFILEIN,MAX_SSSOBS,NOBS,ZDAT,&
           & ZLON,ZLAT,ZBIA)

           SSS%NO = SSS%NO + NOBS
           ISTART = SSS%NO-NOBS+1
           KOS = KOS + NOBS

           IF( SSS%NO .GT. MAX_SSSOBS ) &
           & CALL ABOR1('FATAL IN READSSS: MAX_SSSOBS MUST BE INCREASED')

           SSS%LON(ISTART:SSS%NO) = ZLON(1:NOBS)
           SSS%LAT(ISTART:SSS%NO) = ZLAT(1:NOBS)
           SSS%TIM(ISTART:SSS%NO) = DTSTART+RDAY
           SSS%VAL(ISTART:SSS%NO) = ZDAT(1:NOBS)
           SSS%BCP(ISTART:SSS%NO,1:3) = ZBIA(1:NOBS,1:3)

           SSS%TDIST(ISTART:SSS%NO) = SSS%TIM(ISTART:SSS%NO) - ZANJUL1950
           SSS%KSAT (ISTART:SSS%NO) = ISAT
           WRITE(IOUNLOG,*) CFILEIN,YYYY,MM,DD,&
           & SSS%TIM(ISTART),SSS%TIM(SSS%NO),RDAY

          ENDIF

          RDAY = RDAY + 1._R8
          CALL JU2YMDS (ZJULSTART+RDAY,YYYY,MM,DD,SS)
          WRITE(CSSSDATE,'(I4.4,I2.2,I2.2)') YYYY,MM,DD

       ENDDO

       SSS%ISSSOBS(ISAT) = KOS

   ENDIF

   IF ( SSS%LLAQUARIUS ) THEN

       KOS=0
       ISAT = 1

       CSSSDATE=CDTSTART
       CALL JU2YMDS (ZJULSTART,YYYY,MM,DD,SS)
       WRITE(CSSSDATE,'(I4.4,I2.2,I2.2)') YYYY,MM,DD
       CALL JU2YMDS (ZJULEND,YYYY,MM,DD,SS)
       WRITE(CDTEND,'(I4.4,I2.2,I2.2)') YYYY,MM,DD
       RDAY=0.5_R8

       DO WHILE( CSSSDATE .NE. CDTEND )

          WRITE(CFILEIN,'(A,A,A)') 'AQUARIUS_SSS_',CSSSDATE(1:8),'.NC'
          INQUIRE(FILE=CFILEIN,EXIST=LEX)

          IF(LEX) THEN

           CALL READ_AQUARIUS(CFILEIN,MAX_SSSOBS,NOBS,ZDAT,&
           & ZTIM,ZLON,ZLAT)

           SSS%NO = SSS%NO + NOBS
           ISTART = SSS%NO-NOBS+1
           KOS = KOS + NOBS

           IF( SSS%NO .GT. MAX_SSSOBS ) &
           & CALL ABOR1('FATAL IN READSSS: MAX_SSSOBS MUST BE INCREASED')

           SSS%LON(ISTART:SSS%NO) = ZLON(1:NOBS)
           SSS%LAT(ISTART:SSS%NO) = ZLAT(1:NOBS)
           SSS%TIM(ISTART:SSS%NO) = ZTIM(1:NOBS)
           SSS%VAL(ISTART:SSS%NO) = ZDAT(1:NOBS)

           SSS%TDIST(ISTART:SSS%NO) = SSS%TIM(ISTART:SSS%NO) - ZANJUL1950
           SSS%KSAT (ISTART:SSS%NO) = ISAT

          ENDIF

          RDAY = RDAY + 1._R8
          CALL JU2YMDS (ZJULSTART+RDAY,YYYY,MM,DD,SS)
          WRITE(CSSSDATE,'(I4.4,I2.2,I2.2)') YYYY,MM,DD

       ENDDO

       SSS%ISSSOBS(ISAT) = KOS

   ENDIF

   IF ( SSS%LLSMOS ) THEN

       CALL ABOR1('READSSS: SMOS NOT YET SUPPORTED')

       KOS=0
       ISAT = 2

       CSSSDATE=CDTSTART
       CALL JU2YMDS (ZJULSTART,YYYY,MM,DD,SS)
       WRITE(CSSSDATE,'(I4.4,I2.2,I2.2)') YYYY,MM,DD
       CALL JU2YMDS (ZJULEND,YYYY,MM,DD,SS)
       WRITE(CDTEND,'(I4.4,I2.2,I2.2)') YYYY,MM,DD
       RDAY=0.5_R8

       DO WHILE( CSSSDATE .NE. CDTEND )

          WRITE(CFILEIN,'(A,A,A)') 'AMSRE_',CSSSDATE(1:8),'V5'
          INQUIRE(FILE=CFILEIN,EXIST=LEX)

          IF(LEX) THEN

           CALL DIFFTIME('D',CSSSDATE(1:8),'19500101',IDIFFER1,IRET,.FALSE.)

           CALL READ_AMSR_DAY(CFILEIN,IDIFFER1,MAX_SSSOBS,NOBS,ZDAT,&
           & ZTIM,ZLON,ZLAT)

           SSS%NO = SSS%NO + NOBS
           ISTART = SSS%NO-NOBS+1
           KOS = KOS + NOBS

           IF( SSS%NO .GT. MAX_SSSOBS ) &
           & CALL ABOR1('FATAL IN READSSS: MAX_SSSOBS MUST BE INCREASED')

           SSS%LON(ISTART:SSS%NO) = ZLON(1:NOBS)
           SSS%LAT(ISTART:SSS%NO) = ZLAT(1:NOBS)
           SSS%TIM(ISTART:SSS%NO) = ZTIM(1:NOBS)
           SSS%VAL(ISTART:SSS%NO) = ZDAT(1:NOBS)

           SSS%TDIST(ISTART:SSS%NO) = SSS%TIM(ISTART:SSS%NO) - ZANJUL1950
           SSS%KSAT (ISTART:SSS%NO) = ISAT

          ENDIF

          RDAY = RDAY + 1._R8
          CALL JU2YMDS (ZJULSTART+RDAY,YYYY,MM,DD,SS)
          WRITE(CSSSDATE,'(I4.4,I2.2,I2.2)') YYYY,MM,DD

       ENDDO

       SSS%ISSSOBS(ISAT) = KOS

   ENDIF

WRITE(IOUNLOG,*)
WRITE(IOUNLOG,*) ' ========================================'
WRITE(IOUNLOG,*) ' SSS REPORT : TOTAL NO OF OBS IS',SSS%NO

END SUBROUTINE READSSS
