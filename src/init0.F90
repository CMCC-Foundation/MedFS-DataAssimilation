SUBROUTINE INIT0

!.. READ IN COMMAND-LINE AND ENVIRONMENTAL VARIABLES,
!.. INIT TIME VARS FOR JULIAN CALENDAR COMPUTATIONS
!..
!.. ANDREA STORTO, 2009-03-20


USE SET_KND
USE IOUNITS, ONLY : IOUNLOG,IOUNERR
USE RUN
USE CALENDAR

IMPLICIT NONE

INTEGER(KIND=I4) :: IERR
CHARACTER(LEN=10)  :: CDATEX,CTIMEX,CZONE,CARG,CARG2
INTEGER(KIND=I4) :: ITVALUES(8),JARG
INTEGER,EXTERNAL :: IARGC
INTEGER(KIND=I4) :: IYY,IMM,IDD
REAL(DP) :: ZSS,ZINC1,ZINC2
! NB: RS6K_NO WAS RS6K BUT IT IS NOW TAKEN OUT
!     AS AIX CLOCK_ CONFLICTS WITH DRHOOK.

!... SAVE INITIAL TIME
CALL DATE_AND_TIME(CDATEX,CTIMEX,CZONE,ITVALUES)
WRITE(CINSTRING,'(A)') ' PROCESS STARTED AT '//&
               & CDATEX(1:8)//' '//CTIMEX(1:2)//':'//&
               & CTIMEX(3:4)//':'//CTIMEX(5:6)

IANDATE=999
IANTIME=999
IANTW=999
LL_TWCENTERED = .TRUE.

DO JARG=1,IARGC()

   CALL GETARG(JARG,CARG)
   SELECT CASE(CARG(1:2))

      CASE('-d')
          CALL GETARG(JARG+1,CARG2)
          READ(CARG2,*) IANDATE
      CASE('-t')
          CALL GETARG(JARG+1,CARG2)
          READ(CARG2,*) IANTIME
      CASE('-w')
          CALL GETARG(JARG+1,CARG2)
          READ(CARG2,*) IANTW
      CASE('-c')
          CALL GETARG(JARG+1,CARG2)
          IF( CARG2(1:1) .EQ. '0' ) LL_TWCENTERED = .FALSE.
      CASE DEFAULT
          CALL GETARG(JARG-1,CARG2)
          IF(CARG2(1:2) .NE. '-d' .AND. &
           & CARG2(1:2) .NE. '-t' .AND. &
           & CARG2(1:2) .NE. '-c' .AND. &
           & CARG2(1:2) .NE. '-w' ) &
           & CALL USAGEOV
   ENDSELECT

ENDDO

IF(IANDATE.EQ.999 .OR. IANTIME.EQ.999 .OR. IANTW.EQ.999) &
  & CALL USAGEOV

!... ROUGHLY CHECK DATES
IF(IANDATE.LT.19000000 .OR. IANDATE .GT.21200000 ) CALL USAGEOV
IF(IANTIME.LT.0000 .OR. IANTIME .GT.2359 ) CALL USAGEOV
IF(IANTW  .LT.0    .OR. IANTW   .GT.999  ) CALL USAGEOV

WRITE(CDATETIME,'(I8.8,I4.4)') IANDATE,IANTIME
WRITE(CDATE    ,'(I8.8     )') IANDATE
WRITE(    CTIME,'(     I4.4)')         IANTIME

IANYY = INT(IANDATE/10000)
IANMM = INT((IANDATE - IANYY*10000)/100)
IANDD = IANDATE - IANYY*10000 -IANMM*100
IANHH = INT(IANTIME/100)
IANMI = IANTIME - IANHH*100
!... THE FOLLOWING CODE GENERATES AN 'ILLEGAL OPCODE' ON PWR6/XLF11
#ifndef RS6K
ZANSS = (REAL(IANHH,KIND=DP)*3600._DP + REAL(IANMI,KIND=R8)*60._DP)
#else
! DO NOT CARE OF MINUTS FOR THE TIME BEING
IF(IANHH.EQ.0) THEN
   ZANSS=0._DP
ELSEIF(IANHH.EQ.6) THEN
   ZANSS=21600._DP
ELSEIF(IANHH.EQ.12) THEN
   ZANSS=43200._DP
ELSEIF(IANHH.EQ.18) THEN
   ZANSS=64800._DP
ELSE
   WRITE(IOUNERR,*) 'INIT0 (RS6K SPECIFIC) CANNOT PROCEED'
   CALL ABOR1('INIT0: UNABLE TO DETERMINE SECONDS-OF-DAY')
ENDIF
#endif

CALL IOCONF_CALENDAR('GREGORIAN')

CALL YMDS2JU(IANYY,IANMM,IANDD,ZANSS,ZJULIAN)
CALL YMDS2JU(1950,01,01,0._DP,ZJUL1950)
CALL YMDS2JU(IANYY,IANMM,1,0._DP,ZJULSTARTMM)

ZANJUL1950 = ZJULIAN - ZJUL1950

!... FIND START AND END DATES
IF( LL_TWCENTERED ) THEN
  ZINC1 = REAL(-1800*IANTW)
  ZINC2 = REAL(1800*IANTW)
ELSE
  ZINC1 = REAL(-3600*IANTW)
  ZINC2 = 0._DP
ENDIF

CALL TIME_ADD(IANYY,IANMM,IANDD,ZANSS,ZINC1,&
            & IYY,IMM,IDD,ZSS)
CALL YMDS2JU(IYY,IMM,IDD,ZSS,ZJULSTART)
CALL TIME_ADD(IANYY,IANMM,IANDD,ZANSS,ZINC2,&
            & IYY,IMM,IDD,ZSS)
CALL YMDS2JU(IYY,IMM,IDD,ZSS,ZJULEND  )


CARGS= ' '
DO JARG=0,IARGC()
   CALL GETARG(JARG,CARG)
   CARGS=TRIM(CARGS)//' '//TRIM(CARG)
ENDDO

CALL GETENV('GRIDS_PATH',CGRIDPATH)

IF (LEN_TRIM(CGRIDPATH) .EQ. 0) THEN
    CGRIDPATH='./'
ENDIF

RETURN
END SUBROUTINE INIT0
