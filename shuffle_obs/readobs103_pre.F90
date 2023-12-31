PROGRAM  READOBS103

!.. READ OBS IN BINARY OR NETCDF FORMAT

USE SET_KND
USE OBSDEF
USE OBS_STR
USE ALLOCOBS
USE ALLOCOBS2

IMPLICIT NONE

INTEGER(I4) :: NPROCS, JP
INTEGER(I4),ALLOCATABLE :: KKCNT(:,:)
INTEGER(I4) :: IAUX,IST,IEN,NMAXOBS,JO,NOINT,KOBS,JSAT,K
INTEGER(I4) :: NCID, DIMID, VARID, JOBS, MAXINS, MAXSLA, MAXSST, MAXSSS
CHARACTER(LEN=30) :: CFILE, CAUX
CHARACTER(LEN=4) :: CSTR
LOGICAL :: LL_AB
TYPE (INS_T) :: INSB
TYPE (SLA_T) :: SLAB
TYPE (SST_T) :: SSTB
TYPE (SSS_T) :: SSSB

INTEGER, PARAMETER :: IO_OBSFORMAT = 1
INTEGER :: IARGC

#include "io_obs.h"

IF(IARGC().NE.1 .AND. IARGC().NE.2) CALL ABOR1('Usage: shuffle_obs.x <number_of_domains> [no_abort]')
CALL GETARG(1,CAUX)
READ(CAUX,*) NSUBDOMAINS

LL_AB = .TRUE.
IF( IARGC() .EQ. 2 ) THEN
   CALL GETARG(2,CAUX)
   IF( CAUX(1:1) .EQ. 'T' ) LL_AB=.FALSE.
ENDIF

LLPREP = .TRUE.

WRITE(IOUNOUT,*) ' NUMBER OF DOMAINS IS ',NSUBDOMAINS

CALL SET_DIM

NPROCS = NSUBDOMAINS

CALL SETFGFLAYOUT

!.. GET OBS VECT INDEX

ALLOCATE(KKCNT(NPROCS,NOFAMS))

OPEN(913,FILE='OBS_TAB.DAT',STATUS='OLD')
DO JP=1,NPROCS
    READ(913,*) IAUX,KKCNT(JP,NNINS),KKCNT(JP,NNSLA),&
                   & KKCNT(JP,NNSST),KKCNT(JP,NNSSS)
    WRITE(IOUNOUT,*) '/// ', IAUX,KKCNT(JP,NNINS),KKCNT(JP,NNSLA),&
                   & KKCNT(JP,NNSST),KKCNT(JP,NNSSS)
ENDDO
CLOSE(913)

INS%NO = SUM(KKCNT(:,NNINS))
INS%NC = INS%NO

IF(INS%NO.GT.0) THEN

  WRITE(IOUNLOG,*) ' ** ALLOCATING INSITU OBS VECTOR **'
  WRITE(IOUNLOG,*) '    NUMBER OF OBS IS  ',INS%NO
  CALL ALLOCINS(INS%NO,.FALSE.)
  MAXINS = MAXVAL(KKCNT(:,NNINS))
  WRITE(IOUNLOG,*) ' ** ALLOCATING AUXILIARY INSITU OBS VECTOR **'
  WRITE(IOUNLOG,*) '    NUMBER OF OBS IS  ',MAXINS
  CALL ALLOCINS2(INS%NO,.FALSE.,INSB)

  INS%FLC = 0
  INS%FLG = 0

ENDIF

SLA%NO=SUM(KKCNT(:,NNSLA))
SLA%NC=SLA%NO

IF(SLA%NO.GT.0) THEN
   
  WRITE(IOUNLOG,*) ' ** ALLOCATING SLA OBS VECTOR **'
  WRITE(IOUNLOG,*) '    NUMBER OF OBS IS  ',SLA%NO
  CALL ALLOCSLA(SLA%NO,.TRUE.)
  MAXSLA = MAXVAL(KKCNT(:,NNSLA))
  WRITE(IOUNLOG,*) ' ** ALLOCATING AUXILIARY SLA OBS VECTOR **'
  WRITE(IOUNLOG,*) '    NUMBER OF OBS IS  ',MAXSLA
  CALL ALLOCSLA2(SLA%NO,.FALSE.,SLAB)

  SLA%FLC = 0
  SLA%FLG = 0

ENDIF

SST%NO=SUM(KKCNT(:,NNSST))
SST%NC=SST%NO

IF(SST%NO.GT.0) THEN
  CALL ALLOCSST(SST%NO,.TRUE.)
  MAXSST = MAXVAL(KKCNT(:,NNSST))
  WRITE(IOUNLOG,*) ' ** ALLOCATING AUXILIARY SSTITU OBS VECTOR **'
  WRITE(IOUNLOG,*) '    NUMBER OF OBS IS  ',MAXSST
  CALL ALLOCSST2(SST%NO,.FALSE.,SSTB)
  SST%FLC = 0
  SST%FLG = 0
ENDIF

SSS%NO=SUM(KKCNT(:,NNSSS))
SSS%NC=SSS%NO

IF(SSS%NO.GT.0) THEN
  CALL ALLOCSSS(SSS%NO,.TRUE.)
  MAXSSS = MAXVAL(KKCNT(:,NNSSS))
  WRITE(IOUNLOG,*) ' ** ALLOCATING AUXILIARY SSSITU OBS VECTOR **'
  WRITE(IOUNLOG,*) '    NUMBER OF OBS IS  ',MAXSSS
  CALL ALLOCSSS2(SSS%NO,.FALSE.,SSSB)
  SSS%FLC = 0
  SSS%FLG = 0
ENDIF

CALL FLUSH(IOUNLOG)

DO JP=1,NPROCS

   IST=1
   IEN=KKCNT(JP,NNINS)

   WRITE(IOUNLOG,*) ' JP, KKCNT : ',JP,IEN
   CALL FLUSH(IOUNLOG)

   IF(KKCNT(JP,NNINS) .GT. 0 ) THEN

      WRITE(CFILE,'(A,I4.4,A)') 'INSOBS_',JP-1,'.NC'

      WRITE(IOUNLOG,*) ' FETCHING OBS IN FILE ',TRIM(CFILE)
      CALL FLUSH(IOUNLOG)

      CALL IO_OBS(CFILE=CFILE,IOPT=2,NOBS=KKCNT(JP,NNINS),&
      NLEVS=KM, NPREDS=1,NPQ=NPQ,NPQ2=2*NPQ,&
      FLC=INSB%FLC(IST:IEN),&
      NIND=INSB%NIND(IST:IEN),&
      INO=INSB%INO(IST:IEN),&
      OTYPE=INSB%OTYPE(IST:IEN),&
      PAR=INSB%PAR(IST:IEN),&
      PLNO=INSB%PLNO(IST:IEN),&
      INST=INSB%INST(IST:IEN),&
      EVE=INSB%EVE(IST:IEN),&
      KTY=INSB%KTY(IST:IEN),&
      PRIND=INSB%PRIND(IST:IEN),&
      PROF=INSB%PROF(IST:IEN),&
      TDIST=INSB%TDIST(IST:IEN),&
      LON=INSB%LON(IST:IEN),&
      LAT=INSB%LAT(IST:IEN),&
      DPT=INSB%DPT(IST:IEN),&
      TIM=INSB%TIM(IST:IEN),&
      VAL=INSB%VAL(IST:IEN),&
      BIA=INSB%BIA(IST:IEN),&
      IB=INSB%IB(IST:IEN,:),&
      JB=INSB%JB(IST:IEN,:),&
      KB=INSB%KB(IST:IEN),&
      RB=INSB%RB(IST:IEN),&
      PQ=INSB%PQ(IST:IEN,:))

     CYCINS : DO JOBS=1,KKCNT(JP,NNINS)
        KOBS = INSB%NIND(JOBS)
        IF( KOBS .EQ. 38798 ) WRITE(*,*) 'WARNING',&
            & JP,JOBS,KOBS,INS%NO,INS%FLG( KOBS )
        IF( KOBS .LT. 1 .OR. KOBS .GT. INS%NO .OR. INS%FLG( KOBS ) .EQ. 1 ) THEN
            IF( KOBS .LT. 1 .OR. KOBS .GT. INS%NO ) CSTR='OUTR'
            IF( INS%FLG( KOBS ) .EQ. 1            ) CSTR='DUPL'
            WRITE(332,*) 'ERROR : PROC, JOBS, KOBS, INS%NO, FLG',&
            & JP,JOBS,KOBS,INS%NO,INS%FLG( KOBS ),CSTR
            IF( LL_AB ) CALL ABOR1('INS ERROR')
        ENDIF 

        INS%FLG( KOBS )= 1
        INS%FLC( KOBS )=INSB%FLC(JOBS)
        INS%INO( KOBS )=INSB%INO(JOBS)
        INS%OTYPE( KOBS )=INSB%OTYPE(JOBS)
        INS%PAR( KOBS )=INSB%PAR(JOBS)
        INS%PLNO( KOBS )=INSB%PLNO(JOBS)
        INS%INST( KOBS )=INSB%INST(JOBS)
        INS%EVE( KOBS )=INSB%EVE(JOBS)
        INS%KTY( KOBS )=INSB%KTY(JOBS)
        INS%PROF( KOBS )=INSB%PROF(JOBS)
        INS%TDIST( KOBS )=INSB%TDIST(JOBS)
        INS%LON( KOBS )=INSB%LON(JOBS)
        INS%LAT( KOBS )=INSB%LAT(JOBS)
        INS%DPT( KOBS )=INSB%DPT(JOBS)
        INS%TIM( KOBS )=INSB%TIM(JOBS)
        INS%VAL( KOBS )=INSB%VAL(JOBS)
        INS%BIA( KOBS )=INSB%BIA(JOBS)
        INS%IB( KOBS,: )=INSB%IB(JOBS,:)
        INS%JB( KOBS,: )=INSB%JB(JOBS,:)
        INS%KB( KOBS )=INSB%KB(JOBS)
        INS%RB( KOBS )=INSB%RB(JOBS)
        INS%PQ( KOBS,: )=INSB%PQ(JOBS,:)
     ENDDO CYCINS

   ENDIF

   CALL FLUSH(IOUNLOG)

   IST=1
   IEN=KKCNT(JP,NNSLA)

   IF(KKCNT(JP,NNSLA) .GT. 0 ) THEN

      WRITE(CFILE,'(A,I4.4,A)') 'SLAOBS_',JP-1,'.NC'

      WRITE(IOUNLOG,*) ' SLA CHECK FOR PROCESSOR ',JP
      WRITE(IOUNLOG,*) ' ',ALLOCATED(SLA%TB)
      IF( ALLOCATED(SLA%TB) ) WRITE(IOUNLOG,*) ' DIMS:',SIZE(SLA%TB,1),SIZE(SLA%TB,2)
      WRITE(IOUNLOG,*) ' IST, IEN = ',IST,IEN,IEN-IST+1
      CALL FLUSH(IOUNLOG)

      CALL IO_OBS(CFILE=CFILE,IOPT=2,&
      & NOBS=KKCNT(JP,NNSLA),NLEVS=KM,NPREDS=NPRD_SLA,NPQ=NPQ,NPQ2=NPQ, &
      FLC=SLAB%FLC(IST:IEN),&
      NIND=SLAB%NIND(IST:IEN),&
      INO=SLAB%INO(IST:IEN),&
      TRACK=SLAB%TRACK(IST:IEN),&
      EVE=SLAB%EVE(IST:IEN),&
      KSAT=SLAB%KSAT(IST:IEN),&
      BOT=SLAB%BOT(IST:IEN),&
      TDIST=SLAB%TDIST(IST:IEN),&
      LON=SLAB%LON(IST:IEN),&
      LAT=SLAB%LAT(IST:IEN),&
      TIM=SLAB%TIM(IST:IEN),&
      VAL=SLAB%VAL(IST:IEN),&
      BIA=SLAB%BIA(IST:IEN),&
      IB=SLAB%IB(IST:IEN,:),&
      JB=SLAB%JB(IST:IEN,:),&
      PQ=SLAB%PQ(IST:IEN,:),&
      DPT=SLAB%DPT(IST:IEN))

      CYCSLA : DO JOBS=1,KKCNT(JP,NNSLA)
        KOBS = SLAB%NIND(JOBS)
        IF( KOBS .LT. 1 .OR. KOBS .GT. SLA%NO .OR. SLA%FLG( KOBS ) .EQ. 1 ) THEN
            IF( KOBS .LT. 1 .OR. KOBS .GT. SLA%NO ) CSTR='OUTR'
            IF( SLA%FLG( KOBS ) .EQ. 1            ) CSTR='DUPL'
            WRITE(333,*) 'ERROR : PROC, JOBS, KOBS, SLA%NO, FLG',&
            & JP, JOBS,KOBS,SLA%NO,SLA%FLG( KOBS ),CSTR
            IF( LL_AB ) CALL ABOR1('SLA ERROR')
        ENDIF 
        SLA%FLG( KOBS )= 1
        SLA%FLC( KOBS )=SLAB%FLC(JOBS)
        SLA%INO( KOBS )=SLAB%INO(JOBS)
        SLA%KSAT( KOBS )=SLAB%KSAT(JOBS)
        SLA%TRACK( KOBS )=SLAB%TRACK(JOBS)
        SLA%EVE( KOBS )=SLAB%EVE(JOBS)
        SLA%BOT( KOBS )=SLAB%BOT(JOBS)
        SLA%TDIST( KOBS )=SLAB%TDIST(JOBS)
        SLA%LON( KOBS )=SLAB%LON(JOBS)
        SLA%LAT( KOBS )=SLAB%LAT(JOBS)
        SLA%DPT( KOBS )=SLAB%DPT(JOBS)
        SLA%TIM( KOBS )=SLAB%TIM(JOBS)
        SLA%VAL( KOBS )=SLAB%VAL(JOBS)
        SLA%BIA( KOBS )=SLAB%BIA(JOBS)
        SLA%IB( KOBS,: )=SLAB%IB(JOBS,:)
        SLA%JB( KOBS,: )=SLAB%JB(JOBS,:)
        SLA%PQ( KOBS,: )=SLAB%PQ(JOBS,:)
      ENDDO CYCSLA

   ENDIF

   IST=1
   IEN=KKCNT(JP,NNSST)

   IF(KKCNT(JP,NNSST) .GT. 0 ) THEN

      WRITE(CFILE,'(A,I4.4,A)') 'SSTOBS_',JP-1,'.NC'

      CALL IO_OBS(CFILE=CFILE,IOPT=2,NOBS=KKCNT(JP,NNSST),&
      NLEVS=KM, NPREDS=NPRD_SST,NPQ=NPQ,NPQ2=NPQ,&
      FLC=SSTB%FLC(IST:IEN),&
      NIND=SSTB%NIND(IST:IEN),&
      TRACK=SSTB%TRACK(IST:IEN),&
      EVE=SSTB%EVE(IST:IEN),&
      KSAT=SSTB%KSAT(IST:IEN),&
      TDIST=SSTB%TDIST(IST:IEN),&
      LON=SSTB%LON(IST:IEN),&
      LAT=SSTB%LAT(IST:IEN),&
      TIM=SSTB%TIM(IST:IEN),&
      VAL=SSTB%VAL(IST:IEN),&
      BIA=SSTB%BIA(IST:IEN),&
      IB=SSTB%IB(IST:IEN,:),&
      JB=SSTB%JB(IST:IEN,:),&
      PQ=SSTB%PQ(IST:IEN,:))

      CYCSST : DO JOBS=1,KKCNT(JP,NNSST)
        KOBS = SSTB%NIND(JOBS)
        IF( KOBS .LT. 1 .OR. KOBS .GT. SST%NO .OR. SST%FLG( KOBS ) .EQ. 1 ) THEN
            IF( KOBS .LT. 1 .OR. KOBS .GT. SST%NO ) CSTR='OUTR'
            IF( SST%FLG( KOBS ) .EQ. 1            ) CSTR='DUPL'
            WRITE(334,*) 'ERROR : JOBS, KOBS, SST%NO, FLG',&
            & JOBS,KOBS,SST%NO,SST%FLG( KOBS ),CSTR
            IF( LL_AB ) CALL ABOR1('SST ERROR')
        ENDIF 
        SST%FLG( KOBS )= 1
        SST%FLC( KOBS )=SSTB%FLC(JOBS)
        SST%KSAT( KOBS )=SSTB%KSAT(JOBS)
        SST%TRACK( KOBS )=SSTB%TRACK(JOBS)
        SST%EVE( KOBS )=SSTB%EVE(JOBS)
        SST%TDIST( KOBS )=SSTB%TDIST(JOBS)
        SST%LON( KOBS )=SSTB%LON(JOBS)
        SST%LAT( KOBS )=SSTB%LAT(JOBS)
        SST%TIM( KOBS )=SSTB%TIM(JOBS)
        SST%VAL( KOBS )=SSTB%VAL(JOBS)
        SST%BIA( KOBS )=SSTB%BIA(JOBS)
        SST%IB( KOBS,: )=SSTB%IB(JOBS,:)
        SST%JB( KOBS,: )=SSTB%JB(JOBS,:)
        SST%PQ( KOBS,: )=SSTB%PQ(JOBS,:)
      ENDDO CYCSST

   ENDIF

   IST=1
   IEN=KKCNT(JP,NNSSS)

   IF(KKCNT(JP,NNSSS) .GT. 0 ) THEN

      WRITE(CFILE,'(A,I4.4,A)') 'SSSOBS_',JP-1,'.NC'

      CALL IO_OBS(CFILE=CFILE,IOPT=2,NOBS=KKCNT(JP,NNSSS),&
      NLEVS=KM, NPREDS=1,NPQ=NPQ,NPQ2=NPQ,&
      FLC=SSSB%FLC(IST:IEN),&
      TRACK=SSSB%TRACK(IST:IEN),&
      EVE=SSSB%EVE(IST:IEN),&
      KSAT=SSSB%KSAT(IST:IEN),&
      TDIST=SSSB%TDIST(IST:IEN),&
      LON=SSSB%LON(IST:IEN),&
      LAT=SSSB%LAT(IST:IEN),&
      TIM=SSSB%TIM(IST:IEN),&
      VAL=SSSB%VAL(IST:IEN),&
      BIA=SSSB%BIA(IST:IEN),&
      IB=SSSB%IB(IST:IEN,:),&
      JB=SSSB%JB(IST:IEN,:),&
      PQ=SSSB%PQ(IST:IEN,:))

      CYCSSS : DO JOBS=1,KKCNT(JP,NNSSS)
        KOBS = SSSB%NIND(JOBS)
        IF( KOBS .LT. 1 .OR. KOBS .GT. SSS%NO .OR. SSS%FLG( KOBS ) .EQ. 1 ) THEN
            IF( KOBS .LT. 1 .OR. KOBS .GT. SSS%NO ) CSTR='OUTR'
            IF( SSS%FLG( KOBS ) .EQ. 1            ) CSTR='DUPL'
            WRITE(335,*) 'ERROR : JOBS, KOBS, SSS%NO, FLG',&
            & JOBS,KOBS,SSS%NO,SSS%FLG( KOBS ),CSTR
            IF( LL_AB ) CALL ABOR1('SSS ERROR')
        ENDIF 
        SSS%FLG( KOBS )= 1
        SSS%FLC( KOBS )=SSSB%FLC(JOBS)
        SSS%KSAT( KOBS )=SSSB%KSAT(JOBS)
        SSS%TRACK( KOBS )=SSSB%TRACK(JOBS)
        SSS%EVE( KOBS )=SSSB%EVE(JOBS)
        SSS%TDIST( KOBS )=SSSB%TDIST(JOBS)
        SSS%LON( KOBS )=SSSB%LON(JOBS)
        SSS%LAT( KOBS )=SSSB%LAT(JOBS)
        SSS%TIM( KOBS )=SSSB%TIM(JOBS)
        SSS%VAL( KOBS )=SSSB%VAL(JOBS)
        SSS%BIA( KOBS )=SSSB%BIA(JOBS)
        SSS%IB( KOBS,: )=SSSB%IB(JOBS,:)
        SSS%JB( KOBS,: )=SSSB%JB(JOBS,:)
        SSS%PQ( KOBS,: )=SSSB%PQ(JOBS,:)
      ENDDO CYCSSS

   ENDIF

ENDDO

WRITE(IOUNLOG,*) ' INS READ, FAILURES : ',COUNT(INS%FLG(1:INS%NO).EQ.0)

WRITE(IOUNLOG,*)
WRITE(IOUNLOG,*) ' COUNTS FOR SLA SATELLITES'
DO JSAT=1,NOSLASATS
      SLA%ISLAOBS(JSAT) = COUNT( SLA%KSAT(1:SLA%NO) .EQ. JSAT)
      WRITE(IOUNLOG,*) ' SAT:',JSAT,' ',TRIM(CSLASATID(JSAT)),&
      & '   NO OF OBS:',SLA%ISLAOBS(JSAT)
ENDDO

WRITE(IOUNLOG,*)
WRITE(IOUNLOG,*) ' COUNTS FOR SST SATELLITES'
DO JSAT=1,NOSSTSATS
      SST%ISSTOBS(JSAT) = COUNT( SST%KSAT(1:SST%NO) .EQ. JSAT)
      WRITE(IOUNLOG,*) ' SAT:',JSAT,' ',TRIM(CSSTSATID(JSAT)),&
      & '   NO OF OBS:',SST%ISSTOBS(JSAT)
ENDDO

WRITE(IOUNLOG,*)
WRITE(IOUNLOG,*) ' COUNTS FOR SSS SATELLITES'
DO JSAT=1,NOSSSSATS
      SSS%ISSSOBS(JSAT) = COUNT( SSS%KSAT(1:SSS%NO) .EQ. JSAT)
      WRITE(IOUNLOG,*) ' SAT:',JSAT,' ',TRIM(CSSSSATID(JSAT)),&
      & '   NO OF OBS:',SSS%ISSSOBS(JSAT)
ENDDO

IF( INS%NO .GT. 0 ) THEN
   WHERE( INS%FLG .NE. 1 ) 
       INS%FLG = 0
       INS%FLC = 0
   END WHERE
   IF( SIZE(INS%FLC) .GT. INS%NO ) &
   INS%FLC( (INS%NO+1):SIZE(INS%FLC) ) = 0
ENDIF

CALL FLUSH(IOUNLOG)
CALL FLUSH(IOUNOUT)

IF ( INS%NO .GT. 0 ) THEN
  DEALLOCATE(INSB%INO,INSB%FLG,INSB%FLC,INSB%PAR)
  DEALLOCATE(INSB%LON,INSB%LAT,INSB%DPT,INSB%TIM)
  DEALLOCATE(INSB%VAL,INSB%BAC,INSB%INC,INSB%TDIST)
  DEALLOCATE(INSB%BIA,INSB%ERR,INSB%EVE,INSB%KTY)
  DEALLOCATE(INSB%RES,INSB%B_A,INSB%OTYPE,INSB%PROF)
  DEALLOCATE(INSB%IB,INSB%JB,INSB%KB,INSB%PLNO)
  DEALLOCATE(INSB%RB)
  DEALLOCATE(INSB%PQ)
  DEALLOCATE(INSB%PRIND,INSB%INST)
  DEALLOCATE(INSB%NIND)
ENDIF


IF ( INS%NO .GT. 0 ) THEN
  KOBS=INS%NO
  INS%BAC(:) = 0._R8
  INS%RES(:) = 0._R8
  OPEN(31,FILE='INSOBS_ALL.bin',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
  WRITE(31) INS%FLG(1:KOBS)
  WRITE(31) INS%FLC(1:KOBS)
  WRITE(31) INS%INO(1:KOBS)
  WRITE(31) INS%OTYPE(1:KOBS)
  WRITE(31) INS%PAR(1:KOBS)
  WRITE(31) INS%PLNO(1:KOBS)
  WRITE(31) INS%INST(1:KOBS)
  WRITE(31) INS%EVE(1:KOBS)
  WRITE(31) INS%KTY(1:KOBS)
  WRITE(31) INS%PROF(1:KOBS)
  WRITE(31) INS%TDIST(1:KOBS)
  WRITE(31) INS%LON(1:KOBS)
  WRITE(31) INS%LAT(1:KOBS)
  WRITE(31) INS%DPT(1:KOBS)
  WRITE(31) INS%TIM(1:KOBS)
  WRITE(31) INS%VAL(1:KOBS)
  WRITE(31) INS%BIA(1:KOBS)
  WRITE(31) INS%IB(1:KOBS,:)
  WRITE(31) INS%JB(1:KOBS,:)
  WRITE(31) INS%KB(1:KOBS)
  WRITE(31) INS%RB(1:KOBS)
  WRITE(31) INS%PQ(1:KOBS,:)
  WRITE(31) INS%BAC(1:KOBS)
  WRITE(31) INS%RES(1:KOBS)
  CLOSE(31)

  OPEN(31,FILE='INSOBS_CR.bin',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
  CALL INT_PAR_INS_2(IM,JM,KM)
  WRITE(31) INS%IB(1:KOBS,:)
  WRITE(31) INS%JB(1:KOBS,:)
  WRITE(31) INS%PQ(1:KOBS,:)
  CLOSE(31)

ENDIF

IF ( SLA%NO .GT. 0 ) THEN
  KOBS=SLA%NO

  SLA%BAC = 0._R8
  SLA%TB = 0._R8
  SLA%SB = 0._R8
  SLA%BCP = 0._R8
  SLA%RES = 0._R8
  OPEN(32,FILE='SLAOBS_ALL.bin',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
  WRITE(32) SLA%FLC(1:KOBS)
  WRITE(32) SLA%NIND(1:KOBS)
  WRITE(32) SLA%INO(1:KOBS)
  WRITE(32) SLA%TRACK(1:KOBS)
  WRITE(32) SLA%EVE(1:KOBS)
  WRITE(32) SLA%KSAT(1:KOBS)
  WRITE(32) SLA%BOT(1:KOBS)
  WRITE(32) SLA%TDIST(1:KOBS)
  WRITE(32) SLA%LON(1:KOBS)
  WRITE(32) SLA%LAT(1:KOBS)
  WRITE(32) SLA%TIM(1:KOBS)
  WRITE(32) SLA%VAL(1:KOBS)
  WRITE(32) SLA%BIA(1:KOBS)
  WRITE(32) SLA%IB(1:KOBS,:)
  WRITE(32) SLA%JB(1:KOBS,:)
  WRITE(32) SLA%PQ(1:KOBS,:)
  WRITE(32) SLA%DPT(1:KOBS)
  WRITE(32) SLA%BAC(1:KOBS)
  WRITE(32) SLA%TB(1:KOBS,1:KM)
  WRITE(32) SLA%SB(1:KOBS,1:KM)
  WRITE(32) SLA%BCP(1:KOBS,1:NPRD_SLA)
  WRITE(32) SLA%RES(1:KOBS)
  CLOSE(32)

  OPEN(32,FILE='SLAOBS_CR.bin',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
  WRITE(32) SLA%IB(1:KOBS,:)
  WRITE(32) SLA%JB(1:KOBS,:)
  WRITE(32) SLA%PQ(1:KOBS,:)
  CLOSE(32)
ENDIF

IF ( SST%NO .GT. 0 ) THEN
  SST%BAC = 0._R8
  SST%ERR = 0._R8
  SST%BCP = 0._R8
  SST%RES = 0._R8
  KOBS=SST%NO
  OPEN(33,FILE='SSTOBS_ALL.bin',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
  WRITE(33) SST%FLC(1:KOBS)
  WRITE(33) SST%TRACK(1:KOBS)
  WRITE(33) SST%EVE(1:KOBS)
  WRITE(33) SST%KSAT(1:KOBS)
  WRITE(33) SST%TDIST(1:KOBS)
  WRITE(33) SST%LON(1:KOBS)
  WRITE(33) SST%LAT(1:KOBS)
  WRITE(33) SST%TIM(1:KOBS)
  WRITE(33) SST%VAL(1:KOBS)
  WRITE(33) SST%BIA(1:KOBS)
  WRITE(33) SST%IB(1:KOBS,:)
  WRITE(33) SST%JB(1:KOBS,:)
  WRITE(33) SST%PQ(1:KOBS,:)
  WRITE(33) SST%BAC(1:KOBS)
  WRITE(33) SST%ERR(1:KOBS)
  WRITE(33) SST%BCP(1:KOBS,1:NPRD_SLA)
  WRITE(33) SST%RES(1:KOBS)
  CLOSE(33)

  OPEN(33,FILE='SSTOBS_CR.bin',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
  WRITE(33) SST%IB(1:KOBS,:)
  WRITE(33) SST%JB(1:KOBS,:)
  WRITE(33) SST%PQ(1:KOBS,:)
  CLOSE(33)
ENDIF

IF ( SSS%NO .GT. 0 ) THEN
  SSS%BAC = 0._R8
  SSS%RES = 0._R8
  KOBS=SSS%NO
  OPEN(34,FILE='SSSOBS_ALL.bin',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
  WRITE(34) SSS%FLC(1:KOBS)
  WRITE(34) SSS%TRACK(1:KOBS)
  WRITE(34) SSS%EVE(1:KOBS)
  WRITE(34) SSS%KSAT(1:KOBS)
  WRITE(34) SSS%TDIST(1:KOBS)
  WRITE(34) SSS%LON(1:KOBS)
  WRITE(34) SSS%LAT(1:KOBS)
  WRITE(34) SSS%TIM(1:KOBS)
  WRITE(34) SSS%VAL(1:KOBS)
  WRITE(34) SSS%BIA(1:KOBS)
  WRITE(34) SSS%IB(1:KOBS,:)
  WRITE(34) SSS%JB(1:KOBS,:)
  WRITE(34) SSS%PQ(1:KOBS,:)
  WRITE(34) SSS%BAC(1:KOBS)
  WRITE(34) SSS%RES(1:KOBS)
  CLOSE(34)

  OPEN(34,FILE='SSSOBS_CR.bin',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
  WRITE(34) SSS%IB(1:KOBS,:)
  WRITE(34) SSS%JB(1:KOBS,:)
  WRITE(34) SSS%PQ(1:KOBS,:)
  CLOSE(34)
ENDIF

WRITE(*,*) 'TOTALS : ',INS%NO, SLA%NO, SST%NO, SSS%NO

OPEN(35,FILE='OBS_TAB.BIN')
WRITE(35,*) INS%NO, SLA%NO, SST%NO, SSS%NO
CLOSE(35)

END PROGRAM READOBS103
