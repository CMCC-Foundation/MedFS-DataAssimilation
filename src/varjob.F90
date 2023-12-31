SUBROUTINE VARJOB

!-----------------------------------------------------------------------
!                                                                      !
! THE MAIN DRIVER OF THE VARIATIONAL ANALYSIS                          !
!                                                                      !
! VERSION 1: S.DOBRICIC 2006                                           !
!   HORIZONTAL COVARIANCE WITH RECURSIVE FILTERS, VERTICAL WITH EOFS,  !
!   ASSIMILATION OF SATELLITE OBSERVATIONS OF SLA, IN SITU OBSERVATIONS!
!   BY XBT AND ARGO FLOATS                                             !
!                                                                      !
! VERSION 2: S.DOBRICIC 2006                                           !
!   MULTIGRID METHOD.                                                  !
!                                                                      !
! VERSION 3: S.DOBRICIC 2006                                           !
!   INTERNAL BOUNDARIES FOR HORIZONTAL                                 !
!                                                                      !
! VERSION 4: A.STORTO   2009-2010                                      !
!   REORGANIZATION OF CODE ARCHITECTURE                                !
!   SUPPORT GLOBAL TREATMENT OF INSITU,SLA,SST,SSS                     !
!   SUPPORT NEW OBSERVATION PREPROCESSING                              !
!   SUPPORT SPECIAL TREATMENT OF NON-HOMOG,ANISOTR RF COEFFS           !
!   SUPPORT DIFFERENT CONFIGURATION AND DIFFERENT FORMATS              !
!   EXPERIMENTAL SUPPORT OF SEVERAL MINIMIZERS                         !
!                                                                      !
! VERSION 5: A.STORTO   2011
!   CLEANING AND SUPPRESSION OBSOLETE FEATURES
!-----------------------------------------------------------------------

 USE SET_KND
 USE DRV_STR
 USE EOF_STR, ONLY : NPRINTEOF
 USE CTL_STR
 USE IOUNITS
 USE ICE
 USE RUN
 USE MYFRTPROF, ONLY : MYFRTPROF_WALL, MYFRTPROF_PRINT
 USE OBSHANDLING, ONLY : LLUSEBGERR
 USE EXTGRID, ONLY : LLEXTGRID, DOEXTGRID,UNDOEXTGRID,DOOBSREARR,&
 & UNDOOBSREARR,OLD_IM,NEW_IM
 USE RECFILTER, ONLY : LLRFSR,NNRFSR
 USE MPIREL
 USE BAL_AIRSEA
 USE TRACK_CORREL
 USE PROF_CORREL
 USE TLAD_VARS
 USE TLAD
 USE TLAD_TEST
 USE HYBRID
 USE EOFS3D
 USE CERES

 IMPLICIT NONE

 INTEGER(I4)   ::  KTR
 INTEGER(I4)   ::  IERR

 CALL MYFRTPROF_WALL('VARJOB: VARIATIONAL JOB DRIVER',0)

! OUTER LOOP - MULTIGRID
 DO KTR = 1,DRV%NTR

   DRV%KTR = KTR

! DEFINE GRID PARAMETERS
   CALL RDGRDS(DRV%CGRID(KTR))
   CALL SETFGFLAYOUT
   PRINT*,' GRID READ-IN'

   IF(NCONF.EQ.101) CALL READ_FIRSTGUESS

   IF(NCONF.EQ.101 .OR. NCONF.EQ.102) THEN

! GET OBSERVATIONS
     IF(KTR.EQ.1) THEN
       IF(COBSOPT.EQ.'GLO') THEN
         CALL READOBS
       ELSEIF(COBSOPT.EQ.'ASC') THEN
         CALL READOBSA
       ELSE
         WRITE(IOUNERR,*) 'VARJOB: OPTION ',COBSOPT,' IS NOT RECOGNIZED!'
         CALL ABOR1('CANNOT FIND A WAY TO READ OBSERVATIONS')
       ENDIF
       PRINT*,' OBS READ-IN'
     ENDIF

! DEFINE INTERPOLATION PARAMETERS
     IF(COBSOPT.EQ.'GLO' .OR. COBSOPT.EQ.'ASC') THEN
       CALL INT_PAR2
     ENDIF
     PRINT*,' INTERPOLATION SET-UP'

! COMPUTE FIRST GUESS DEPARTURES OR EXIT

     IF(NCONF.EQ.102) THEN
       IF( CINS_FMT(1:3) .EQ. 'NRT' ) CALL REDCHECK
       CALL REPEVE
       CALL WRITEOBS102
       CALL MYFRTPROF_WALL('VARJOB: VARIATIONAL JOB DRIVER',1)
       RETURN
     ENDIF

     IF(NCONF.EQ.101) THEN
       CALL INT_BACVAL_SLA
       CALL OBS_INS
       CALL OBS_SST
       CALL OBS_SSS
       CALL UNBIAS_SLA
       IF( LLOBSSTAT ) CALL OBSSTATS('CONF 101: AFTER CALLING OBSOP',.FALSE.)
       CALL DEALLBG
       CALL REO4MIN
     ENDIF

   ENDIF

   IF(NCONF.EQ.103) THEN
     IF(IDOUBLEDOM .EQ. 2 ) THEN
       CALL READOBS103DD
     ELSE
       IF(COBSOPT.EQ.'ASC') THEN
         CALL READOBSA103
       ELSE
         IF( LL_READALLBIN103 ) THEN
           CALL READOBS103BIN
         ELSE
           CALL ABOR1('READOBS103 NOT SUPPORTED ANYMORE')
           CALL READOBS103
         ENDIF
         IF( LLOBSSTAT ) CALL OBSSTATS('CONF 103: AFTER READING MISFITS FILES',.TRUE.)
         CALL PREPPROF
       ENDIF
     ENDIF
   ENDIF

   IF( NCONF .EQ. 202 ) THEN
     CALL PREP_C202
     CALL ASSIGN_SLAERR
   ENDIF

   IF(NCONF.EQ.101 .OR. NCONF.EQ.103 .OR. NCONF.EQ.105 .OR. NCONF.EQ.401) THEN
     WRITE(IOUNOUT,*) ' EOF SETUP'
     CALL SUEOF
     IF(NPRINTEOF .GT. 0) CALL PRINT_EOFS
     IF( LL_HYBRID ) CALL SUEOFH
     IF( LL_EOFS3D ) CALL SUEOF3D
     WRITE(IOUNOUT,*) ' DONE'
   ENDIF

! OBS PREPROC CHAIN

   IF( ( COBSOPT.EQ.'GLO' .OR. COBSOPT.EQ.'ASC' ) .AND. &
       ( NCONF .EQ. 101 .OR. NCONF .EQ. 103 ) .AND. IDOUBLEDOM .LT. 2) THEN
     WRITE(IOUNOUT,*) ' OBS PREPROCESSING CHAIN'
     CALL CHECK_TIME
     CALL BLACKLIST
     IF( SSS%NO .GT. 0) THEN
       CALL PREPSSS
       CALL BIASCORR
     ENDIF
     CALL COASTREJ
     CALL OBSERRORS

     !... SLA BIAS CORRECTION
     IF(LL_SLABCORR .AND. NCONF .EQ. 103) CALL SLA_BCORR
     !... INDEPENDENT DECISIONS
     CALL CLIMCHECK
     IF(LLUSEBGERR) THEN
       CALL BGSDE_OBS
       CALL BGQCHECKP
     ELSE
       CALL BGQCHECK
     ENDIF

     IF(LLSLAATCORR) THEN
       WRITE(IOUNOUT,*) ' CALLING SLAATCORR'
       CALL SLAATCORR
       CALL MYFRTPROF_WALL('VARJOB: VARIATIONAL JOB DRIVER',1)
       RETURN
     ENDIF

     IF( NN_REG_REJ .GT. 0 ) CALL REG_REJ(NN_REG_REJ)

     !... CHECK BACKGROUND OK FOR SLA
     CALL CHECK_SLABG(1)
     !... REJECTION UNDER MODEL SEA-ICE
     IF(LL_ICEREJ) CALL ICEREJECT
     !... DEPENDENT DECISIONS
     !... REJECTIONS ABOVE
     IF(LL_VERTCHECK) CALL VERTCHECK
     !... THINNING
     CALL THINNING
     WRITE(IOUNOUT,*) ' DONE'
   ENDIF

   IF( NCONF .EQ. 202 ) THEN
     IF(LLUSEBGERR) THEN
       CALL BGSDE_OBS
       CALL BGQCHECKP
     ELSE
       CALL BGQCHECK
     ENDIF
   ENDIF

   IF( NCONF.EQ.101) THEN
     IF( LL_UNBIAS_SLA_AT ) THEN
       CALL UNBIAS_SLA_AT
     ELSE
       CALL UNBIAS_SLA
     ENDIF
   ENDIF


   IF( ( COBSOPT.EQ.'GLO' .OR. COBSOPT.EQ.'ASC' ) .AND. IDOUBLEDOM .LT. 2 &
     & .AND. ( NCONF .EQ. 101 .OR. NCONF .EQ. 103 .OR. NCONF .EQ. 202) ) THEN
      !... REPORT EVENTS
     WRITE(IOUNOUT,*) ' EVENTS REPORTING'
     CALL REPEVE
     WRITE(IOUNOUT,*) ' DONE'
   ENDIF

   ! HERE TO WRITE STATS
   IF( LL_WRITE_OBS_SCREEN .AND. IDOUBLEDOM .LT. 2 .AND. NCONF .NE. 301) &
       & CALL WRT_OBS_SCREEN
    !Jenny change code in order to evaluate misfit  in simulation
   IF (LL_MISFIT_SIMU) THEN
     CALL TERMIN0
     STOP
   ENDIF
   IF( LL_HYBRID .AND. LL_ALPHA_ONLINE ) CALL ALPHA_ONLINE

   ! STOP IF ONLY PREPROC
   IF(LLONLYOBSPP) THEN
     CALL MYFRTPROF_WALL('VARJOB: VARIATIONAL JOB DRIVER',1)
     RETURN
   ENDIF
    
   IF (NCONF.NE.104 .AND. NCONF.NE.105 .AND. NCONF.NE.301 .AND. NCONF.NE.401) THEN

     IF( LLOBSSTAT .AND. IDOUBLEDOM .LT. 2 ) CALL OBSSTATS('BEFORE MINIMIZATION',.FALSE.)

       ! DEFINE OBSERVATIONAL VECTOR
       IF(COBSOPT.EQ.'MFS') THEN
         CALL OBS_VEC
       ELSEIF(COBSOPT.EQ.'GLO' .OR. COBSOPT.EQ.'ASC') THEN
         WRITE(IOUNOUT,*) ' OBS REARRAINGING'
         IF(LLEXTGRID) THEN
           IF( IDOUBLEDOM .LT. 2 ) CALL REO4MIN
           CALL DOOBSREARR
         ENDIF
         CALL OBS_VEC2
         IF( SLA%NO .EQ. 0) NN_TRACK_CORREL = 0
         IF( NN_TRACK_CORREL .GT. 0 ) THEN
           CALL PREP_SLARINV
         ENDIF
         IF( INS%NO .EQ. 0) NN_PROF_CORREL = 0
         IF( NN_PROF_CORREL .GT. 0 ) THEN
           CALL PREP_INSRINV
         ENDIF
       ENDIF

       WRITE(IOUNOUT,*) ' DONE'
     ENDIF

! DEFINE CONSTANTS FOR BACKGROUND COVARIANCES
! AND EVENTUALLY EXTEND THE DOMAIN

     IF ( NCONF.NE.104 .AND. NCONF.NE.202 .AND. NCONF.NE.301) THEN
       IF( LL_MDEOFS ) THEN 
         CALL MDEOFS
         IF(NPRINTEOF .GT. 0) CALL PRINT_EOFS
       ENDIF
     ENDIF

     IF(NCONF.EQ.105) THEN
       CALL MYFRTPROF_WALL('VARJOB: VARIATIONAL JOB DRIVER',1)
     RETURN
   ENDIF

   CALL PREP_EXTGRID
   IF(LL_HYBRID_CR ) CALL PREP_EXTGRIDF

   CALL TLAD_INIT
   IF(LL_ADTEST) CALL TLADMOD_TEST

! INITIALIZE SIMPLIFIED AIR-SEA BALANCE MODEL
   IF(DRV%BA2(DRV%KTR) .EQ. 1) THEN
     CALL INIT_BAL_AS
     PRINT*,' BAL_AS INITIALISED'
   ENDIF

   CALL WEAKLY_INIT

   IF(NCONF .EQ. 301) CALL PREP301
   IF(LLEXTGRID) CALL DOEXTGRID

   IF(NCONF .EQ. 301) THEN
     CALL TLMOD(1)
     CALL MYFRTPROF_WALL('VARJOB: VARIATIONAL JOB DRIVER',1)
     RETURN
   ENDIF

   CALL TLAD_INIT_OBS

   CALL RF_INIT

   IF(LL_HYBRID_CR ) CALL SURF_LVF
   CALL SURF_LV

   IF(NCONF.EQ.104) THEN
     CALL MYFRTPROF_WALL('VARJOB: VARIATIONAL JOB DRIVER',1)
     RETURN
   ENDIF

   IF( NCONF .NE. 202 ) CALL CHECK_SLABG(2)

! INITIALISE BAROTROPIC MODEL
   IF(DRV%BMD(DRV%KTR) .EQ. 1) THEN
     CALL INI_BMD
     PRINT*,' BM INITIALISED'
   ENDIF

! INITIALIZE SIMPLIFIED BALANCE MODEL
   IF(DRV%BAL(DRV%KTR) .EQ. 1) THEN
     CALL INI_BAL
     PRINT*,' BAL INITIALISED'
   ENDIF

   IF(LL_CERES) CALL CERES_INIT

! INITIALIZE COST FUNCTION AND ITS GRADIENT
   CALL INI_CFN
   PRINT*,' J INITIALISED'

! INITIALISE FROM OLD ITERRATION
   IF(KTR.GT.1) THEN
     CALL INI_ITR
     PRINT*,' J INITIALISED FROM COARSER GRID'
   ENDIF

! BACKGROUND PERTURBATION
   IF(NCONF.EQ.401) THEN
     CALL BGPERT
     CALL MYFRTPROF_WALL('VARJOB: VARIATIONAL JOB DRIVER',1)
     RETURN
   ENDIF


! MINIMIZE THE COST FUNCTION (INNER LOOP)
   CALL MIN_CFN
   PRINT*,' J MINIMIZED'
 
#ifndef NOMPI
   IF( LLUSEMPI ) CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
#endif

! CONVERT TO INNOVATIONS
   CALL CNV_INN
   PRINT*,' CONTROL VECTOR MOVED INTO PHYSICAL SPACE'

! FROM EXTENDED TO PHYSICAL GRID
   IF(LLEXTGRID) THEN
     CALL UNDOEXTGRID
     CALL UNDOOBSREARR
   ENDIF

! WRITE OUTPUTS AND DIAGNOSTICS
   IF( LL_ICEANFILT ) CALL ICEANFILT
   CALL WRT_OUT
   PRINT*,' POST-PROCESSING COMPLETED'

! SAVE OLD ITERRATION
!      CALL SAV_ITR
!      PRINT*,' MULTI-GRID ITERATION SAVED'

! END OF OUTER LOOP
 ENDDO

!-----------------------------------------------------------------

 CALL MYFRTPROF_WALL('VARJOB: VARIATIONAL JOB DRIVER',1)
END SUBROUTINE VARJOB
