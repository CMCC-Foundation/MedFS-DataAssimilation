SUBROUTINE THINOUT_INS_A(ITER,ZIN1,ZIN2,ZIN3)

!.. THINNNING OF INS IN-SITU OBSERVATIONS
!..
!.. METHOD A: FIND REPORTS FROM THE SAME
!.. PLATFORMS AND REJECT PROFILES FURTHER
!.. FROM ANALYSIS TIME.
!.. THIS ROUTINE SHOULD BE CALLED TWICE
!.. (ITER=1,2) TO ADD A LAT/LON OFFSET
!.. IN ORDER TO PERMIT REJECTION CLOSE
!.. TO THE THINNING BOX BOUNDARIES.
!.. THE PROCEDURE IS PERFORMED REGARDLESS
!.. OF THE NUMBER OF ACTIVE OBS FOR EACH
!.. PROFILE, ASSUMING THAT THE QUALITY
!.. OF EACH PLATFORM DOES NOT CHANGE
!.. SIGNIFICANTLY WITHIN THE ASSIMILATION
!.. WINDOW.
!.. ROUGHLY PERFORM VERTICAL THINNING.
!..
!.. ANDREA STORTO *INGV* 2009-04-6

USE SET_KND
USE OBSDEF
USE OBS_STR
USE GRD_STR
USE MYNETCDF
USE IOUNITS
USE MYFRTPROF

IMPLICIT NONE

 INTEGER(I4),INTENT(IN) :: ITER
 REAL(R8),INTENT(IN) :: ZIN1,ZIN2,ZIN3
 REAL(R8) :: XRES,YRES,ZAUX,ZM,ZRES,LAST_DEP
 INTEGER(I4)   ::  K,JOBS,IIND,JIND,KOKTH,KKOTH,NX,NY,NZ,MAXOBS,OIND
 INTEGER(I4)   ::  JX,JY,JZ,JTYPE,ZIND,JPAR,I1,J1,K1,JINCR,PKN
 REAL(R8) :: ZDIST2,P1,Q1,R1
 REAL(DP) :: ZDIST
 INTEGER(I4), ALLOCATABLE :: KCOUNTTH(:,:),KOBS(:,:,:)
 INTEGER(I4) :: JDUP, KDUPL, KK, JPROF, K2, IND1, IND2,KHORTH,&
              & NODP,KGOOD,KS,KE,KVERTH,KEVE_HTHN
 LOGICAL :: LVERTTHIN, LLDBG
 LOGICAL, PARAMETER :: LLDEBUG = .FALSE.
 INTEGER(I4),PARAMETER :: MAXDUPL=1000
 CHARACTER(LEN=8) :: CDUPL(MAXDUPL)
 INTEGER(I4) :: LAST_GOOD, JCTRL, IIND2, JIND2

#include "obs_events.h"

 CALL MYFRTPROF_WALL('THINOUT_INS_A: INS DATA THINNING A',0)

  WRITE(IOUNOUT,*) ' *** INS THINNING - ITERATION ',ITER
  WRITE(IOUNLOG,*)
  WRITE(IOUNLOG,*) ' *** INS THINNING - ITERATION ',ITER

  LVERTTHIN = ( ITER == 2 )
  KOKTH=0
  KKOTH=0
  KHORTH=0
  KVERTH=0
  NODP=0
  XRES=ZIN1
  YRES=ZIN2
  ZRES=ZIN3

  IF(ITER.EQ.1) THEN
     PKN=1
     KEVE_HTHN = KEVE_HTH1
     LLDBG = LLDEBUG
  ELSEIF(ITER.EQ.2) THEN
     PKN=2
     KEVE_HTHN = KEVE_HTH2
     LLDBG = .FALSE.
  ELSE
     CALL ABOR1('THINOUT: UNSUPPORTED NUMBER OF ITERATION')
  ENDIF

  NX = INT(GRD%IM/ZIN1) + 1
  NY = INT(GRD%JM/ZIN2) + 1

  MAXOBS = 5000

  WRITE(IOUNLOG,*)
  WRITE(IOUNLOG,*) ' THINNING CONFIGURATION'
  WRITE(IOUNLOG,*) ' ITERATION  ', ITER
  WRITE(IOUNLOG,*) ' NX,  NY    ', NX, NY
  WRITE(IOUNLOG,*) ' XRES, YRES ', XRES, YRES
  WRITE(IOUNLOG,*) '       ZRES ', ZRES
  WRITE(IOUNLOG,*) '       PKN  ', PKN
  WRITE(IOUNLOG,*) '      MAXOBS', MAXOBS

  ALLOCATE(KCOUNTTH(NX,NY),  &
         & KOBS(NX,NY,MAXOBS) )

!... STEP 1. IDENTIFY STATIONS IN SAME REGION

  KCOUNTTH=0

  WRITE(IOUNLOG,*) ' NUMBER OF PROFILES :', INS%NPROFS
  CALL FLUSH(IOUNLOG)


  EPROF : DO JPROF=1,INS%NPROFS

    K =  INS%PRIND(JPROF)

    IF(JPROF.LT.INS%NPROFS) THEN
     K2 = INS%PRIND(JPROF+1)-1
    ELSE
     K2 = INS%NO
    ENDIF

    IF( ALL( INS%FLC(K:K2) .EQ. 0) ) THEN
         NODP=NODP+1
         CYCLE EPROF
    ENDIF

    JINCR=0
    DO WHILE ( INS%FLC(K+JINCR).NE.1)
       JINCR=JINCR+1
    ENDDO

    IF(K+JINCR.GE.K2) THEN
       INS%FLC(K:K2)=0
       CYCLE EPROF
    ENDIF

    IIND = INT( (REAL(INS%IB(K+JINCR,PKN),KIND=R8))/XRES) + 1
    JIND = INT( (REAL(INS%JB(K+JINCR,PKN),KIND=R8))/YRES) + 1

    IF(IIND.GT.NX .OR. JIND.GT.NY .OR. IIND.LT.1 .OR. JIND.LT.1) THEN
       WRITE(IOUNERR,*) 'INDEXES: ',IIND, JIND
       CALL ABOR1('INDEX ERROR IN THINOUT_INS_A')
    ENDIF

    KCOUNTTH(IIND,JIND) = KCOUNTTH(IIND,JIND) + 1
    OIND = KCOUNTTH(IIND,JIND)
    IF(OIND.GT. MAXOBS) &
    & CALL ABOR1('THINOUT: MAXOBS MUST BE INCREASED')
    KOBS(IIND,JIND,OIND) = K

  ENDDO EPROF

  WRITE(IOUNLOG,*) ' PROFS, NODP, SUM(KCOUNT)',INS%NPROFS,NODP,SUM(KCOUNTTH)
  CALL FLUSH(IOUNLOG)

  IF(LLDBG) CALL WRITE_VARNCDF('KCOUNTTH.NC',NX,NY,KCOUNTTH,'KCOUNTTH','n',&
             & 'OBS COUNTS')

!... LOOP OVER THE GRIDS

   IF(LLDBG) OPEN(129,FILE='thinning_dupl.txt')

   DO JY=1,NY
    CYX : DO JX=1,NX

      IF(KCOUNTTH(JX,JY) .EQ. 0 ) CYCLE CYX

!... STEP 2. IDENTIFY DUPLICATED REPORTS WITHIN THE NX*NY THINNING BOXES

            KDUPL=0
            OUTER : DO K= 1, KCOUNTTH(JX,JY)
               IND1 = KOBS(JX,JY,K)
               IF(KDUPL.GT.0) THEN
                  DO KK=1,KDUPL
                     IF(CDUPL(KK)(1:8) .EQ. INS%PLNO( IND1 )(1:8) ) &
                     CYCLE OUTER
                  ENDDO
               ENDIF

               INNER : DO K2=1, KCOUNTTH(JX,JY)
                 IND2 = KOBS(JX,JY,K2)
                 IF(IND2.EQ.IND1) CYCLE INNER

                 IF(INS%PLNO( IND1 )(1:8) .EQ. INS%PLNO( IND2 )(1:8) ) THEN

                    KDUPL=KDUPL+1

                    IF(KDUPL.GT.MAXDUPL) &
                    & CALL ABOR1('THINOUT : MAXDUPL MUST BE INCREASED')

                    CDUPL(KDUPL)(1:8) = INS%PLNO( IND1 )(1:8)
                    CYCLE OUTER

                 ENDIF

               ENDDO INNER
            ENDDO OUTER

            IF(LLDBG) THEN
               WRITE(129,*) JX, JY, KDUPL
               IF(KDUPL.GT.0) THEN
                  DO JDUP=1,KDUPL
                   WRITE(129,*) JDUP, ' >>'//CDUPL(JDUP)(1:8)//'<<'
                  ENDDO
               ENDIF
            ENDIF

!... STEP 3. DECISION: ONLY THE PROFILE CLOSER TO ANTIME IS RETAINED

            DUPL : DO JDUP=1,KDUPL

               ZDIST=10000000._DP
               KGOOD=0
               CALL FLUSH(6)
               DO K2=1, KCOUNTTH(JX,JY)
                  IND2 = KOBS(JX,JY,K2)
                  IF (CDUPL(JDUP)(1:8) .EQ. INS%PLNO( IND2 )(1:8)) THEN
                      IF(ABS(INS%TDIST( IND2 )).LT.ZDIST) THEN
                         KGOOD=INS%PROF( IND2 )
                         ZDIST=ABS(INS%TDIST( IND2 ))
                      ENDIF
                  ENDIF
               ENDDO
               IF(KGOOD.EQ.0) CALL ABOR1('INS THINNING')

               IF(LLDBG) THEN
                   WRITE(129,*) JDUP,CDUPL(JDUP)(1:8),ZDIST, KGOOD
               ENDIF

!... STEP 4. CORRECTION: REJECT OBSERVATIONS ACCORDINGLY

               DO JPROF=1,INS%NPROFS

                  KS=INS%PRIND(JPROF)

                  IF(JPROF.LT.INS%NPROFS) THEN
                    KE = INS%PRIND(JPROF+1)-1
                  ELSE
                    KE = INS%NO
                  ENDIF

                  IF(CDUPL(JDUP)(1:8) .EQ. INS%PLNO(KS)(1:8))THEN
                     IF(LLDBG) THEN
                       WRITE(129,*) '-- ',JDUP,KS,INS%PROF(KS)
                       CALL FLUSH(129)
                     ENDIF
                     IIND2 = INT( (REAL(INS%IB(KS,PKN),KIND=R8))/XRES) + 1
                     JIND2 = INT( (REAL(INS%JB(KS,PKN),KIND=R8))/YRES) + 1
                     IF(INS%PROF(KS).NE.KGOOD.AND.JX.EQ.IIND2.AND.JY.EQ.JIND2) THEN
                       DO K2=KS,KE
                          ! IF ALREADY REJECTED DO NOT CHANGE EVE
                          IF(INS%FLC(K2).EQ.1) THEN
                             INS%FLC(K2) = 0
                             INS%EVE(K2) = KEVE_HTHN
                             KHORTH = KHORTH + 1
                          ENDIF
                       ENDDO
                     ELSE
                          KOKTH = KOKTH+SUM(INS%FLC(KS:KE))

!... STEP 5. ROUGHLY PERFORM VERTICAL THINNING

                        IF (LVERTTHIN) THEN
                          VTHIN : DO K2=KS+1,KE
                             IF(INS%FLC(K2).EQ.0) CYCLE VTHIN
                             LAST_DEP=-1._R8
                             LAST_GOOD=0
                             CONTROL_LOOP: DO JCTRL=K2-1,KS,-1
                                IF( INS%FLC(JCTRL).EQ.1 .AND. INS%PAR(JCTRL).EQ.INS%PAR(K2) ) THEN
                                   LAST_DEP=INS%KB(JCTRL)+INS%RB(JCTRL)
                                   LAST_GOOD=JCTRL
                                   EXIT CONTROL_LOOP
                                ENDIF
                             ENDDO CONTROL_LOOP
                             IF( LAST_GOOD .EQ. 0 ) CYCLE VTHIN
                             ZDIST2=ABS( REAL(INS%KB(K2))+INS%RB(K2) - &
                                      & LAST_DEP )
                             IF(ZDIST2.LT.ZRES) THEN
                                INS%FLC(K2) = 0
                                INS%EVE(K2) = KEVE_VTHN
                                KVERTH = KVERTH + 1
                             ENDIF
                          ENDDO VTHIN
                        ENDIF


                     ENDIF
                  ENDIF
               ENDDO
            ENDDO DUPL

     ENDDO CYX
  ENDDO

  IF(LLDBG) CLOSE(129)

  WRITE(IOUNLOG,*) ' KHORTH' ,KHORTH
  WRITE(IOUNOUT,*) ' END INS THINNING - ITERATION ',ITER
  WRITE(IOUNLOG,*)
  WRITE(IOUNLOG,*) ' *** THINNING REPORT - ITERATION ',ITER
  WRITE(IOUNLOG,*) ' THINNING GRID : ',NX,' X',NY
  WRITE(IOUNLOG,*) ' FILTERED-OUT HORIZONTALLY : ',KHORTH
  IF (ITER.EQ.2) &
  & WRITE(IOUNLOG,*) ' FILTERED-OUT VERTICALLY   : ',KVERTH
  WRITE(IOUNLOG,*) ' RETAINED      : ',KOKTH
  WRITE(IOUNLOG,*) ' TOTAL         : ',KOKTH + KKOTH
  WRITE(IOUNLOG,*) ' FLAGS COUNT   : ',SUM(INS%FLC(1:INS%NO))

  INS%NC = KOKTH

  DEALLOCATE(KCOUNTTH, KOBS )

CALL FLUSH(IOUNLOG)
CALL FLUSH(IOUNOUT)

CALL MYFRTPROF_WALL('THINOUT_INS_A: INS DATA THINNING A',1)
END SUBROUTINE THINOUT_INS_A
