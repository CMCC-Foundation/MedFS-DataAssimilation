SUBROUTINE READWODFILE(CFIN,IOUN,NMAX,NOB,KTYP,PAR,DEP,VAL,TIM,LON,LAT,&
& PRF,PLA,INS)

! READ OBSERVATION IN REFORMATTED WOD FORMAT
!
! A.S. - 26.11.2013

  USE SET_KND
  USE CALENDAR
  USE OBSDEF
  USE RUN , ONLY : ZJULSTART, ZJULEND

  IMPLICIT NONE

  INTEGER(I4), INTENT(IN) :: IOUN,NMAX
  CHARACTER(LEN=*), INTENT(IN) :: CFIN
  CHARACTER(LEN=*), INTENT(IN) :: INS
  INTEGER(I4), INTENT(OUT) :: NOB
  INTEGER(I4), DIMENSION(NMAX),INTENT(OUT) :: KTYP, PAR, PRF
  REAL   (DP), DIMENSION(NMAX),INTENT(OUT) :: TIM
  REAL   (R8), DIMENSION(NMAX),INTENT(OUT) :: DEP,  VAL, LON, LAT
  CHARACTER(LEN=8 ), DIMENSION(NMAX),INTENT(OUT) :: PLA

  INTEGER(I4), PARAMETER :: MAXVARS=14
  INTEGER(I4) :: KOBS, KPRF, NOTW, NFLR, NUSP, KKTP, KKOP
  INTEGER(I4) :: IERR, JK, J, KPAR, J2
  REAL(DP)    :: ZTIME, TIME
  REAL(R8)    :: RLAT, RLON, VALT(MAXVARS), &
  & DEPTH, TM_SAL
  INTEGER(I4) :: ICRUISE, IY, IM, ID, JJ, NLEVS, NVARS
  INTEGER(I4) :: IP2(MAXVARS),MSIG(MAXVARS),IDERROR(0:MAXVARS),IORIGFLAG(0:MAXVARS)
  LOGICAL     :: LL_INIT, LLTIME, LCOO
  CHARACTER(LEN=2) :: CC
 
  WRITE(IOUN,*) ' ^^^ READ_WODRF PROCESSING FILE ',TRIM(CFIN)

  KOBS=0
  KKOP=0
  KPRF=0
  NOTW=0
  NUSP=0
  NFLR=0

  IERR=0

  IF( INS(1:4) .EQ. 'OSDO' ) THEN
    KKTP = KKOSD
  ELSE
    CALL ABOR1('UNRECOGNIZED INPUT INSTRUMENT IN READWODFILE')
  ENDIF

800   FORMAT(1X,A2,I8,1X,F7.3,1X,F8.3,1X,I4,1X,I2,1X,I2,&
     &      1X,F7.2,1X,I8,1X,I6,1X,I6)
809   FORMAT(F7.1,1X,I1,I1,14(2X,I1,1X,F9.3,1X,I1,1X,I1,I1))

  OPEN(71,FILE=CFIN,STATUS='OLD',IOSTAT=IERR)
  IF( IERR .NE. 0 ) CALL ABOR1('READWODFILE: OPEN FAILED')
  PROFILES : DO WHILE (IERR .EQ. 0)
     LCOO = .TRUE.
     READ(71,800,IOSTAT=IERR) &
     & CC, ICRUISE, RLAT, RLON, IY, IM, ID, TIME, JJ, NLEVS, NVARS
     IF( IERR .NE. 0 ) EXIT PROFILES
     IF(RLON .LT. -180._R8 .OR. RLON .GT. 360._R8 &
     .OR. RLAT .LT. -90._R8 .OR. RLAT .GT. 180.-R8 ) THEN
         KKOP=KKOP+1
         LCOO = .FALSE.
     ENDIF
     IF( TIME < 0 ) TIME = 12._R8
     CALL YMDS2JU(IY,IM,ID,TIME/24._DP,ZTIME)
     LLTIME = (ZTIME .GE. ZJULSTART .AND. ZTIME .LE. ZJULEND)
     IF(.NOT.LLTIME) NOTW = NOTW + 1
     IF( RLON .GT. 180._R8 ) RLON = RLON -360._R8
     IF(.NOT.LLTIME .OR. .NOT.LCOO) THEN
        DO JK=1,NLEVS
           READ(71,*)
        ENDDO
     ELSE
       LL_INIT = .FALSE.
       DO JK=1,NLEVS
         READ(71,809,IOSTAT=IERR) DEPTH,IDERROR(0),IORIGFLAG(0),&
         & (IP2(J),VALT(J),MSIG(J),IDERROR(J),IORIGFLAG(J),J=1,NVARS)
         PARAM : DO J=1,NVARS
           IF( IDERROR(0)+IDERROR(J)+IORIGFLAG(0)+IORIGFLAG(J) .GT. 0 &
               .OR. VALT(J) .LT. -2._R8 .OR. VALT(J) .GT. 45._R8 ) THEN
               NFLR = NFLR + 1
               CYCLE PARAM
           ENDIF
           IF( IP2(J) .EQ. 1 ) THEN
               KPAR = KKTEMP
           ELSEIF ( IP2(J) .EQ. 2 ) THEN
               KPAR = KKSAL
           ELSE
               NUSP = NUSP + 1
               CYCLE PARAM
           ENDIF

           KOBS=KOBS+1
           IF(.NOT. LL_INIT) THEN
               LL_INIT=.TRUE.
               KPRF = KPRF + 1
           ENDIF
           LON(KOBS) = RLON
           LAT(KOBS) = RLAT
           PAR(KOBS) = KPAR
           KTYP(KOBS) = KKTP
           DEP(KOBS) = DEPTH
           VAL(KOBS) = VALT(J)
           TIM(KOBS) = ZTIME
           PRF(KOBS) = KPRF
           WRITE(PLA(KOBS),'(I8.8)') ICRUISE
           IF( KPAR .EQ. KKTEMP ) THEN
               TM_SAL = -9._R8
               DO J2=1,NVARS
                  IF( IP2(J2) .EQ. 2 .AND. VALT(J2) .GT. 0._R8 .AND. &
                & VALT(J2) .LE. 40._R8 ) TM_SAL = VALT(J2)
               ENDDO
               IF( TM_SAL .GT. 0._R8 ) THEN
                   VAL(KOBS) = POTEMP( TM_SAL, VAL(KOBS), DEP_TO_P(DEPTH,RLAT),0._R8)
               ELSE
                   VAL(KOBS) = VAL(KOBS) + 200._R8
               ENDIF 
           ENDIF
         ENDDO PARAM
       ENDDO
     ENDIF
  ENDDO PROFILES
  CLOSE(71)

  NOB = KOBS

  WRITE(IOUN,*)
  WRITE(IOUN,*) ' ^^^ REPORT '
  WRITE(IOUN,*) ' NUMBER OF BAD  PROFILES :', KKOP
  WRITE(IOUN,*) ' NUMBER OF GOOD PROFILES :', KPRF
  WRITE(IOUN,*) ' NUMBER OF GOOD OBSERV.  :', KOBS
  WRITE(IOUN,*) ' NUMBER OF P OUT OF TIME :', NOTW
  WRITE(IOUN,*) ' NUMBER OF UNSUPPOR. PAR.:', NUSP
  WRITE(IOUN,*) ' NUMBER OF BAD FLAGS OBS :', NFLR
  WRITE(IOUN,*) 

CONTAINS

   REAL(KIND=R8) FUNCTION dep_to_p( p_dep, p_lat )
      REAL(KIND=R8), INTENT(IN) :: p_dep    ! Depth in meters
      REAL(KIND=R8), INTENT(IN) :: p_lat    ! Latitude in degrees
      REAL(KIND=R8) :: z_x
      REAL(KIND=R8) :: z_c1
      REAL(KIND=R8) :: z_c2
      REAL(KIND=R8) :: z_d

      z_x = SIN( p_lat / 57.29578 )
      z_x = z_x * z_x
      z_c1 = ( 5.92  + 5.25 * z_x ) * 1e-3
      z_c2 = 2.21e-6
      z_d = ( z_c1 - 1 ) * ( z_c1 - 1  ) - 4 * z_c2 * p_dep
      dep_to_p = (( 1 - z_c1 ) - SQRT( z_d )) / ( 2 * z_c2 )
   END FUNCTION dep_to_p

   REAL(KIND=R8) FUNCTION potemp( ps, pt, pp, ppr )

      REAL(KIND=R8), INTENT(IN) :: ps
      REAL(KIND=R8), INTENT(IN) :: pt
      REAL(KIND=R8), INTENT(IN) :: pp
      REAL(KIND=R8), INTENT(IN) :: ppr

      REAL(KIND=R8) :: zpol
      REAL(KIND=R8), PARAMETER :: a1 =  1.067610e-05
      REAL(KIND=R8), PARAMETER :: a2 = -1.434297e-06
      REAL(KIND=R8), PARAMETER :: a3 = -7.566349e-09
      REAL(KIND=R8), PARAMETER :: a4 = -8.535585e-06
      REAL(KIND=R8), PARAMETER :: a5 =  3.074672e-08
      REAL(KIND=R8), PARAMETER :: a6 =  1.918639e-08
      REAL(KIND=R8), PARAMETER :: a7 =  1.788718e-10

      zpol = a1 + a2 * ps + a3 * ( pp + ppr ) + a4 * pt &
         & + a5 * ps * pt + a6 * pt * pt + a7 * pt * ( pp + ppr )

      potemp = pt + ( pp - ppr ) * zpol

   END FUNCTION potemp

END SUBROUTINE READWODFILE
