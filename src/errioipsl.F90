!$ID: ERRIOIPSL.F90 11 2007-03-12 16:01:04Z BELLIER $
!-
MODULE ERRIOIPSL
!---------------------------------------------------------------------
IMPLICIT NONE
!-
PRIVATE
!-
PUBLIC :: IPSLNLF, IPSLERR, IPSLERR_ACT, IPSLERR_INQ, HISTERR, IPSLDBG
!-
  INTEGER :: N_L=6, ILV_CUR=0, ILV_MAX=0
  LOGICAL :: IOIPSL_DEBUG=.FALSE., LACT_MODE=.TRUE.
!-
!===
CONTAINS
!===
SUBROUTINE IPSLNLF (NEW_NUMBER,OLD_NUMBER)
!!--------------------------------------------------------------------
!! THE "IPSLNLF" ROUTINE ALLOWS TO KNOW AND MODIFY
!! THE CURRENT LOGICAL NUMBER FOR THE MESSAGES.
!!
!! SUBROUTINE IPSLNLF (NEW_NUMBER,OLD_NUMBER)
!!
!! OPTIONAL INPUT ARGUMENT
!!
!! (I) NEW_NUMBER : NEW LOGICAL NUMBER OF THE FILE
!!
!! OPTIONAL OUTPUT ARGUMENT
!!
!! (I) OLD_NUMBER : CURRENT LOGICAL NUMBER OF THE FILE
!!--------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,OPTIONAL,INTENT(IN)  :: NEW_NUMBER
  INTEGER,OPTIONAL,INTENT(OUT) :: OLD_NUMBER
!---------------------------------------------------------------------
  IF (PRESENT(OLD_NUMBER)) THEN
    OLD_NUMBER = N_L
  ENDIF
  IF (PRESENT(NEW_NUMBER)) THEN
    N_L = NEW_NUMBER
  ENDIF
!---------------------
END SUBROUTINE IPSLNLF
!===
SUBROUTINE IPSLERR (PLEV,PCNAME,PSTR1,PSTR2,PSTR3)
!---------------------------------------------------------------------
!! THE "IPSLERR" ROUTINE
!! ALLOWS TO HANDLE THE MESSAGES TO THE USER.
!!
!! INPUT
!!
!! PLEV   : CATEGORY OF MESSAGE TO BE REPORTED TO THE USER
!!          1 = NOTE TO THE USER
!!          2 = WARNING TO THE USER
!!          3 = FATAL ERROR
!! PCNAME : NAME OF SUBROUTINE WHICH HAS CALLED IPSLERR
!! PSTR1
!! PSTR2  : STRINGS CONTAINING THE EXPLANATIONS TO THE USER
!! PSTR3
!---------------------------------------------------------------------
   IMPLICIT NONE
!-
   INTEGER :: PLEV
   CHARACTER(LEN=*) :: PCNAME,PSTR1,PSTR2,PSTR3
!-
   CHARACTER(LEN=30),DIMENSION(3) :: PEMSG = &
  &  (/ "NOTE TO THE USER FROM ROUTINE ", &
  &     "WARNING FROM ROUTINE          ", &
  &     "FATAL ERROR FROM ROUTINE      " /)
!---------------------------------------------------------------------
   IF ( (PLEV >= 1).AND.(PLEV <= 3) ) THEN
     ILV_CUR = PLEV
     ILV_MAX = MAX(ILV_MAX,PLEV)
     WRITE(N_L,'(/,A," ",A)') TRIM(PEMSG(PLEV)),TRIM(PCNAME)
     WRITE(N_L,'(3(" --> ",A,/))') TRIM(PSTR1),TRIM(PSTR2),TRIM(PSTR3)
   ENDIF
   IF ( (PLEV == 3).AND.LACT_MODE) THEN
     STOP 'FATAL ERROR FROM IOIPSL. SEE STDOUT FOR MORE DETAILS'
   ENDIF
!---------------------
END SUBROUTINE IPSLERR
!===
SUBROUTINE IPSLERR_ACT (NEW_MODE,OLD_MODE)
!!--------------------------------------------------------------------
!! THE "IPSLERR_ACT" ROUTINE ALLOWS TO KNOW AND MODIFY
!! THE CURRENT "ACTION MODE" FOR THE ERROR MESSAGES,
!! AND REINITIALIZE THE ERROR LEVEL VALUES.
!!
!! SUBROUTINE IPSLERR_ACT (NEW_MODE,OLD_MODE)
!!
!! OPTIONAL INPUT ARGUMENT
!!
!! (I) NEW_MODE : NEW ERROR ACTION MODE
!!                .TRUE.  -> STOP     IN CASE OF FATAL ERROR
!!                .FALSE. -> CONTINUE IN CASE OF FATAL ERROR
!!
!! OPTIONAL OUTPUT ARGUMENT
!!
!! (I) OLD_MODE : CURRENT ERROR ACTION MODE
!!--------------------------------------------------------------------
  IMPLICIT NONE
!-
  LOGICAL,OPTIONAL,INTENT(IN)  :: NEW_MODE
  LOGICAL,OPTIONAL,INTENT(OUT) :: OLD_MODE
!---------------------------------------------------------------------
  IF (PRESENT(OLD_MODE)) THEN
    OLD_MODE = LACT_MODE
  ENDIF
  IF (PRESENT(NEW_MODE)) THEN
    LACT_MODE = NEW_MODE
  ENDIF
  ILV_CUR = 0
  ILV_MAX = 0
!-------------------------
END SUBROUTINE IPSLERR_ACT
!===
SUBROUTINE IPSLERR_INQ (CURRENT_LEVEL,MAXIMUM_LEVEL)
!!--------------------------------------------------------------------
!! THE "IPSLERR_INQ" ROUTINE ALLOWS TO KNOW
!! THE CURRENT LEVEL OF THE ERROR MESSAGES
!! AND THE MAXIMUM LEVEL ENCOUNTERED SINCE THE
!! LAST CALL TO "IPSLERR_ACT".
!!
!! SUBROUTINE IPSLERR_INQ (CURRENT_LEVEL,MAXIMUM_LEVEL)
!!
!! OPTIONAL OUTPUT ARGUMENT
!!
!! (I) CURRENT_LEVEL : CURRENT ERROR LEVEL
!! (I) MAXIMUM_LEVEL : MAXIMUM ERROR LEVEL
!!--------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,OPTIONAL,INTENT(OUT) :: CURRENT_LEVEL,MAXIMUM_LEVEL
!---------------------------------------------------------------------
  IF (PRESENT(CURRENT_LEVEL)) THEN
    CURRENT_LEVEL = ILV_CUR
  ENDIF
  IF (PRESENT(MAXIMUM_LEVEL)) THEN
    MAXIMUM_LEVEL = ILV_MAX
  ENDIF
!-------------------------
END SUBROUTINE IPSLERR_INQ
!===
SUBROUTINE HISTERR (PLEV,PCNAME,PSTR1,PSTR2,PSTR3)
!---------------------------------------------------------------------
!- INPUT
!- PLEV   : CATEGORY OF MESSAGE TO BE REPORTED TO THE USER
!-          1 = NOTE TO THE USER
!-          2 = WARNING TO THE USER
!-          3 = FATAL ERROR
!- PCNAME : NAME OF SUBROUTINE WHICH HAS CALLED HISTERR
!- PSTR1
!- PSTR2  : STRING CONTAINING THE EXPLANATIONS TO THE USER
!- PSTR3
!---------------------------------------------------------------------
   IMPLICIT NONE
!-
   INTEGER :: PLEV
   CHARACTER(LEN=*) :: PCNAME,PSTR1,PSTR2,PSTR3
!-
   CHARACTER(LEN=30),DIMENSION(3) :: PEMSG = &
  &  (/ "NOTE TO THE USER FROM ROUTINE ", &
  &     "WARNING FROM ROUTINE          ", &
  &     "FATAL ERROR FROM ROUTINE      " /)
!---------------------------------------------------------------------
   IF ( (PLEV >= 1).AND.(PLEV <= 3) ) THEN
     WRITE(*,'("     ")')
     WRITE(*,'(A," ",A)') TRIM(PEMSG(PLEV)),TRIM(PCNAME)
     WRITE(*,'(" --> ",A)') PSTR1
     WRITE(*,'(" --> ",A)') PSTR2
     WRITE(*,'(" --> ",A)') PSTR3
   ENDIF
   IF (PLEV == 3) THEN
     STOP 'FATAL ERROR FROM IOIPSL. SEE STDOUT FOR MORE DETAILS'
   ENDIF
!---------------------
END SUBROUTINE HISTERR
!===
SUBROUTINE IPSLDBG (NEW_STATUS,OLD_STATUS)
!!--------------------------------------------------------------------
!! THE "IPSLDBG" ROUTINE
!! ALLOWS TO ACTIVATE OR DEACTIVATE THE DEBUG,
!! AND TO KNOW THE CURRENT STATUS OF THE DEBUG.
!!
!! SUBROUTINE IPSLDBG (NEW_STATUS,OLD_STATUS)
!!
!! OPTIONAL INPUT ARGUMENT
!!
!! (L) NEW_STATUS : NEW STATUS OF THE DEBUG
!!
!! OPTIONAL OUTPUT ARGUMENT
!!
!! (L) OLD_STATUS : CURRENT STATUS OF THE DEBUG
!!--------------------------------------------------------------------
  IMPLICIT NONE
!-
  LOGICAL,OPTIONAL,INTENT(IN)  :: NEW_STATUS
  LOGICAL,OPTIONAL,INTENT(OUT) :: OLD_STATUS
!---------------------------------------------------------------------
  IF (PRESENT(OLD_STATUS)) THEN
    OLD_STATUS = IOIPSL_DEBUG
  ENDIF
  IF (PRESENT(NEW_STATUS)) THEN
    IOIPSL_DEBUG = NEW_STATUS
  ENDIF
!---------------------
END SUBROUTINE IPSLDBG
!===
!-------------------
END MODULE ERRIOIPSL
