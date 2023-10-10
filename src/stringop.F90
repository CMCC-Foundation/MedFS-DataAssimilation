MODULE STRINGOP
!-
!$ID: STRINGOP.F90 122 2007-08-03 13:42:20Z BELLIER $
!---------------------------------------------------------------------
!-
  INTEGER,DIMENSION(30) :: &
 & PRIME=(/1,2,3,5,7,11,13,17,19,23,29,31,37,41,43, &
 & 47,53,59,61,67,71,73,79,83,89,97,101,103,107,109/)
!-
!---------------------------------------------------------------------
CONTAINS
!=
SUBROUTINE CMPBLANK (STR)
!---------------------------------------------------------------------
!- COMPACT BLANKS
!---------------------------------------------------------------------
  CHARACTER(LEN=*),INTENT(INOUT) :: STR
!-
  INTEGER :: LCC,IPB
!---------------------------------------------------------------------
  LCC = LEN_TRIM(STR)
  IPB = 1
  DO
    IF (IPB >= LCC)   EXIT
    IF (STR(IPB:IPB+1) == '  ') THEN
      STR(IPB+1:) = STR(IPB+2:LCC)
      LCC = LCC-1
    ELSE
      IPB = IPB+1
    ENDIF
  ENDDO
!----------------------
END SUBROUTINE CMPBLANK
!===
INTEGER FUNCTION CNTPOS (C_C,L_C,C_R,L_R)
!---------------------------------------------------------------------
!- FINDS NUMBER OF OCCURENCES OF C_R IN C_C
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=*),INTENT(IN) :: C_C
  INTEGER,INTENT(IN) :: L_C
  CHARACTER(LEN=*),INTENT(IN) :: C_R
  INTEGER,INTENT(IN) :: L_R
!-
  INTEGER :: IPOS,INDX
!---------------------------------------------------------------------
  CNTPOS = 0
  IPOS   = 1
  DO
    INDX = INDEX(C_C(IPOS:L_C),C_R(1:L_R))
    IF (INDX > 0) THEN
      CNTPOS = CNTPOS+1
      IPOS   = IPOS+INDX+L_R-1
    ELSE
      EXIT
    ENDIF
  ENDDO
!------------------
END FUNCTION CNTPOS
!===
INTEGER FUNCTION FINDPOS (C_C,L_C,C_R,L_R)
!---------------------------------------------------------------------
!- FINDS POSITION OF C_R IN C_C
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=*),INTENT(IN) :: C_C
  INTEGER,INTENT(IN) :: L_C
  CHARACTER(LEN=*),INTENT(IN) :: C_R
  INTEGER,INTENT(IN) :: L_R
!---------------------------------------------------------------------
  FINDPOS = INDEX(C_C(1:L_C),C_R(1:L_R))
  IF (FINDPOS == 0)  FINDPOS=-1
!-------------------
END FUNCTION FINDPOS
!===
SUBROUTINE FIND_STR (STR_TAB,STR,POS)
!---------------------------------------------------------------------
!- THIS SUBROUTINE LOOKS FOR A STRING IN A TABLE
!---------------------------------------------------------------------
!- INPUT
!-   STR_TAB  : TABLE  OF STRINGS
!-   STR      : TARGET WE ARE LOOKING FOR
!- OUTPUT
!-   POS      : -1 IF STR NOT FOUND, ELSE VALUE IN THE TABLE
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=*),DIMENSION(:),INTENT(IN) :: STR_TAB
  CHARACTER(LEN=*),INTENT(IN) :: STR
  INTEGER,INTENT(OUT) :: POS
!-
  INTEGER :: NB_STR,I
!---------------------------------------------------------------------
  POS = -1
  NB_STR=SIZE(STR_TAB)
  IF ( NB_STR > 0 ) THEN
    DO I=1,NB_STR
      IF ( TRIM(STR_TAB(I)) == TRIM(STR) ) THEN
        POS = I
        EXIT
      ENDIF
    ENDDO
  ENDIF
!----------------------
END SUBROUTINE FIND_STR
!===
SUBROUTINE NOCOMMA (STR)
!---------------------------------------------------------------------
!- REPLACE COMMAS WITH BLANKS
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=*) :: STR
!-
  INTEGER :: I
!---------------------------------------------------------------------
  DO I=1,LEN_TRIM(STR)
    IF (STR(I:I) == ',')   STR(I:I) = ' '
  ENDDO
!---------------------
END SUBROUTINE NOCOMMA
!===
SUBROUTINE STRLOWERCASE (STR)
!---------------------------------------------------------------------
!- CONVERTS A STRING INTO LOWERCASE
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=*) :: STR
!-
  INTEGER :: I,IC
!---------------------------------------------------------------------
  DO I=1,LEN_TRIM(STR)
    IC = IACHAR(STR(I:I))
    IF ( (IC >= 65).AND.(IC <= 90) )  STR(I:I) = ACHAR(IC+32)
  ENDDO
!--------------------------
END SUBROUTINE STRLOWERCASE
!===
SUBROUTINE STRUPPERCASE (STR)
!---------------------------------------------------------------------
!- CONVERTS A STRING INTO UPPERCASE
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=*) :: STR
!-
  INTEGER :: I,IC
!---------------------------------------------------------------------
  DO I=1,LEN_TRIM(STR)
    IC = IACHAR(STR(I:I))
    IF ( (IC >= 97).AND.(IC <= 122) )  STR(I:I) = ACHAR(IC-32)
  ENDDO
!--------------------------
END SUBROUTINE STRUPPERCASE
!===
SUBROUTINE GENSIG (STR,SIG)
!---------------------------------------------------------------------
!- GENERATE A SIGNATURE FROM THE FIRST 30 CHARACTERS OF THE STRING
!- THIS SIGNATURE IS NOT UNIQUE AND THUS WHEN ONE LOOKS FOR THE
!- ONE NEEDS TO ALSO VERIFY THE STRING.
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=*) :: STR
  INTEGER          :: SIG
!-
  INTEGER :: I
!---------------------------------------------------------------------
  SIG = 0
  DO I=1,MIN(LEN_TRIM(STR),30)
    SIG = SIG + PRIME(I)*IACHAR(STR(I:I))
  ENDDO
!--------------------
END SUBROUTINE GENSIG
!===
SUBROUTINE FIND_SIG (NB_SIG,STR_TAB,STR,SIG_TAB,SIG,POS)
!---------------------------------------------------------------------
!- FIND THE STRING SIGNATURE IN A LIST OF SIGNATURES
!---------------------------------------------------------------------
!- INPUT
!-   NB_SIG      : LENGTH OF TABLE OF SIGNATURES
!-   STR_TAB     : TABLE OF STRINGS
!-   STR         : TARGET STRING WE ARE LOOKING FOR
!-   SIG_TAB     : TABLE OF SIGNATURES
!-   SIG         : TARGET SIGNATURE WE ARE LOOKING FOR
!- OUTPUT
!-   POS         : -1 IF STR NOT FOUND, ELSE VALUE IN THE TABLE
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER :: NB_SIG
  CHARACTER(LEN=*),DIMENSION(NB_SIG) :: STR_TAB
  CHARACTER(LEN=*) :: STR
  INTEGER,DIMENSION(NB_SIG) :: SIG_TAB
  INTEGER :: SIG
!-
  INTEGER :: POS
  INTEGER,DIMENSION(NB_SIG) :: LOCZEROS
!-
  INTEGER :: IL,LEN
  INTEGER,DIMENSION(1) :: MINPOS
!---------------------------------------------------------------------
  POS = -1
  IL = LEN_TRIM(STR)
!-
  IF ( NB_SIG > 0 ) THEN
    LOCZEROS = ABS(SIG_TAB(1:NB_SIG)-SIG)
    IF ( COUNT(LOCZEROS < 1) == 1 ) THEN
      MINPOS = MINLOC(LOCZEROS)
      LEN = LEN_TRIM(STR_TAB(MINPOS(1)))
      IF (     (INDEX(STR_TAB(MINPOS(1)),STR(1:IL)) > 0) &
          .AND.(LEN == IL) ) THEN
        POS = MINPOS(1)
      ENDIF
    ELSE IF ( COUNT(LOCZEROS < 1) > 1 ) THEN
      DO WHILE (COUNT(LOCZEROS < 1) >= 1 .AND. POS < 0 )
        MINPOS = MINLOC(LOCZEROS)
        LEN = LEN_TRIM(STR_TAB(MINPOS(1)))
        IF (     (INDEX(STR_TAB(MINPOS(1)),STR(1:IL)) > 0) &
            .AND.(LEN == IL) ) THEN
          POS = MINPOS(1)
        ELSE
          LOCZEROS(MINPOS(1)) = 99999
        ENDIF
      ENDDO
    ENDIF
  ENDIF
!-----------------------
 END SUBROUTINE FIND_SIG
!===
!------------------
END MODULE STRINGOP
