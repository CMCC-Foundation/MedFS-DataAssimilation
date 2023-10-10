   SUBROUTINE INDEX_SORT( PVAL, KINDX, KVALS )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE INDEX_SORT  ***
      !!
      !! ** PURPOSE : GET INDICIES FOR ASCENDING ORDER FOR A DOUBLE PRECISION ARRAY
      !!
      !! ** METHOD  : HEAPSORT
      !!
      !! ** ACTION  :
      !!
      !! REFERENCES : HTTP://EN.WIKIPEDIA.ORG/WIKI/HEAPSORT
      !!
      !! HISTORY :
      !!        !  06-05  (K. MOGENSEN)  ORIGINAL CODE
      !!        !  06-10  (A. WEAVER) CLEANING
      !!        !  09-10  (A. STORTO) ADAPTED FOR OCEANVAR
      !!----------------------------------------------------------------------

      USE SET_KND
      
      IMPLICIT NONE

      !! * ARGUMENTS
      INTEGER(I4), INTENT(IN) :: &
         & KVALS                           ! NUMBER OF VALUES
      REAL(R8), DIMENSION(KVALS), INTENT(IN) :: &
         & PVAL                            ! ARRAY TO BE SORTED
      INTEGER(I4), DIMENSION(KVALS), INTENT(INOUT) :: &
         & KINDX                           ! INDICIES FOR ORDERING

      !! * LOCAL DECLARATIONS
      INTEGER(I4) :: &
         & JI,      &
         & JJ,      &
         & JT,      &
         & JN,      &
         & JPARENT, &
         & JCHILD

      DO JI = 1, KVALS
         KINDX(JI) = JI
      END DO

      JI = KVALS/2 + 1
      JN = KVALS

      MAIN_LOOP : DO

         IF ( JI > 1 ) THEN
            JI = JI-1
            JT = KINDX(JI)
         ELSE
            JT = KINDX(JN)
            KINDX(JN) = KINDX(1)
            JN = JN-1
            IF ( JN <= 1 ) THEN
               KINDX(1) = JT
               EXIT MAIN_LOOP
            ENDIF
         ENDIF

         JPARENT = JI
         JCHILD =  2 * JI

         INNER_LOOP : DO

            IF ( JCHILD > JN ) EXIT INNER_LOOP
            IF ( JCHILD < JN ) THEN
               IF ( PVAL(KINDX(JCHILD)) < PVAL(KINDX(JCHILD+1)) ) THEN
                 JCHILD = JCHILD+1
               ENDIF
            ENDIF
            IF  ( PVAL(JT) < PVAL(KINDX(JCHILD))) THEN
               KINDX(JPARENT) = KINDX(JCHILD)
               JPARENT = JCHILD
               JCHILD  = JCHILD*2
            ELSE
               JCHILD = JN + 1
            ENDIF

         END DO INNER_LOOP

         KINDX(JPARENT) = JT

      END DO MAIN_LOOP

   END SUBROUTINE INDEX_SORT

