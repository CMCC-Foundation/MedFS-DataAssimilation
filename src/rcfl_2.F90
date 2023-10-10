SUBROUTINE RCFL_2(JPJ,G,B1,B2,B3,B)
USE SET_KND
IMPLICIT NONE

     INTEGER(I4)  C,JPJ
     REAL(R8)   W(JPJ)
     REAL(R8)   G(JPJ)
     REAL(R8)   B1,B2,B3,B

         W(1) = B * G(1)
         W(2) = B * G(2) + B1 * W(1)
         W(3) = B * G(3) + B1 * W(2) + B2 * W(1)
          DO  C = 4,JPJ
            W(C) = B * G(C) +B1 * W(C-1) + B2 * W(C-2) + B3 * W(C-3)
         ENDDO
     
        G(JPJ)   = B * W(JPJ)
        G(JPJ-1) = B * W(JPJ-1) + B1 * G(JPJ)
        G(JPJ-2) = B * W(JPJ-2) + B1 * G(JPJ-1) + B2 * G(JPJ)
         DO  C = JPJ-3,1,-1
              G(C) = B * W(C) +  B1 * G(C+1) + B2 * G(C+2) + B3 * G(C+3)
        ENDDO

END SUBROUTINE RCFL_2
