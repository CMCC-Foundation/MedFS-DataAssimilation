   SUBROUTINE RCFL_2_AD(JPJ,G,B1,B2,B3,B)

USE SET_KND

IMPLICIT NONE

       INTEGER(I4)       C,JPJ
       REAL(R8)        W(JPJ)
       REAL(R8)        G(JPJ)
       REAL(R8)        B1,B2,B3,B
       REAL(R8)        A(JPJ)

        W(:)=0._R8
        A(:)=0._R8
        DO  C =1,JPJ-3
           G(C+1)=G(C+1)+B1*G(C)
           G(C+2)=G(C+2)+B2*G(C)
           G(C+3)=G(C+3)+B3*G(C)
           W(C  )=W(C  )+B *G(C)
        ENDDO

      W(JPJ-2)=W(JPJ-2)+B*G(JPJ-2)
      G(JPJ-1)= G(JPJ-1)+B1*G(JPJ-2)
      G(JPJ)=G(JPJ)+B2*G(JPJ-2)

      W(JPJ-1)= W(JPJ-1)+B *G(JPJ-1)
      G(JPJ  )= G(JPJ  )+B1*G(JPJ-1)

      W(JPJ)   = W(JPJ)+ B * G(JPJ)

      G(:) = 0._R8

        DO C =JPJ,4,-1
           G(C  )=G(C  )+B *W(C)
           W(C-1)=W(C-1)+B1*W(C)
           W(C-2)=W(C-2)+B2*W(C)
           W(C-3)=W(C-3)+B3*W(C)
        END DO

        G(3)=G(3)+B *W(3)
        W(2)=W(2)+B1*W(3)
        W(1)=W(1)+B2*W(3)
  
       G(2)=G(2)+ B*W(2)
       W(1)=W(1)+B1*W(2)
 
       G(1) =G(1)+B * W(1)


    END SUBROUTINE RCFL_2_AD
