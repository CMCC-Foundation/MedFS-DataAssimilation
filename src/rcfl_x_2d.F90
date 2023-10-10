SUBROUTINE RCFL_X_2D(IM,JM,FLD)
USE RECFILTER

!-- RECURSIVE FILTER IN X-DIR

 USE SET_KND

 IMPLICIT NONE

 INTEGER(I4)    :: IM, JM

 REAL(R8)       :: FLD(IM,JM)
 REAL(R8)       :: A(JM,IMAX), B(JM,IMAX), C(JM,IMAX)
 REAL(R8)       :: ALP(JM,IMAX), BTA(JM,IMAX)

 INTEGER(I4)    :: I,J,KTR
 INTEGER(I4)    :: NIT

 RFKCOUNT(5) = RFKCOUNT(5)+1
 IF(LLPPT) WRITE(IOUNOUT,*) 'RCFL_X_2D ',RFKCOUNT(5)

 NIT=RF_NTR

        A(:,:) = 0.0
        B(:,:) = 0.0
        C(:,:) = 0.0

        DO J=1,JM
        DO I=1,IM
           A(J,INX(I,J,1)) = FLD(I,J)
        ENDDO
        ENDDO

       ALP(:,:) = RF_AEX(:,:,1)
       BTA(:,:) = RF_BEX(:,:,1)

       DO KTR = 1,NIT

! POSITIVE DIRECTION
         IF( KTR.EQ.1 )THEN
           B(:,1) = (1.-ALP(:,1)) * A(:,1)
         ELSEIF( KTR.EQ.2 )THEN
            B(:,1) = A(:,1) / (1.+ALP(:,1))
         ELSE
            B(:,1) = (1.-ALP(:,1)) * (A(:,1)-ALP(:,1)**3 * A(:,2)) / &
            & (1.-ALP(:,1)**2)**2
         ENDIF

        DO J=2,IMX(1)
             B(:,J) = ALP(:,J)*B(:,J-1) + (1.-ALP(:,J))*A(:,J)
        ENDDO

! NEGATIVE DIRECTION
         IF( KTR.EQ.1 )THEN
           C(:,IMX(1)) = B(:,IMX(1)) / (1.+BTA(:,IMX(1)))
         ELSE
           C(:,IMX(1)) = (1.-BTA(:,IMX(1))) * (B(:,IMX(1))-BTA(:,IMX(1))**3 * &
           & B(:,IMX(1)-1)) / (1.-BTA(:,IMX(1))**2)**2
         ENDIF

         DO J=IMX(1)-1,1,-1
          C(:,J) = BTA(:,J)*C(:,J+1) + (1.-BTA(:,J))*B(:,J)
         ENDDO

         A(:,:) = C(:,:)

       ENDDO

        DO J=1,JM
        DO I=1,IM
         FLD(I,J) = A(J,INX(I,J,1))
        ENDDO
        ENDDO

END SUBROUTINE RCFL_X_2D
