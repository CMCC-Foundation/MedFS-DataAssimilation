SUBROUTINE RCFL_X_AD_2D(IM,JM,FLD)
USE RECFILTER


 USE SET_KND

 IMPLICIT NONE

 INTEGER(I4)    :: IM, JM

 REAL(R8)       :: FLD(IM,JM)
 REAL(R8)       :: A(JM,IMAX), B(JM,IMAX), C(JM,IMAX)
 REAL(R8)       :: ALP(JM,IMAX), BTA(JM,IMAX)

 INTEGER(I4)    :: I,J,KTR
 INTEGER(I4)    :: NIT

 RFKCOUNT(7) = RFKCOUNT(7)+1
 IF(LLPPT) WRITE(IOUNOUT,*) 'RCFL_X_AD_2D ',RFKCOUNT(7)

 NIT=RF_NTR

        A(:,:) = 0.0
        B(:,:) = 0.0
        C(:,:) = 0.0

        DO J=1,JM
        DO I=1,IM
           C(J,INX(I,J,1)) = FLD(I,J)
        ENDDO
        ENDDO

       ALP(:,:) = RF_AEX(:,:,1)
       BTA(:,:) = RF_BEX(:,:,1)

       DO KTR = 1,NIT

! NEGATIVE DIRECTION (INVERSE OF)
        B(:,:) = 0.0

         DO J=1,IMX(1)-1
          C(:,J+1) = C(:,J+1) + BTA(:,J)*C(:,J)
          B(:,J)   = (1.-BTA(:,J))*C(:,J)
         ENDDO


         IF( KTR.EQ.1 )THEN
           B(:,IMX(1)) = B(:,IMX(1)) + C(:,IMX(1)) / (1.+BTA(:,IMX(1)))
         ELSE
           B(:,IMX(1)  ) = B(:,IMX(1)  ) + (1.-BTA(:,IMX(1))) * C(:,IMX(1)) /&
           &  (1.-BTA(:,IMX(1))**2)**2
           B(:,IMX(1)-1) = B(:,IMX(1)-1) - (1.-BTA(:,IMX(1))) * BTA(:,IMX(1))**3 &
           & * C(:,IMX(1)) / (1.-BTA(:,IMX(1))**2)**2
         ENDIF

! POSITIVE DIRECTION (INVERSE OF)
        A(:,:) = 0.0

         DO J=IMX(1),2,-1
          B(:,J-1) = B(:,J-1) + ALP(:,J)*B(:,J)
          A(:,J) = A(:,J) + (1.-ALP(:,J))*B(:,J)
         ENDDO


         IF( KTR.EQ.1 )THEN
           A(:,1) = A(:,1) + (1.-ALP(:,1)) * B(:,1)
         ELSEIF( KTR.EQ.2 )THEN
           A(:,1) = A(:,1) + B(:,1) / (1.+ALP(:,1))
         ELSE
           A(:,1) = A(:,1) + (1.-ALP(:,1)) * B(:,1) / (1.-ALP(:,1)**2)**2
           A(:,2) = A(:,2) - (1.-ALP(:,1)) * ALP(:,1)**3 * B(:,1) /&
           &  (1.-ALP(:,1)**2)**2
         ENDIF

         C(:,:) = A(:,:)

       ENDDO

        DO J=1,JM
        DO I=1,IM
         FLD(I,J) = C(J,INX(I,J,1))
        ENDDO
        ENDDO

END SUBROUTINE RCFL_X_AD_2D
