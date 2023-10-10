SUBROUTINE RCFL_Y_2D(IM,JM,FLD)
USE RECFILTER

 USE SET_KND

 IMPLICIT NONE

 INTEGER(I4)    :: IM, JM

 REAL(R8)       :: FLD(IM,JM)
 REAL(R8)       :: A(IM,JMAX), B(IM,JMAX), C(IM,JMAX)
 REAL(R8)       :: ALP(IM,JMAX), BTA(IM,JMAX)

 INTEGER(I4)    :: I,J,KTR
 INTEGER(I4)    :: NIT

 RFKCOUNT(6) = RFKCOUNT(6)+1
 IF(LLPPT) WRITE(IOUNOUT,*) 'RCFL_Y_2D ',RFKCOUNT(6)

 NIT=RF_NTR

        A(:,:) = 0.0
        B(:,:) = 0.0
        C(:,:) = 0.0

        DO J=1,JM
         DO I=1,IM
            A(I,JNX(I,J,1)) = FLD(I,J)
         ENDDO
        ENDDO

       ALP(:,:) = RF_AEY(:,:,1)
       BTA(:,:) = RF_BEY(:,:,1)

       DO KTR = 1,NIT

! POSITIVE DIRECTION
         IF( KTR.EQ.1 )THEN
           B(:,1) = (1.-ALP(:,1)) * A(:,1)
         ELSEIF( KTR.EQ.2 )THEN
            B(:,1) = A(:,1) / (1.+ALP(:,1))
         ELSE
            B(:,1) = (1.-ALP(:,1)) * (A(:,1)-ALP(:,1)**3 * A(:,2)) /&
            &  (1.-ALP(:,1)**2)**2
         ENDIF

        DO J=2,JMX(1)
             B(:,J) = ALP(:,J)*B(:,J-1) + (1.-ALP(:,J))*A(:,J)
        ENDDO
! NEGATIVE DIRECTION
         IF( KTR.EQ.1 )THEN
           C(:,JMX(1)) = B(:,JMX(1)) / (1.+BTA(:,JMX(1)))
         ELSE
           C(:,JMX(1)) = (1.-BTA(:,JMX(1))) * (B(:,JMX(1))-BTA(:,JMX(1))**3 *&
           &  B(:,JMX(1)-1)) / (1.-BTA(:,JMX(1))**2)**2
         ENDIF

         DO J=JMX(1)-1,1,-1
          C(:,J) = BTA(:,J)*C(:,J+1) + (1.-BTA(:,J))*B(:,J)
         ENDDO

         A(:,:) = C(:,:)

       ENDDO

        DO J=1,JM
         DO I=1,IM
          FLD(I,J) = A(I,JNX(I,J,1))
         ENDDO
        ENDDO

END SUBROUTINE RCFL_Y_2D
