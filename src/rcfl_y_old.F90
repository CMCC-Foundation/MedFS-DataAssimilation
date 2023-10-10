SUBROUTINE RCFL_Y_OLD(IM,JM,KM,L_JMAX,L_AL,L_BT,FLD,L_JNX,L_JMX)
USE RECFILTER

 USE SET_KND

 IMPLICIT NONE

 INTEGER(I4)    :: NTR
 INTEGER(I4)    :: IM, JM, KM, L_JMAX

 REAL(R8)       :: FLD(IM,JM)
 REAL(R8)       :: L_AL(IM,L_JMAX), L_BT(IM,L_JMAX)
 INTEGER(I4)    :: L_JNX(IM,JM), L_JMX
 REAL(R8)       :: A(IM,L_JMAX), B(IM,L_JMAX), C(IM,L_JMAX)
 REAL(R8)       :: L_ALP(IM,L_JMAX), L_BTA(IM,L_JMAX)

 INTEGER(I4)    :: I,J,K, KTR

 NTR=RF_NTR

        A(:,:) = 0.0
        B(:,:) = 0.0
        C(:,:) = 0.0

        DO J=1,JM
         DO I=1,IM
            A(I,L_JNX(I,J)) = FLD(I,J)
         ENDDO
        ENDDO
         L_ALP(:,:) = L_AL(:,:)
         L_BTA(:,:) = L_BT(:,:)

       DO KTR = 1,NTR

! POSITIVE DIRECTION
         IF( KTR.EQ.1 )THEN
           B(:,1) = (1.-L_ALP(:,1)) * A(:,1)
         ELSEIF( KTR.EQ.2 )THEN
            B(:,1) = A(:,1) / (1.+L_ALP(:,1))
         ELSE
            B(:,1) = (1.-L_ALP(:,1)) * (A(:,1)-L_ALP(:,1)**3 * A(:,2)) /&
            &  (1.-L_ALP(:,1)**2)**2
         ENDIF

        DO J=2,L_JMX
             B(:,J) = L_ALP(:,J)*B(:,J-1) + (1.-L_ALP(:,J))*A(:,J)
        ENDDO

! NEGATIVE DIRECTION
         IF( KTR.EQ.1 )THEN
           C(:,L_JMX) = B(:,L_JMX) / (1.+L_BTA(:,L_JMX))
         ELSE
           C(:,L_JMX) = (1.-L_BTA(:,L_JMX)) * &
           & (B(:,L_JMX)-L_BTA(:,L_JMX)**3 *&
           &  B(:,L_JMX-1)) / (1.-L_BTA(:,L_JMX)**2)**2
         ENDIF

         DO J=L_JMX-1,1,-1
          C(:,J) = L_BTA(:,J)*C(:,J+1) + (1.-L_BTA(:,J))*B(:,J)
         ENDDO

         A(:,:) = C(:,:)

       ENDDO

        DO J=1,JM
         DO I=1,IM
          FLD(I,J) = A(I,L_JNX(I,J))
         ENDDO
        ENDDO

END SUBROUTINE RCFL_Y_OLD
