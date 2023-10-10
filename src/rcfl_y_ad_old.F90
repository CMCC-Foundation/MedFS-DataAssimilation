SUBROUTINE RCFL_Y_AD_OLD(IM,JM,KM,L_JMAX,L_AL,L_BT,FLD,L_JNX,L_JMX)
USE RECFILTER

 USE SET_KND

 IMPLICIT NONE

 INTEGER(I4)    :: NTR
 INTEGER(I4)    :: IM, JM, KM, L_JMAX

 REAL(R8)       :: FLD(IM,JM)
 REAL(R8)       :: L_AL(IM,L_JMAX), L_BT(IM,L_JMAX)
 INTEGER(I4)    :: L_JNX(IM,JM),L_JMX
 REAL(R8)       :: A(IM,L_JMAX), B(IM,L_JMAX), C(IM,L_JMAX)
 REAL(R8)       :: L_ALP(IM,L_JMAX), L_BTA(IM,L_JMAX)

 INTEGER(I4)    :: I,J,K, KTR

 NTR=RF_NTR

        A(:,:) = 0.0
        B(:,:) = 0.0
        C(:,:) = 0.0

        DO J=1,JM
         DO I=1,IM
            C(I,L_JNX(I,J)) = FLD(I,J)
         ENDDO
        ENDDO
         L_ALP(:,:) = L_AL(:,:)
         L_BTA(:,:) = L_BT(:,:)

       DO KTR = 1,NTR

! NEGATIVE DIRECTION (INVERSE OF)
        B(:,:) = 0.0

         DO J=1,L_JMX-1
          C(:,J+1) = C(:,J+1) + L_BTA(:,J)*C(:,J)
          B(:,J)   = (1.-L_BTA(:,J))*C(:,J)
         ENDDO


         IF( KTR.EQ.1 )THEN
           B(:,L_JMX) = B(:,L_JMX) + C(:,L_JMX) / &
           & (1.+L_BTA(:,L_JMX))
         ELSE
           B(:,L_JMX) = B(:,L_JMX) + (1.-L_BTA(:,L_JMX)) &
           & * C(:,L_JMX) /&
           &  (1.-L_BTA(:,L_JMX)**2)**2
           B(:,L_JMX-1) = B(:,L_JMX-1) - (1.-L_BTA(:,L_JMX)) * &
           & L_BTA(:,L_JMX)**3&
           & * C(:,L_JMX) / (1.-L_BTA(:,L_JMX)**2)**2
         ENDIF

! POSITIVE DIRECTION (INVERSE OF)
        A(:,:) = 0.0

         DO J=L_JMX,2,-1
          B(:,J-1) = B(:,J-1) + L_ALP(:,J)*B(:,J)
          A(:,J) = A(:,J) + (1.-L_ALP(:,J))*B(:,J)
         ENDDO


         IF( KTR.EQ.1 )THEN
           A(:,1) = A(:,1) + (1.-L_ALP(:,1)) * B(:,1)
         ELSEIF( KTR.EQ.2 )THEN
           A(:,1) = A(:,1) + B(:,1) / (1.+L_ALP(:,1))
         ELSE
           A(:,1) = A(:,1) + (1.-L_ALP(:,1)) * B(:,1) / (1.-L_ALP(:,1)**2)**2
           A(:,2) = A(:,2) - (1.-L_ALP(:,1)) * L_ALP(:,1)**3 * B(:,1) /&
           &  (1.-L_ALP(:,1)**2)**2
         ENDIF

         C(:,:) = A(:,:)

       ENDDO

        DO J=1,JM
         DO I=1,IM
          FLD(I,J) = C(I,L_JNX(I,J))
         ENDDO
        ENDDO

END SUBROUTINE RCFL_Y_AD_OLD
