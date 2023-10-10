SUBROUTINE RCFL_X_AD(IM,JM,KM,FLD)
USE RECFILTER
USE MYFRTPROF


 USE SET_KND

 IMPLICIT NONE

 INTEGER(I4)    :: IM, JM, KM

 REAL(R8)       :: FLD(IM,JM,2*KM)
 REAL(R8)       :: A(JM,IMAX), B(JM,IMAX), C(JM,IMAX)
 REAL(R8)       :: ALP(JM,IMAX), BTA(JM,IMAX)

 INTEGER(I4)    :: I,J,K, KTR,K2
 INTEGER(I4)    :: NIT

CALL MYFRTPROF_WALL('RCFL_X_AD: ADJOINT RECURSIVE FILTER X-DIR',0)

 NIT=RF_NTR

#ifdef SHARED_MEMORY
!$OMP PARALLEL DEFAULT(SHARED),PRIVATE(K2,K,A,B,C,I,J,KTR,ALP,BTA)
!$OMP DO SCHEDULE(DYNAMIC,1)
#endif
OUTER :   DO K2=1,2*KM

        K=K2
        IF(K>KM.AND..NOT.LLTSAPART.AND..NOT.LLRFSR) K=K-KM

        A(:,:) = 0.0
        B(:,:) = 0.0
        C(:,:) = 0.0

        DO J=1,JM
        DO I=1,IM
           C(J,INX(I,J,K)) = FLD(I,J,K2)
        ENDDO
        ENDDO

         ALP(:,:) = RF_AEX(:,:,K)
         BTA(:,:) = RF_BEX(:,:,K)

         DO J=1,IMX(K)-1
          C(:,J+1) = C(:,J+1) + BTA(:,J)*C(:,J)
          B(:,J)   = (1.-BTA(:,J))*C(:,J)
         ENDDO
           B(:,IMX(K)) = B(:,IMX(K)) + C(:,IMX(K)) / (1.+BTA(:,IMX(K)))

         DO J=IMX(K),2,-1
          B(:,J-1) = B(:,J-1) + ALP(:,J)*B(:,J)
          A(:,J) = A(:,J) + (1.-ALP(:,J))*B(:,J)
         ENDDO
           A(:,1) = A(:,1) + (1.-ALP(:,1)) * B(:,1)
         C(:,:) = A(:,:)

        B(:,:) = 0.0
         DO J=1,IMX(K)-1
          C(:,J+1) = C(:,J+1) + BTA(:,J)*C(:,J)
          B(:,J)   = (1.-BTA(:,J))*C(:,J)
         ENDDO
           B(:,IMX(K)  ) = B(:,IMX(K)  ) + (1.-BTA(:,IMX(K))) * C(:,IMX(K)) /&
           &  (1.-BTA(:,IMX(K))**2)**2
           B(:,IMX(K)-1) = B(:,IMX(K)-1) - (1.-BTA(:,IMX(K))) * BTA(:,IMX(K))**3 &
           & * C(:,IMX(K)) / (1.-BTA(:,IMX(K))**2)**2

        A(:,:) = 0.0
         DO J=IMX(K),2,-1
          B(:,J-1) = B(:,J-1) + ALP(:,J)*B(:,J)
          A(:,J) = A(:,J) + (1.-ALP(:,J))*B(:,J)
         ENDDO
           A(:,1) = A(:,1) + B(:,1) / (1.+ALP(:,1))
         C(:,:) = A(:,:)


       DO KTR = 3,NIT

        B(:,:) = 0.0
         DO J=1,IMX(K)-1
          C(:,J+1) = C(:,J+1) + BTA(:,J)*C(:,J)
          B(:,J)   = (1.-BTA(:,J))*C(:,J)
         ENDDO

           B(:,IMX(K)  ) = B(:,IMX(K)  ) + (1.-BTA(:,IMX(K))) * C(:,IMX(K)) /&
           &  (1.-BTA(:,IMX(K))**2)**2
           B(:,IMX(K)-1) = B(:,IMX(K)-1) - (1.-BTA(:,IMX(K))) * BTA(:,IMX(K))**3 &
           & * C(:,IMX(K)) / (1.-BTA(:,IMX(K))**2)**2

        A(:,:) = 0.0
         DO J=IMX(K),2,-1
          B(:,J-1) = B(:,J-1) + ALP(:,J)*B(:,J)
          A(:,J) = A(:,J) + (1.-ALP(:,J))*B(:,J)
         ENDDO

           A(:,1) = A(:,1) + (1.-ALP(:,1)) * B(:,1) / (1.-ALP(:,1)**2)**2
           A(:,2) = A(:,2) - (1.-ALP(:,1)) * ALP(:,1)**3 * B(:,1) /&
           &  (1.-ALP(:,1)**2)**2

         C(:,:) = A(:,:)

       ENDDO

        DO J=1,JM
        DO I=1,IM
         FLD(I,J,K2) = C(J,INX(I,J,K))
        ENDDO
        ENDDO

   ENDDO OUTER
#ifdef SHARED_MEMORY
!$OMP END DO
!$OMP END PARALLEL
#endif

CALL MYFRTPROF_WALL('RCFL_X_AD: ADJOINT RECURSIVE FILTER X-DIR',1)

END SUBROUTINE RCFL_X_AD
