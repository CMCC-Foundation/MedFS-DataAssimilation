SUBROUTINE RCFL_Y_AD(IM,JM,KM,FLD)

 USE EXTGRID
 USE RECFILTER
 USE SET_KND
 USE MYFRTPROF

 IMPLICIT NONE

 INTEGER(I4)    :: IM, JM, KM

 REAL(R8)       :: FLD(IM,JM,2*KM)
 REAL(R8)       :: A(IM,JMAX), B(IM,JMAX), C(IM,JMAX)
 REAL(R8)       :: ALP(IM,JMAX), BTA(IM,JMAX)

 INTEGER(I4)    :: I,J,K, KTR,K2
 INTEGER(I4)    :: NIT

CALL MYFRTPROF_WALL('RCFL_Y_AD: ADJOINT RECURSIVE FILTER Y-DIR',0)

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
            C(I,JNX(I,J,K)) = FLD(I,J,K2)
         ENDDO
        ENDDO

       IF(LLEXTGRID .AND. NNSTEP_EXT .EQ. 2) THEN
         ALP(FRST:LST,:) = RF_AEY(2:OLD_IM-1,:,K)
         ALP(LST2:NEW_IM,:) = RF_AEY(2:NWRAPS+1,:,K)
         ALP(1:FRST2,:) = RF_AEY(OLD_IM-NWRAPS:OLD_IM,:,K)
         BTA(FRST:LST,:) = RF_BEY(2:OLD_IM-1,:,K)
         BTA(LST2:NEW_IM,:) = RF_BEY(2:NWRAPS+1,:,K)
         BTA(1:FRST2,:) = RF_BEY(OLD_IM-NWRAPS:OLD_IM,:,K)
       ELSE
         ALP(:,:) = RF_AEY(:,:,K)
         BTA(:,:) = RF_BEY(:,:,K)
       ENDIF

         DO J=1,JMX(K)-1
          C(:,J+1) = C(:,J+1) + BTA(:,J)*C(:,J)
          B(:,J)   = (1.-BTA(:,J))*C(:,J)
         ENDDO
         B(:,JMX(K)) = B(:,JMX(K)) + C(:,JMX(K)) / (1.+BTA(:,JMX(K)))

         DO J=JMX(K),2,-1
          B(:,J-1) = B(:,J-1) + ALP(:,J)*B(:,J)
          A(:,J) = A(:,J) + (1.-ALP(:,J))*B(:,J)
         ENDDO
         A(:,1) = A(:,1) + (1.-ALP(:,1)) * B(:,1)
         C(:,:) = A(:,:)

         B(:,:) = 0.0
         DO J=1,JMX(K)-1
          C(:,J+1) = C(:,J+1) + BTA(:,J)*C(:,J)
          B(:,J)   = (1.-BTA(:,J))*C(:,J)
         ENDDO
         B(:,JMX(K)  ) = B(:,JMX(K)  ) + (1.-BTA(:,JMX(K))) * C(:,JMX(K)) /&
         &  (1.-BTA(:,JMX(K))**2)**2
         B(:,JMX(K)-1) = B(:,JMX(K)-1) - (1.-BTA(:,JMX(K))) * BTA(:,JMX(K))**3&
         &  * C(:,JMX(K)) / (1.-BTA(:,JMX(K))**2)**2

         A(:,:) = 0.0
         DO J=JMX(K),2,-1
          B(:,J-1) = B(:,J-1) + ALP(:,J)*B(:,J)
          A(:,J) = A(:,J) + (1.-ALP(:,J))*B(:,J)
         ENDDO
         A(:,1) = A(:,1) + B(:,1) / (1.+ALP(:,1))
         C(:,:) = A(:,:)


       DO KTR = 3,NIT
        B(:,:) = 0.0

         DO J=1,JMX(K)-1
          C(:,J+1) = C(:,J+1) + BTA(:,J)*C(:,J)
          B(:,J)   = (1.-BTA(:,J))*C(:,J)
         ENDDO


           B(:,JMX(K)  ) = B(:,JMX(K)  ) + (1.-BTA(:,JMX(K))) * C(:,JMX(K)) /&
           &  (1.-BTA(:,JMX(K))**2)**2
           B(:,JMX(K)-1) = B(:,JMX(K)-1) - (1.-BTA(:,JMX(K))) * BTA(:,JMX(K))**3&
           &  * C(:,JMX(K)) / (1.-BTA(:,JMX(K))**2)**2

! POSITIVE DIRECTION (INVERSE OF)
        A(:,:) = 0.0

         DO J=JMX(K),2,-1
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
          FLD(I,J,K2) = C(I,JNX(I,J,K))
         ENDDO
        ENDDO

   ENDDO OUTER
#ifdef SHARED_MEMORY
!$OMP END DO
!$OMP END PARALLEL
#endif

CALL MYFRTPROF_WALL('RCFL_Y_AD: ADJOINT RECURSIVE FILTER Y-DIR',1)
END SUBROUTINE RCFL_Y_AD
