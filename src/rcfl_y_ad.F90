SUBROUTINE RCFL_Y_AD(IM,JM,KM,FLD)

 USE EXTGRID
 USE RECFILTER
 USE SET_KND
 USE MYFRTPROF

 IMPLICIT NONE

 INTEGER(I4)    :: IM, JM, KM

 REAL(R8)       :: FLD(IM,JM,KM)
 REAL(R8)       :: A(IM,JMAX)
 REAL(R8)       :: B(IM,-2:JMAX)
 REAL(R8)       :: C(IM,JMAX+3)
 REAL(R8)       :: ALP(IM,JMAX), BTA(IM,JMAX)
 REAL(R8)       :: GMA(IM,JMAX), DEL(IM,JMAX)

 INTEGER(I4)    :: I,J,K, KTR,K2
 INTEGER(I4)    :: NIT

CALL MYFRTPROF_WALL('RCFL_Y_AD: ADJOINT RECURSIVE FILTER Y-DIR',0)

 NIT=RF_NTR

#ifdef SHARED_MEMORY
!$OMP PARALLEL DEFAULT(SHARED),PRIVATE(K2,K,A,B,C,I,J,KTR,ALP,BTA,GMA,DEL)
!$OMP DO SCHEDULE(DYNAMIC,1)
#endif
OUTER :   DO K2=1,KM
        K=K2

        A(:,:) = 0._R8
        B(:,:) = 0._R8
        C(:,:) = 0._R8

        DO J=1,JM
         DO I=1,IM
            C(I,JNX(I,J,K)) = FLD(I,J,K2)
         ENDDO
        ENDDO

        ALP(:,:) = RF_AEY(:,:,K)
        BTA(:,:) = RF_BEY(:,:,K)
        IF( RCF_TYPE .EQ. 2 ) THEN
          GMA(:,:) = RF_GEY(:,:,K)
          DEL(:,:) = RF_DEY(:,:,K)
        ENDIF

       DO KTR = 1,NIT

         IF( RCF_TYPE .EQ. 1 ) THEN
           DO J=1,JMX(K)-1
             C(:,J+1) = C(:,J+1) + BTA(:,J)*C(:,J)
           ENDDO
           DO J=1,JMX(K)-1
             B(:,J)   = (1._R8-BTA(:,J))*C(:,J)
           ENDDO
           IF( KTR .EQ. 1 ) THEN
              B(:,JMX(K)) = B(:,JMX(K)) + C(:,JMX(K)) / (1._R8+BTA(:,JMX(K)))
           ELSE
              B(:,JMX(K)  ) = B(:,JMX(K)  ) + (1._R8-BTA(:,JMX(K))) * C(:,JMX(K)) /&
              &  (1._R8-BTA(:,JMX(K))**2)**2
              B(:,JMX(K)-1) = B(:,JMX(K)-1) - (1._R8-BTA(:,JMX(K))) * BTA(:,JMX(K))**3&
              &  * C(:,JMX(K)) / (1._R8-BTA(:,JMX(K))**2)**2
           ENDIF
         ELSE
           DO J=1,JMX(K)
             C(:,J+1) = C(:,J+1) + ALP(:,J)*C(:,J)
             C(:,J+2) = C(:,J+2) + BTA(:,J)*C(:,J)
             C(:,J+3) = C(:,J+3) + GMA(:,J)*C(:,J)
             B(:,J)   = B(:,J)   + DEL(:,J)*C(:,J)
           ENDDO
           B(:,JMX(K)-2) = B(:,JMX(K)-2) + DEL(:,JMX(K)-2)*C(:,JMX(K)-2)
           C(:,JMX(K)-1) = C(:,JMX(K)-1) + ALP(:,JMX(K)-2)*C(:,JMX(K)-2)
           C(:,JMX(K)  ) = C(:,JMX(K)  ) + BTA(:,JMX(K)-2)*C(:,JMX(K)-2)

           B(:,JMX(K)-1) = B(:,JMX(K)-1) + DEL(:,JMX(K)-1)*C(:,JMX(K)-1)
           C(:,JMX(K)  ) = C(:,JMX(K)  ) + ALP(:,JMX(K)-1)*C(:,JMX(K)-1)

           B(:,JMX(K)) = B(:,JMX(K)) + DEL(:,JMX(K))*C(:,JMX(K))
         ENDIF

        IF( RCF_TYPE .EQ. 1 ) THEN
          DO J=JMX(K),2,-1
            B(:,J-1) = B(:,J-1) + ALP(:,J)*B(:,J)
          ENDDO
          DO J=JMX(K),2,-1
            A(:,J) = A(:,J) + (1._R8-ALP(:,J))*B(:,J)
          ENDDO
          IF( KTR .EQ. 1 ) THEN
            A(:,1) = A(:,1) + (1._R8-ALP(:,1)) * B(:,1)
          ELSEIF( KTR .EQ. 2 ) THEN
            A(:,1) = A(:,1) +  B(:,1) / (1._R8+ALP(:,1))
          ELSE
            A(:,1) = A(:,1) + (1._R8-ALP(:,1)) * B(:,1) / (1._R8-ALP(:,1)**2)**2
            A(:,2) = A(:,2) - (1._R8-ALP(:,1)) * ALP(:,1)**3 * B(:,1) /&
            &  (1._R8-ALP(:,1)**2)**2
          ENDIF
        ELSE
          DO J=JMX(K),1,-1
            B(:,J-1) = B(:,J-1) + ALP(:,J)*B(:,J)
            B(:,J-2) = B(:,J-2) + BTA(:,J)*B(:,J)
            B(:,J-3) = B(:,J-3) + GMA(:,J)*B(:,J)
            A(:,J  ) = A(:,J  ) + DEL(:,J)*B(:,J)
          ENDDO
          A(:,3) = A(:,3) + DEL(:,3)*B(:,3)
          B(:,2) = B(:,2) + ALP(:,3)*B(:,3)
          B(:,1) = B(:,1) + BTA(:,3)*B(:,3)

          A(:,2) = A(:,2) + DEL(:,2)*B(:,2)
          B(:,1) = B(:,1) + ALP(:,2)*B(:,2)

          A(:,1) = A(:,1) + DEL(:,1)*B(:,1)
        ENDIF

        C(:,1:JMAX) = A(:,1:JMAX)
        B(:,:) = 0._R8
        A(:,:) = 0._R8

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
