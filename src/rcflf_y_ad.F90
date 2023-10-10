SUBROUTINE RCFLF_Y_AD(IM,JM,KM,FLD)

 USE EXTGRID
 USE RECFILTER
 USE SET_KND
 USE MYFRTPROF

 IMPLICIT NONE

 INTEGER(I4)    :: IM, JM, KM

 REAL(R8)       :: FLD(IM,JM,KM)
 REAL(R8)       :: A(IM,JMAXF)
 REAL(R8)       :: B(IM,-2:JMAXF)
 REAL(R8)       :: C(IM,JMAXF+3)
 REAL(R8)       :: ALP(IM,JMAXF), BTA(IM,JMAXF)
 REAL(R8)       :: GMA(IM,JMAXF), DEL(IM,JMAXF)

 INTEGER(I4)    :: I,J,K, KTR,K2
 INTEGER(I4)    :: NIT

CALL MYFRTPROF_WALL('RCFLF_Y_AD: ADJOINT RECURSIVE FILTER Y-DIR',0)

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
            C(I,JNXF(I,J,K)) = FLD(I,J,K2)
         ENDDO
        ENDDO

        ALP(:,:) = RFF_AEY(:,:,K)
        BTA(:,:) = RFF_BEY(:,:,K)
        IF( RCF_TYPE .EQ. 2 ) THEN
          GMA(:,:) = RFF_GEY(:,:,K)
          DEL(:,:) = RFF_DEY(:,:,K)
        ENDIF

       DO KTR = 1,NIT

         IF( RCF_TYPE .EQ. 1 ) THEN
           DO J=1,JMXF(K)-1
             C(:,J+1) = C(:,J+1) + BTA(:,J)*C(:,J)
           ENDDO
           DO J=1,JMXF(K)-1
             B(:,J)   = (1._R8-BTA(:,J))*C(:,J)
           ENDDO
           IF( KTR .EQ. 1 ) THEN
              B(:,JMXF(K)) = B(:,JMXF(K)) + C(:,JMXF(K)) / (1._R8+BTA(:,JMXF(K)))
           ELSE
              B(:,JMXF(K)  ) = B(:,JMXF(K)  ) + (1._R8-BTA(:,JMXF(K))) * C(:,JMXF(K)) /&
              &  (1._R8-BTA(:,JMXF(K))**2)**2
              B(:,JMXF(K)-1) = B(:,JMXF(K)-1) - (1._R8-BTA(:,JMXF(K))) * BTA(:,JMXF(K))**3&
              &  * C(:,JMXF(K)) / (1._R8-BTA(:,JMXF(K))**2)**2
           ENDIF
         ELSE
           DO J=1,JMXF(K)
             C(:,J+1) = C(:,J+1) + ALP(:,J)*C(:,J)
             C(:,J+2) = C(:,J+2) + BTA(:,J)*C(:,J)
             C(:,J+3) = C(:,J+3) + GMA(:,J)*C(:,J)
             B(:,J)   = B(:,J)   + DEL(:,J)*C(:,J)
           ENDDO
           B(:,JMXF(K)-2) = B(:,JMXF(K)-2) + DEL(:,JMXF(K)-2)*C(:,JMXF(K)-2)
           C(:,JMXF(K)-1) = C(:,JMXF(K)-1) + ALP(:,JMXF(K)-2)*C(:,JMXF(K)-2)
           C(:,JMXF(K)  ) = C(:,JMXF(K)  ) + BTA(:,JMXF(K)-2)*C(:,JMXF(K)-2)

           B(:,JMXF(K)-1) = B(:,JMXF(K)-1) + DEL(:,JMXF(K)-1)*C(:,JMXF(K)-1)
           C(:,JMXF(K)  ) = C(:,JMXF(K)  ) + ALP(:,JMXF(K)-1)*C(:,JMXF(K)-1)

           B(:,JMXF(K)) = B(:,JMXF(K)) + DEL(:,JMXF(K))*C(:,JMXF(K))
         ENDIF

        IF( RCF_TYPE .EQ. 1 ) THEN
          DO J=JMXF(K),2,-1
            B(:,J-1) = B(:,J-1) + ALP(:,J)*B(:,J)
          ENDDO
          DO J=JMXF(K),2,-1
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
          DO J=JMXF(K),1,-1
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
          FLD(I,J,K2) = C(I,JNXF(I,J,K))
        ENDDO
      ENDDO

   ENDDO OUTER
#ifdef SHARED_MEMORY
!$OMP END DO
!$OMP END PARALLEL
#endif

CALL MYFRTPROF_WALL('RCFLF_Y_AD: ADJOINT RECURSIVE FILTER Y-DIR',1)
END SUBROUTINE RCFLF_Y_AD
