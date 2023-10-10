SUBROUTINE RCFLF_Y(IM,JM,KM,FLD)

USE RECFILTER
USE MYFRTPROF

 USE SET_KND

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

CALL MYFRTPROF_WALL('RCFLF_Y: RECURSIVE FILTER Y-DIR',0)

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
            A(I,JNXF(I,J,K)) = FLD(I,J,K2)
         ENDDO
        ENDDO

        ALP(:,:) = RFF_AEY(:,:,K)
        BTA(:,:) = RFF_BEY(:,:,K)
        IF( RCF_TYPE .EQ. 2 ) THEN
          GMA(:,:) = RFF_GEY(:,:,K)
          DEL(:,:) = RFF_DEY(:,:,K)
        ENDIF

       DO KTR = 1,NIT


! POSITIVE DIRECTION
        IF( RCF_TYPE .EQ. 1 ) THEN
          IF( KTR.EQ.1 )THEN
            B(:,1) = (1._R8-ALP(:,1)) * A(:,1)
          ELSEIF( KTR.EQ.2 )THEN
            B(:,1) = A(:,1) / (1._R8+ALP(:,1))
          ELSE
             B(:,1) = (1._R8-ALP(:,1)) * (A(:,1)-ALP(:,1)**3 * A(:,2)) /&
             &  (1._R8-ALP(:,1)**2)**2
          ENDIF
          DO J=2,JMXF(K)
             B(:,J) = ALP(:,J)*B(:,J-1) + (1._R8-ALP(:,J))*A(:,J)
          ENDDO
        ELSE
          B(:,-2) = 0._R8
          B(:,-1) = 0._R8
          B(:, 0) = 0._R8
          B(:, 1) = DEL(:,1)*A(:,1)
          DO J=1,JMXF(K)
             B(:,J) = ALP(:,J)*B(:,J-1) +  BTA(:,J)*B(:,J-2) + &
             GMA(:,J)*B(:,J-3) +  DEL(:,J)*A(:,J)
          ENDDO
        ENDIF

! NEGATIVE DIRECTION
        IF( RCF_TYPE .EQ. 1 ) THEN
          IF( KTR.EQ.1 )THEN
            C(:,JMXF(K)) = B(:,JMXF(K)) / (1._R8+BTA(:,JMXF(K)))
          ELSE
            C(:,JMXF(K)) = (1._R8-BTA(:,JMXF(K))) * (B(:,JMXF(K))-BTA(:,JMXF(K))**3 *&
            &  B(:,JMXF(K)-1)) / (1._R8-BTA(:,JMXF(K))**2)**2
          ENDIF
          DO J=JMXF(K)-1,1,-1
           C(:,J) = BTA(:,J)*C(:,J+1) + (1._R8-BTA(:,J))*B(:,J)
          ENDDO
        ELSE
          C(:,JMXF(K)) = DEL(:,JMXF(K)) *B(:,JMXF(K))
          C(:,JMXF(K)+1) = 0._R8
          C(:,JMXF(K)+2) = 0._R8
          C(:,JMXF(K)+3) = 0._R8
          DO J=JMXF(K),1,-1
           C(:,J) = ALP(:,J)*C(:,J+1) + BTA(:,J)*C(:,J+2) +&
           GMA(:,J)*C(:,J+3) + DEL(:,J)*B(:,J)
          ENDDO
        ENDIF

         A(:,1:JMAX) = C(:,1:JMAX)
         B(:,:) = 0._R8
         C(:,:) = 0._R8

       ENDDO

        DO J=1,JM
         DO I=1,IM
          FLD(I,J,K2) = A(I,JNXF(I,J,K))
         ENDDO
        ENDDO

   ENDDO OUTER
#ifdef SHARED_MEMORY
!$OMP END DO
!$OMP END PARALLEL
#endif

CALL MYFRTPROF_WALL('RCFLF_Y: RECURSIVE FILTER Y-DIR',1)

END SUBROUTINE RCFLF_Y
