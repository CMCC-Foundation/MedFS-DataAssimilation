SUBROUTINE RCFLF_X_AD(IM,JM,KM,FLD)
USE RECFILTER
USE MYFRTPROF


 USE SET_KND

 IMPLICIT NONE

 INTEGER(I4)    :: IM, JM, KM

 REAL(R8)       :: FLD(IM,JM,KM)
 REAL(R8)       :: A(JM,IMAXF)
 REAL(R8)       :: B(JM,-2:IMAXF)
 REAL(R8)       :: C(JM,IMAXF+3)
 REAL(R8)       :: ALP(JM,IMAXF), BTA(JM,IMAXF)
 REAL(R8)       :: GMA(JM,IMAXF), DEL(JM,IMAXF)

 INTEGER(I4)    :: I,J,K, KTR,K2
 INTEGER(I4)    :: NIT

CALL MYFRTPROF_WALL('RCFLF_X_AD: ADJOINT RECURSIVE FILTER X-DIR',0)

 NIT=RF_NTR

#ifdef SHARED_MEMORY
!$OMP PARALLEL DEFAULT(SHARED),PRIVATE(K2,K,A,B,C,I,J,KTR,ALP,BTA,GMA,DEL)
!$OMP DO SCHEDULE(DYNAMIC,1)
#endif
OUTER : DO K2=1,KM

        K=K2

        A(:,:) = 0._R8
        B(:,:) = 0._R8
        C(:,:) = 0._R8

        DO J=1,JM
         DO I=1,IM
           C(J,INXF(I,J,K)) = FLD(I,J,K2)
         ENDDO
        ENDDO

        ALP(:,:) = RFF_AEX(:,:,K)
        BTA(:,:) = RFF_BEX(:,:,K)
        IF( RCF_TYPE .EQ. 2 ) THEN
          GMA(:,:) = RFF_GEX(:,:,K)
          DEL(:,:) = RFF_DEX(:,:,K)
        ENDIF

        DO KTR=1,NIT

        IF( RCF_TYPE .EQ. 1 ) THEN
           DO J=1,IMXF(K)-1
             C(:,J+1) = C(:,J+1) + BTA(:,J)*C(:,J)
           ENDDO
           DO J=1,IMXF(K)-1
             B(:,J)   = (1._R8-BTA(:,J))*C(:,J)
           ENDDO
           IF( KTR .EQ. 1 ) THEN
               B(:,IMXF(K)) = B(:,IMXF(K)) + C(:,IMXF(K)) / (1._R8+BTA(:,IMXF(K)))
           ELSE
               B(:,IMXF(K)  ) = B(:,IMXF(K)) + (1._R8-BTA(:,IMXF(K))) * C(:,IMXF(K)) /&
               &  (1._R8-BTA(:,IMXF(K))**2)**2
               B(:,IMXF(K)-1) = B(:,IMXF(K)-1) - (1._R8-BTA(:,IMXF(K))) * BTA(:,IMXF(K))**3 &
               & * C(:,IMXF(K)) / (1._R8-BTA(:,IMXF(K))**2)**2
           ENDIF
        ELSE
           DO J=1,IMXF(K)
             C(:,J+1) = C(:,J+1) + ALP(:,J)*C(:,J)
             C(:,J+2) = C(:,J+2) + BTA(:,J)*C(:,J)
             C(:,J+3) = C(:,J+3) + GMA(:,J)*C(:,J)
             B(:,J)   = B(:,J)   + DEL(:,J)*C(:,J)
           ENDDO
           B(:,IMXF(K)-2) = B(:,IMXF(K)-2) + DEL(:,IMXF(K)-2)*C(:,IMXF(K)-2)
           C(:,IMXF(K)-1) = C(:,IMXF(K)-1) + ALP(:,IMXF(K)-2)*C(:,IMXF(K)-2)
           C(:,IMXF(K)  ) = C(:,IMXF(K)  ) + BTA(:,IMXF(K)-2)*C(:,IMXF(K)-2)

           B(:,IMXF(K)-1) = B(:,IMXF(K)-1) + DEL(:,IMXF(K)-1)*C(:,IMXF(K)-1)
           C(:,IMXF(K)  ) = C(:,IMXF(K)  ) + ALP(:,IMXF(K)-1)*C(:,IMXF(K)-1)

           B(:,IMXF(K)) = B(:,IMXF(K)) + DEL(:,IMXF(K))*C(:,IMXF(K))
        ENDIF

        IF( RCF_TYPE .EQ. 1 ) THEN
          DO J=IMXF(K),2,-1
            B(:,J-1) = B(:,J-1) + ALP(:,J)*B(:,J)
          ENDDO
          DO J=IMXF(K),2,-1
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
          DO J=IMXF(K),1,-1
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

        C(:,1:IMAX) = A(:,1:IMAX)
        B(:,:) = 0._R8
        A(:,:) = 0._R8

       ENDDO

       DO J=1,JM
         DO I=1,IM
           FLD(I,J,K2) = C(J,INXF(I,J,K))
         ENDDO
       ENDDO

   ENDDO OUTER
#ifdef SHARED_MEMORY
!$OMP END DO
!$OMP END PARALLEL
#endif

CALL MYFRTPROF_WALL('RCFLF_X_AD: ADJOINT RECURSIVE FILTER X-DIR',1)

END SUBROUTINE RCFLF_X_AD
