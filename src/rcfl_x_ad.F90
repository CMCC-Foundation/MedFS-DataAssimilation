SUBROUTINE RCFL_X_AD(IM,JM,KM,FLD)
USE RECFILTER
USE MYFRTPROF


 USE SET_KND

 IMPLICIT NONE

 INTEGER(I4)    :: IM, JM, KM

 REAL(R8)       :: FLD(IM,JM,KM)
 REAL(R8)       :: A(JM,IMAX)
 REAL(R8)       :: B(JM,-2:IMAX)
 REAL(R8)       :: C(JM,IMAX+3)
 REAL(R8)       :: ALP(JM,IMAX), BTA(JM,IMAX)
 REAL(R8)       :: GMA(JM,IMAX), DEL(JM,IMAX)

 INTEGER(I4)    :: I,J,K, KTR,K2
 INTEGER(I4)    :: NIT

CALL MYFRTPROF_WALL('RCFL_X_AD: ADJOINT RECURSIVE FILTER X-DIR',0)

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
           C(J,INX(I,J,K)) = FLD(I,J,K2)
         ENDDO
        ENDDO

        ALP(:,:) = RF_AEX(:,:,K)
        BTA(:,:) = RF_BEX(:,:,K)
        IF( RCF_TYPE .EQ. 2 ) THEN
          GMA(:,:) = RF_GEX(:,:,K)
          DEL(:,:) = RF_DEX(:,:,K)
        ENDIF

        DO KTR=1,NIT

        IF( RCF_TYPE .EQ. 1 ) THEN
           DO J=1,IMX(K)-1
             C(:,J+1) = C(:,J+1) + BTA(:,J)*C(:,J)
           ENDDO
           DO J=1,IMX(K)-1
             B(:,J)   = (1._R8-BTA(:,J))*C(:,J)
           ENDDO
           IF( KTR .EQ. 1 ) THEN
               B(:,IMX(K)) = B(:,IMX(K)) + C(:,IMX(K)) / (1._R8+BTA(:,IMX(K)))
           ELSE
               B(:,IMX(K)  ) = B(:,IMX(K)) + (1._R8-BTA(:,IMX(K))) * C(:,IMX(K)) /&
               &  (1._R8-BTA(:,IMX(K))**2)**2
               B(:,IMX(K)-1) = B(:,IMX(K)-1) - (1._R8-BTA(:,IMX(K))) * BTA(:,IMX(K))**3 &
               & * C(:,IMX(K)) / (1._R8-BTA(:,IMX(K))**2)**2
           ENDIF
        ELSE
           DO J=1,IMX(K)
             C(:,J+1) = C(:,J+1) + ALP(:,J)*C(:,J)
             C(:,J+2) = C(:,J+2) + BTA(:,J)*C(:,J)
             C(:,J+3) = C(:,J+3) + GMA(:,J)*C(:,J)
             B(:,J)   = B(:,J)   + DEL(:,J)*C(:,J)
           ENDDO
           B(:,IMX(K)-2) = B(:,IMX(K)-2) + DEL(:,IMX(K)-2)*C(:,IMX(K)-2)
           C(:,IMX(K)-1) = C(:,IMX(K)-1) + ALP(:,IMX(K)-2)*C(:,IMX(K)-2)
           C(:,IMX(K)  ) = C(:,IMX(K)  ) + BTA(:,IMX(K)-2)*C(:,IMX(K)-2)

           B(:,IMX(K)-1) = B(:,IMX(K)-1) + DEL(:,IMX(K)-1)*C(:,IMX(K)-1)
           C(:,IMX(K)  ) = C(:,IMX(K)  ) + ALP(:,IMX(K)-1)*C(:,IMX(K)-1)

           B(:,IMX(K)) = B(:,IMX(K)) + DEL(:,IMX(K))*C(:,IMX(K))
        ENDIF

        IF( RCF_TYPE .EQ. 1 ) THEN
          DO J=IMX(K),2,-1
            B(:,J-1) = B(:,J-1) + ALP(:,J)*B(:,J)
          ENDDO
          DO J=IMX(K),2,-1
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
          DO J=IMX(K),1,-1
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
