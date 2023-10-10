SUBROUTINE RCFLF_X(IM,JM,KM,FLD)
USE MYFRTPROF, ONLY : MYFRTPROF_WALL
USE RECFILTER

!-- RECURSIVE FILTER IN X-DIR

 USE SET_KND
 USE MYNETCDF

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

 CALL MYFRTPROF_WALL('RCFLF_X: RECURSIVE FILTER X-DIR',0)

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
           A(J,INXF(I,J,K)) = FLD(I,J,K2)
         ENDDO
        ENDDO

        ALP(:,:) = RFF_AEX(:,:,K)
        BTA(:,:) = RFF_BEX(:,:,K)
        IF (RCF_TYPE .EQ. 2 ) THEN
          GMA(:,:) = RFF_GEX(:,:,K)
          DEL(:,:) = RFF_DEX(:,:,K)
        ENDIF

        DO KTR = 1,NIT

! POSITIVE DIRECTION
        IF( RCF_TYPE .EQ.1 ) THEN
          IF( KTR.EQ.1 )THEN
             B(:,1) = (1._R8-ALP(:,1)) * A(:,1)
          ELSEIF( KTR.EQ.2 )THEN
             B(:,1) = A(:,1) / (1._R8+ALP(:,1))
          ELSE
             B(:,1) = (1._R8-ALP(:,1)) * (A(:,1)-ALP(:,1)**3 * A(:,2)) / &
             & (1._R8-ALP(:,1)**2)**2
          ENDIF
          DO J=2,IMXF(K)
             B(:,J) = ALP(:,J)*B(:,J-1) + (1._R8-ALP(:,J))*A(:,J)
          ENDDO
        ELSE
           B(:,-2) = 0._R8
           B(:,-1) = 0._R8
           B(:, 0) = 0._R8
           B(:, 1) = DEL(:,1)*A(:,1)
           DO J=1,IMXF(K)
             B(:,J) = ALP(:,J)*B(:,J-1) + BTA(:,J)*B(:,J-2) + &
             GMA(:,J)*B(:,J-3)+DEL(:,J)*A(:,J)
           ENDDO
        ENDIF

! NEGATIVE DIRECTION
        IF( RCF_TYPE .EQ.1 ) THEN
          IF( KTR.EQ.1 )THEN
            C(:,IMXF(K)) = B(:,IMXF(K)) / (1._R8+BTA(:,IMXF(K)))
          ELSE
            C(:,IMXF(K)) = (1._R8-BTA(:,IMXF(K))) * (B(:,IMXF(K))-BTA(:,IMXF(K))**3 * &
            & B(:,IMXF(K)-1)) / (1._R8-BTA(:,IMXF(K))**2)**2
          ENDIF
          DO J=IMXF(K)-1,1,-1
            C(:,J) = BTA(:,J)*C(:,J+1) + (1._R8-BTA(:,J))*B(:,J)
          ENDDO
        ELSE
          C(:,IMXF(K)) = DEL(:,IMXF(K))*B(:,IMXF(K))
          C(:,IMXF(K)+1 ) = 0._R8
          C(:,IMXF(K)+2 ) = 0._R8
          C(:,IMXF(K)+3 ) = 0._R8
          DO J=IMXF(K),1,-1
            C(:,J) = ALP(:,J)*C(:,J+1)+BTA(:,J)*C(:,J+2) + &
            GMA(:,J)*C(:,J+3)+DEL(:,J)*B(:,J)
          ENDDO
        ENDIF

        IF(LLDEBRCFLX) THEN
           WRITE(IOUNLOG,*) '-- ',K,KTR,SUM(BTA(:,1:IMXF(K)))/REAL(JM*IMXF(K)),&
                            & SUM(C(:,1:IMXF(K)))/REAL(JM*IMXF(K)),&
                            & SUM(B(:,1:IMXF(K)))/REAL(JM*IMXF(K))
           CALL FLUSH(IOUNLOG)
        ENDIF

        A(:,1:IMXF(K)) = C(:,1:IMXF(K))
        B(:,:) = 0._R8
        C(:,:) = 0._R8

       ENDDO

       DO J=1,JM
        DO I=1,IM
         FLD(I,J,K2) = A(J,INXF(I,J,K))
        ENDDO
       ENDDO

   ENDDO OUTER
#ifdef SHARED_MEMORY
!$OMP END DO
!$OMP END PARALLEL
#endif
CALL FLUSH(IOUNLOG)
IF(LLDEBRCFLX) LLDEBRCFLX=.FALSE.
 CALL MYFRTPROF_WALL('RCFLF_X: RECURSIVE FILTER X-DIR',1)
END SUBROUTINE RCFLF_X
