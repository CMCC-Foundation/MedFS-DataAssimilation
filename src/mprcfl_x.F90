SUBROUTINE RCFL_X(IM,JM,KM,FLD)
USE MYFRTPROF, ONLY : MYFRTPROF_WALL
USE RECFILTER

!-- RECURSIVE FILTER IN X-DIR

 USE SET_KND
 USE MYNETCDF

 IMPLICIT NONE

 INTEGER(I4)    :: IM, JM, KM

 REAL(R8)       :: FLD(IM,JM,2*KM)
 REAL(R8)       :: A(JM,IMAX), B(JM,IMAX), C(JM,IMAX)
 REAL(R8)       :: ALP(JM,IMAX), BTA(JM,IMAX)

 INTEGER(I4)    :: I,J,K, KTR,K2
 INTEGER(I4)    :: NIT

 CALL MYFRTPROF_WALL('RCFL_X: RECURSIVE FILTER X-DIR',0)

 NIT=RF_NTR

#ifdef SHARED_MEMORY
!$OMP PARALLEL DEFAULT(SHARED),PRIVATE(K2,K,A,B,C,I,J,KTR,ALP,BTA)
!$OMP DO SCHEDULE(DYNAMIC,1)
#endif
OUTER :   DO K2=1,2*KM
        K=K2
        IF(K>KM.AND..NOT.LLTSAPART.AND..NOT.LLRFSR) K=K-KM

        A(MJS:MJE,:) = 0._R8
        B(MJS:MJE,:) = 0._R8
        C(MJS:MJE,:) = 0._R8

        DO J=MJS,MJE
        DO I=1,IM
           A(J,INX(I,J,K)) = FLD(I,J,K2)
        ENDDO
        ENDDO

       ALP(MJS:MJE,:) = RF_AEX(MJS:MJE,:,K)
       BTA(MJS:MJE,:) = RF_BEX(MJS:MJE,:,K)

       DO KTR = 1,NIT

! POSITIVE DIRECTION
         IF( KTR.EQ.1 )THEN
           B(MJS:MJE,1) = (1.-ALP(MJS:MJE,1)) * A(MJS:MJE,1)
         ELSEIF( KTR.EQ.2 )THEN
            B(MJS:MJE,1) = A(MJS:MJE,1) / (1.+ALP(MJS:MJE,1))
         ELSE
            B(MJS:MJE,1) = (1.-ALP(MJS:MJE,1)) * &
            & (A(MJS:MJE,1)-ALP(MJS:MJE,1)**3 * A(:,2)) / &
            & (1.-ALP(MJS:MJE,1)**2)**2
         ENDIF

        DO J=2,IMX(K)
             B(MJS:MJE,J) = ALP(MJS:MJE,J)*B(MJS:MJE,J-1) + &
             & (1.-ALP(MJS:MJE,J))*A(MJS:MJE,J)
        ENDDO

! NEGATIVE DIRECTION
         IF( KTR.EQ.1 )THEN
           C(MJS:MJE,IMX(K)) = B(MJS:MJE,IMX(K)) / (1.+BTA(MJS:MJE,IMX(K)))
         ELSE
           C(MJS:MJE,IMX(K)) = (1.-BTA(MJS:MJE,IMX(K))) * &
           & (B(MJS:MJE,IMX(K))-BTA(MJS:MJE,IMX(K))**3 * &
           & B(MJS:MJE,IMX(K)-1)) / (1.-BTA(MJS:MJE,IMX(K))**2)**2
         ENDIF

        IF(LLDEBRCFLX) THEN
           WRITE(IOUNLOG,*) '-- ',K,KTR,SUM(BTA(:,1:IMX(K)))/REAL(JM*IMX(K)),&
                            & SUM(C(:,1:IMX(K)))/REAL(JM*IMX(K)),&
                            & SUM(B(:,1:IMX(K)))/REAL(JM*IMX(K))
           CALL FLUSH(IOUNLOG)
        ENDIF

         DO J=IMX(K)-1,1,-1
          C(MJS:MJE,J) = BTA(MJS:MJE,J)*C(MJS:MJE,J+1) + (1.-BTA(MJS:MJE,J))*B(MJS:MJE,J)
         ENDDO

         A(MJS:MJE,:) = C(MJS:MJE,:)

       ENDDO

        DO J=MJS:MJE
        DO I=1,IM
         FLD(I,J,K2) = A(J,INX(I,J,K))
        ENDDO
        ENDDO

   ENDDO OUTER
#ifdef SHARED_MEMORY
!$OMP END DO
!$OMP END PARALLEL
#endif
CALL FLUSH(IOUNLOG)
IF(LLDEBRCFLX) LLDEBRCFLX=.FALSE.

IF(LL_RFMPI) CALL RF_TO_ALL( 'RCFL_X',FLD)

 CALL MYFRTPROF_WALL('RCFL_X: RECURSIVE FILTER X-DIR',1)
END SUBROUTINE RCFL_X
