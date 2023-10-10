SUBROUTINE RCFL_YZ(IM,JM,FLD)
USE RECFILTER
USE MYFRTPROF, ONLY : MYFRTPROF_WALL

 USE SET_KND

 IMPLICIT NONE

 INTEGER(I4)    :: IM, JM

 REAL(R8)       :: FLD(IM,JM)
 REAL(R8)       :: A(IM,JMAXZ), B(IM,JMAXZ), C(IM,JMAXZ)
 REAL(R8)       :: ALP(IM,JMAXZ), BTA(IM,JMAXZ)

 INTEGER(I4)    :: I,J,KTR
 INTEGER(I4)    :: NIT

 CALL MYFRTPROF_WALL('RCFL_YZ: RECURSIVE FILTER Y-DIR',0)

        A(:,:) = 0.0
        B(:,:) = 0.0
        C(:,:) = 0.0

        DO J=1,JM
         DO I=1,IM
            A(I,JNXZ(I,J)) = FLD(I,J)
         ENDDO
        ENDDO

         ALP(:,:) = RF_AEYZ(:,:)
         BTA(:,:) = RF_BEYZ(:,:)

       DO KTR = 1,NIT


! POSITIVE DIRECTION
         IF( KTR.EQ.1 )THEN
           B(:,1) = (1.-ALP(:,1)) * A(:,1)
         ELSEIF( KTR.EQ.2 )THEN
            B(:,1) = A(:,1) / (1.+ALP(:,1))
         ELSE
            B(:,1) = (1.-ALP(:,1)) * (A(:,1)-ALP(:,1)**3 * A(:,2)) /&
            &  (1.-ALP(:,1)**2)**2
         ENDIF

        DO J=2,JMXZ
             B(:,J) = ALP(:,J)*B(:,J-1) + (1.-ALP(:,J))*A(:,J)
        ENDDO
! NEGATIVE DIRECTION
         IF( KTR.EQ.1 )THEN
           C(:,JMXZ) = B(:,JMXZ) / (1.+BTA(:,JMXZ))
         ELSE
           C(:,JMXZ) = (1.-BTA(:,JMXZ)) * (B(:,JMXZ)-BTA(:,JMXZ)**3 *&
           &  B(:,JMXZ-1)) / (1.-BTA(:,JMXZ)**2)**2
         ENDIF

         DO J=JMXZ-1,1,-1
          C(:,J) = BTA(:,J)*C(:,J+1) + (1.-BTA(:,J))*B(:,J)
         ENDDO

         A(:,:) = C(:,:)

       ENDDO

        DO J=1,JM
         DO I=1,IM
          FLD(I,J) = A(I,JNXZ(I,J))
         ENDDO
        ENDDO

 CALL MYFRTPROF_WALL('RCFL_YZ: RECURSIVE FILTER Y-DIR',1)

END SUBROUTINE RCFL_YZ
