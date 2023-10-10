SUBROUTINE RCFL_XZ(IM,JM,FLD)
USE MYFRTPROF, ONLY : MYFRTPROF_WALL
USE RECFILTER

!-- RECURSIVE FILTER IN X-DIR

 USE SET_KND
 USE MYNETCDF

 IMPLICIT NONE

 INTEGER(I4)    :: IM, JM

 REAL(R8)       :: FLD(IM,JM)
 REAL(R8)       :: A(JM,IMAXZ), B(JM,IMAXZ), C(JM,IMAXZ)
 REAL(R8)       :: ALP(JM,IMAXZ), BTA(JM,IMAXZ)

 INTEGER(I4)    :: I,J,KTR
 INTEGER(I4)    :: NIT

 CALL MYFRTPROF_WALL('RCFL_XZ: RECURSIVE FILTER X-DIR',0)

 NIT=RF_NTR

        A(:,:) = 0.0
        B(:,:) = 0.0
        C(:,:) = 0.0

        DO J=1,JM
        DO I=1,IM
           A(J,INXZ(I,J)) = FLD(I,J)
        ENDDO
        ENDDO

       ALP(:,:) = RF_AEXZ(:,:)
       BTA(:,:) = RF_BEXZ(:,:)

       DO KTR = 1,NIT

! POSITIVE DIRECTION
         IF( KTR.EQ.1 )THEN
           B(:,1) = (1.-ALP(:,1)) * A(:,1)
         ELSEIF( KTR.EQ.2 )THEN
            B(:,1) = A(:,1) / (1.+ALP(:,1))
         ELSE
            B(:,1) = (1.-ALP(:,1)) * (A(:,1)-ALP(:,1)**3 * A(:,2)) / &
            & (1.-ALP(:,1)**2)**2
         ENDIF

        DO J=2,IMXZ
             B(:,J) = ALP(:,J)*B(:,J-1) + (1.-ALP(:,J))*A(:,J)
        ENDDO

! NEGATIVE DIRECTION
         IF( KTR.EQ.1 )THEN
           C(:,IMXZ) = B(:,IMXZ) / (1.+BTA(:,IMXZ))
         ELSE
           C(:,IMXZ) = (1.-BTA(:,IMXZ)) * (B(:,IMXZ)-BTA(:,IMXZ)**3 * &
           & B(:,IMXZ-1)) / (1.-BTA(:,IMXZ)**2)**2
         ENDIF

         DO J=IMXZ-1,1,-1
          C(:,J) = BTA(:,J)*C(:,J+1) + (1.-BTA(:,J))*B(:,J)
         ENDDO

         A(:,:) = C(:,:)

       ENDDO

        DO J=1,JM
        DO I=1,IM
         FLD(I,J) = A(J,INXZ(I,J))
        ENDDO
        ENDDO

 CALL MYFRTPROF_WALL('RCFL_XZ: RECURSIVE FILTER X-DIR',1)
END SUBROUTINE RCFL_XZ
