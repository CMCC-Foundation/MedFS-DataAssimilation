SUBROUTINE RCFL_XZ_AD(IM,JM,FLD)
USE RECFILTER


 USE SET_KND

 IMPLICIT NONE

 INTEGER(I4)    :: IM, JM

 REAL(R8)       :: FLD(IM,JM)
 REAL(R8)       :: A(JM,IMAXZ), B(JM,IMAXZ), C(JM,IMAXZ)
 REAL(R8)       :: ALP(JM,IMAXZ), BTA(JM,IMAXZ)

 INTEGER(I4)    :: I,J,KTR
 INTEGER(I4)    :: NIT

 NIT=RF_NTR

        A(:,:) = 0.0
        B(:,:) = 0.0
        C(:,:) = 0.0

        DO J=1,JM
        DO I=1,IM
           C(J,INXZ(I,J)) = FLD(I,J)
        ENDDO
        ENDDO

         ALP(:,:) = RF_AEXZ(:,:)
         BTA(:,:) = RF_BEXZ(:,:)

         DO J=1,IMXZ-1
          C(:,J+1) = C(:,J+1) + BTA(:,J)*C(:,J)
          B(:,J)   = (1.-BTA(:,J))*C(:,J)
         ENDDO
           B(:,IMXZ) = B(:,IMXZ) + C(:,IMXZ) / (1.+BTA(:,IMXZ))

         DO J=IMXZ,2,-1
          B(:,J-1) = B(:,J-1) + ALP(:,J)*B(:,J)
          A(:,J) = A(:,J) + (1.-ALP(:,J))*B(:,J)
         ENDDO
           A(:,1) = A(:,1) + (1.-ALP(:,1)) * B(:,1)
         C(:,:) = A(:,:)

        B(:,:) = 0.0
         DO J=1,IMXZ-1
          C(:,J+1) = C(:,J+1) + BTA(:,J)*C(:,J)
          B(:,J)   = (1.-BTA(:,J))*C(:,J)
         ENDDO
           B(:,IMXZ  ) = B(:,IMXZ  ) + (1.-BTA(:,IMXZ)) * C(:,IMXZ) /&
           &  (1.-BTA(:,IMXZ)**2)**2
           B(:,IMXZ-1) = B(:,IMXZ-1) - (1.-BTA(:,IMXZ)) * BTA(:,IMXZ)**3 &
           & * C(:,IMXZ) / (1.-BTA(:,IMXZ)**2)**2

        A(:,:) = 0.0
         DO J=IMXZ,2,-1
          B(:,J-1) = B(:,J-1) + ALP(:,J)*B(:,J)
          A(:,J) = A(:,J) + (1.-ALP(:,J))*B(:,J)
         ENDDO
           A(:,1) = A(:,1) + B(:,1) / (1.+ALP(:,1))
         C(:,:) = A(:,:)


       DO KTR = 3,NIT

        B(:,:) = 0.0
         DO J=1,IMXZ-1
          C(:,J+1) = C(:,J+1) + BTA(:,J)*C(:,J)
          B(:,J)   = (1.-BTA(:,J))*C(:,J)
         ENDDO

           B(:,IMXZ  ) = B(:,IMXZ  ) + (1.-BTA(:,IMXZ)) * C(:,IMXZ) /&
           &  (1.-BTA(:,IMXZ)**2)**2
           B(:,IMXZ-1) = B(:,IMXZ-1) - (1.-BTA(:,IMXZ)) * BTA(:,IMXZ)**3 &
           & * C(:,IMXZ) / (1.-BTA(:,IMXZ)**2)**2

        A(:,:) = 0.0
         DO J=IMXZ,2,-1
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
         FLD(I,J) = C(J,INXZ(I,J))
        ENDDO
        ENDDO

END SUBROUTINE RCFL_XZ_AD
