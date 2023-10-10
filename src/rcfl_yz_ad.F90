SUBROUTINE RCFL_YZ_AD(IM,JM,FLD)

 USE RECFILTER
 USE SET_KND

 IMPLICIT NONE

 INTEGER(I4)    :: IM, JM

 REAL(R8)       :: FLD(IM,JM)
 REAL(R8)       :: A(IM,JMAXZ), B(IM,JMAXZ), C(IM,JMAXZ)
 REAL(R8)       :: ALP(IM,JMAXZ), BTA(IM,JMAXZ)

 INTEGER(I4)    :: I,J,KTR
 INTEGER(I4)    :: NIT

 NIT=RF_NTR

        A(:,:) = 0.0
        B(:,:) = 0.0
        C(:,:) = 0.0

        DO J=1,JM
         DO I=1,IM
            C(I,JNXZ(I,J)) = FLD(I,J)
         ENDDO
        ENDDO

         ALP(:,:) = RF_AEYZ(:,:)
         BTA(:,:) = RF_BEYZ(:,:)

         DO J=1,JMXZ-1
          C(:,J+1) = C(:,J+1) + BTA(:,J)*C(:,J)
          B(:,J)   = (1.-BTA(:,J))*C(:,J)
         ENDDO
         B(:,JMXZ) = B(:,JMXZ) + C(:,JMXZ) / (1.+BTA(:,JMXZ))

         DO J=JMXZ,2,-1
          B(:,J-1) = B(:,J-1) + ALP(:,J)*B(:,J)
          A(:,J) = A(:,J) + (1.-ALP(:,J))*B(:,J)
         ENDDO
         A(:,1) = A(:,1) + (1.-ALP(:,1)) * B(:,1)
         C(:,:) = A(:,:)

         B(:,:) = 0.0
         DO J=1,JMXZ-1
          C(:,J+1) = C(:,J+1) + BTA(:,J)*C(:,J)
          B(:,J)   = (1.-BTA(:,J))*C(:,J)
         ENDDO
         B(:,JMXZ  ) = B(:,JMXZ  ) + (1.-BTA(:,JMXZ)) * C(:,JMXZ) /&
         &  (1.-BTA(:,JMXZ)**2)**2
         B(:,JMXZ-1) = B(:,JMXZ-1) - (1.-BTA(:,JMXZ)) * BTA(:,JMXZ)**3&
         &  * C(:,JMXZ) / (1.-BTA(:,JMXZ)**2)**2

         A(:,:) = 0.0
         DO J=JMXZ,2,-1
          B(:,J-1) = B(:,J-1) + ALP(:,J)*B(:,J)
          A(:,J) = A(:,J) + (1.-ALP(:,J))*B(:,J)
         ENDDO
         A(:,1) = A(:,1) + B(:,1) / (1.+ALP(:,1))
         C(:,:) = A(:,:)


       DO KTR = 3,NIT
        B(:,:) = 0.0

         DO J=1,JMXZ-1
          C(:,J+1) = C(:,J+1) + BTA(:,J)*C(:,J)
          B(:,J)   = (1.-BTA(:,J))*C(:,J)
         ENDDO


           B(:,JMXZ  ) = B(:,JMXZ  ) + (1.-BTA(:,JMXZ)) * C(:,JMXZ) /&
           &  (1.-BTA(:,JMXZ)**2)**2
           B(:,JMXZ-1) = B(:,JMXZ-1) - (1.-BTA(:,JMXZ)) * BTA(:,JMXZ)**3&
           &  * C(:,JMXZ) / (1.-BTA(:,JMXZ)**2)**2

! POSITIVE DIRECTION (INVERSE OF)
        A(:,:) = 0.0

         DO J=JMXZ,2,-1
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
          FLD(I,J) = C(I,JNXZ(I,J))
         ENDDO
        ENDDO

END SUBROUTINE RCFL_YZ_AD
