SUBROUTINE INI_ITR

!-----------------------------------------------------------------------
!                                                                      !
! DEFINE THE GRID                                                      !
!                                                                      !
! VERSION 1: S.DOBRICIC 2006                                           !
!-----------------------------------------------------------------------


  USE SET_KND
  USE DRV_STR
  USE GRD_STR
  USE EOF_STR
  USE CTL_STR

  IMPLICIT NONE

  INTEGER(I4)               :: I, J, K , II, JJ, KK
  REAL(R8)                  :: RI, RJ, P, Q
  REAL(R8)                  :: DIV_X, DIV_Y
  REAL(R8),    ALLOCATABLE  :: PQ1(:,:), PQ2(:,:), PQ3(:,:), PQ4(:,:)
  INTEGER(I4), ALLOCATABLE  :: I1(:,:), J1(:,:)

  INTEGER(I4) :: IZS, INDIC, NITER
  REAL(R8) :: RZS
  REAL(R8) :: DZS
  REAL(R8), DIMENSION(CTL%N)    :: ZY, X, XGJ
  REAL(R8)       :: XJS, XJE,Z


   ALLOCATE( PQ1(GRD%IM,GRD%JM) )
   ALLOCATE( PQ2(GRD%IM,GRD%JM) )
   ALLOCATE( PQ3(GRD%IM,GRD%JM) )
   ALLOCATE( PQ4(GRD%IM,GRD%JM) )
   ALLOCATE(  I1(GRD%IM,GRD%JM) )
   ALLOCATE(  J1(GRD%IM,GRD%JM) )

! ---
! INTERPOLATE BETWEEN GRIDS

     DO JJ=1,GRD%JM
     DO II=1,GRD%IM
       RI=MAX(1.,MIN(REAL(DRV%IM-1),REAL(II-1)/REAL(DRV%RATIO(DRV%KTR)) + 1.))
        I=INT(RI)
       P=RI-I
       RJ=MAX(1.,MIN(REAL(DRV%JM-1),REAL(JJ-1)/REAL(DRV%RATIO(DRV%KTR)) + 1.))
        J=INT(RJ)
       Q=RJ-J

           I1(II,JJ) = I
           J1(II,JJ) = J

         DIV_Y =  (1.-Q) * MAX(DRV%MSK(I,J  ),DRV%MSK(I+1,J  ))      &
                 +    Q  * MAX(DRV%MSK(I,J+1),DRV%MSK(I+1,J+1))
         DIV_X =  (1.-P) * DRV%MSK(I  ,J) + P * DRV%MSK(I+1,J)
          PQ1(II,JJ) = DRV%MSK(I,J)                                  &
                      * MAX(DRV%MSK(I,J),DRV%MSK(I+1,J))             &
                       * (1.-P) * (1.-Q)                             &
                      /( DIV_X * DIV_Y + 1.E-16 )
          PQ2(II,JJ) = DRV%MSK(I+1,J)                                &
                      * MAX(DRV%MSK(I,J),DRV%MSK(I+1,J))             &
                      *     P  * (1.-Q)                              &
                      /( DIV_X * DIV_Y + 1.E-16 )
         DIV_X =  (1.-P) * DRV%MSK(I  ,J+1) + P * DRV%MSK(I+1,J+1)
          PQ3(II,JJ) = DRV%MSK(I,J+1)                                &
                      * MAX(DRV%MSK(I,J+1),DRV%MSK(I+1,J+1))         &
                      * (1.-P) *     Q                               &
                      /( DIV_X * DIV_Y + 1.E-16 )
          PQ4(II,JJ) = DRV%MSK(I+1,J+1)                              &
                      * MAX(DRV%MSK(I,J+1),DRV%MSK(I+1,J+1))         &
                      *     P  *     Q                               &
                      /( DIV_X * DIV_Y + 1.E-16 )

     ENDDO
     ENDDO


    DO K = 1,ROS%NEOF
     DO JJ=1,GRD%JM
     DO II=1,GRD%IM
        I=I1(II,JJ)
        J=J1(II,JJ)
       GRD%RO(II,JJ,K) = PQ1(II,JJ) * DRV%RO(I,J  ,K) + PQ2(II,JJ) * DRV%RO(I+1,J  ,K)  &
                       + PQ3(II,JJ) * DRV%RO(I,J+1,K) + PQ4(II,JJ) * DRV%RO(I+1,J+1,K)
     ENDDO
     ENDDO
    ENDDO

! ---
! RECONSTRUCT THE CONTROL VECTOR
       KK = 0
   DO K=1,ROS%NEOF
    DO J=1,GRD%JM
     DO I=1,GRD%IM
       KK = KK+1
       CTL%X_C(KK) = GRD%RO(I,J,K)/DRV%RATIO(DRV%KTR)
     ENDDO
    ENDDO
   ENDDO


   DEALLOCATE( DRV%RO, DRV%RO_AD, DRV%MSK)
   DEALLOCATE( PQ1 )
   DEALLOCATE( PQ2 )
   DEALLOCATE( PQ3 )
   DEALLOCATE( PQ4 )
   DEALLOCATE(  I1 )
   DEALLOCATE(  J1 )

! CALCULATE THE COST FUNCTION AND ITS GRADIENT

   CALL COSTF(INDIC,CTL%N,CTL%X_C,XJS,XGJ,NITER,IZS,RZS,DZS)


END SUBROUTINE INI_ITR
