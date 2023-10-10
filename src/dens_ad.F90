SUBROUTINE DENS_AD

!-----------------------------------------------------------------------
!                                                                      !
! CALCULATE DENSITY - ADJOINT
!                                                                      !
!-----------------------------------------------------------------------


  USE SET_KND
  USE GRD_STR
  USE BAL
  USE OCEANTOOLS

  IMPLICIT NONE

  INTEGER(I4)    :: K,I,J
  REAL(R8)       :: SAD, TAD

   IF( BAL_NNEOS .EQ. 1 ) THEN
     DO K=1,GRD%KM
       GRD%TEM_AD(:,:,K) = GRD%TEM_AD(:,:,K) - 0.24*GRD%DNS(:,:,K) * GRD%MSK(:,:,K)
       GRD%SAL_AD(:,:,K) = GRD%SAL_AD(:,:,K) + 0.74*GRD%DNS(:,:,K) * GRD%MSK(:,:,K)
     ENDDO
   ELSEIF( BAL_NNEOS .EQ. 2 ) THEN
     DO K=1,GRD%KM
        DO J=1,GRD%JM
           DO I=1,GRD%IM
              CALL RHO_UNESCOAD(GRD%DNS(I,J,K),BAL_SALB(I,J,K),BAL_TEMB(I,J,K),SAD, TAD)
              GRD%TEM_AD(I,J,K) = GRD%TEM_AD(I,J,K) + TAD*GRD%MSK(I,J,K)
              GRD%SAL_AD(I,J,K) = GRD%SAL_AD(I,J,K) + SAD*GRD%MSK(I,J,K)
           ENDDO
        ENDDO
     ENDDO
   ELSE
       CALL ABOR1('DENS_AD: UNSUPPORTED E.O.S. OPTION')
   ENDIF

END SUBROUTINE DENS_AD
