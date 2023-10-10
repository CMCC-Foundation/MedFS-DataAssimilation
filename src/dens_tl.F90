SUBROUTINE DENS_TL

!-----------------------------------------------------------------------
! CALCULATE DENSITY INCREMENTS
!-----------------------------------------------------------------------


  USE SET_KND
  USE GRD_STR
  USE BAL
  USE OCEANTOOLS

  IMPLICIT NONE

  INTEGER(I4)    :: K,I,J

   IF( BAL_NNEOS .EQ. 1 ) THEN
     DO K=1,GRD%KM
       GRD%DNS(:,:,K) = (-0.24*GRD%TEM(:,:,K) + 0.74*GRD%SAL(:,:,K))* &
                      & GRD%MSK(:,:,K)
     ENDDO
   ELSEIF( BAL_NNEOS .EQ. 2 ) THEN
     DO K=1,GRD%KM
        DO J=1,GRD%JM
           DO I=1,GRD%IM
              GRD%DNS(I,J,K) = RHO_UNESCOTL(BAL_SALB(I,J,K),BAL_TEMB(I,J,K),&
              & GRD%SAL(I,J,K),GRD%TEM(I,J,K),0._R8,.FALSE.)*GRD%MSK(I,J,K)
           ENDDO
        ENDDO
     ENDDO
   ELSE
       CALL ABOR1('DENS_TL: UNSUPPORTED E.O.S. OPTION')
   ENDIF

END SUBROUTINE DENS_TL
