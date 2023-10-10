SUBROUTINE GET_BYG

!-----------------------------------------------------------------------
!                                                                      !
! CALCULATE VERTICAL INTEGRAL OF BOUYANCY GRADIENT                     !
!                                                                      !
! VERSION 1: S.DOBRICIC 2007                                           !
!-----------------------------------------------------------------------


  USE SET_KND
  USE DRV_STR
  USE GRD_STR
  USE BMD_STR
  USE OCEANTOOLS
  USE MYFRTPROF, ONLY :MYFRTPROF_WALL
  USE MOD_BMD
  IMPLICIT NONE

  INTEGER(I4)    :: K,IX,I,J

CALL MYFRTPROF_WALL('GET_BYG: BUOYANCY GRADIENT',0)

IF( BMD_NNEOS .EQ. 1 ) THEN
     DO K=1,GRD%KM
       GRD%DNS(:,:,K) = (-0.24*GRD%TEM(:,:,K) + 0.74*GRD%SAL(:,:,K))* &
                      & 9.81/1020._R8* GRD%MSK(:,:,K)
     ENDDO
   ELSEIF( BMD_NNEOS .EQ. 2 ) THEN
     DO K=1,GRD%KM
        DO J=1,GRD%JM
           DO I=1,GRD%IM
              GRD%DNS(I,J,K) = RHO_UNESCOTL(BMD_SALB(I,J,K),BMD_TEMB(I,J,K),&
              & GRD%SAL(I,J,K),GRD%TEM(I,J,K),0._R8,.FALSE.)*9.81/1020._R8* GRD%MSK(I,J,K)
           ENDDO
        ENDDO
     ENDDO
   ELSE
       CALL ABOR1('GET_BYG: UNSUPPORTED E.O.S. OPTION')
   ENDIF


      GRD%B_X = 0.0
      GRD%B_Y = 0.0
      GRD%B_X(2:GRD%IM,:,1)=(GRD%DNS(2:GRD%IM,:,1)-GRD%DNS(1:GRD%IM-1,:,1))*&
       & GRD%DZ(1)*GRD%MSK(2:GRD%IM,:,1)*GRD%MSK(1:GRD%IM-1,:,1)
      GRD%B_Y(:,2:GRD%JM,1)=(GRD%DNS(:,2:GRD%JM,1)-GRD%DNS(:,1:GRD%JM-1,1))*&
       & GRD%DZ(1)*GRD%MSK(:,2:GRD%JM,1)*GRD%MSK(:,1:GRD%JM-1,1)

     DO K=2,GRD%KM

      GRD%B_X(2:GRD%IM,:,K) = GRD%B_X(2:GRD%IM,:,K-1) + &
      & (GRD%DNS(2:GRD%IM,:,K)-GRD%DNS(1:GRD%IM-1,:,K))*&
      & GRD%DZ(K)*GRD%MSK(2:GRD%IM,:,K)*GRD%MSK(1:GRD%IM-1,:,K)

      GRD%B_Y(:,2:GRD%JM,K) = GRD%B_Y(:,2:GRD%JM,K-1) + &
      & (GRD%DNS(:,2:GRD%JM,K)-GRD%DNS(:,1:GRD%JM-1,K))*&
      & GRD%DZ(K)*GRD%MSK(:,2:GRD%JM,K)*GRD%MSK(:,1:GRD%JM-1,K)

     ENDDO

      GRD%BX = 0.0
      GRD%BY = 0.0

     DO K=1,GRD%KM
      GRD%BX(2:GRD%IM,:) = GRD%BX(2:GRD%IM,:) + GRD%B_X(2:GRD%IM,:,K)*GRD%DZ(K) * &
      & GRD%MSK(2:GRD%IM,:,K)*GRD%MSK(1:GRD%IM-1,:,K)
      GRD%BY(:,2:GRD%JM) = GRD%BY(:,2:GRD%JM) + GRD%B_Y(:,2:GRD%JM,K)*GRD%DZ(K) * &
      & GRD%MSK(:,2:GRD%JM,K)*GRD%MSK(:,1:GRD%JM-1,K)
     ENDDO


#ifdef SHARED_MEMORY
!$OMP PARALLEL DEFAULT(SHARED), PRIVATE(IX)
!$OMP DO SCHEDULE(DYNAMIC)
#endif
    DO IX=1,GRD%IM-1
        GRD%BX(1+IX,:) = GRD%BX(1+IX,:)/GRD%DX(1+IX,:)
    ENDDO
#ifdef SHARED_MEMORY
!$OMP END DO
!$OMP DO SCHEDULE(DYNAMIC)
#endif
    DO IX=1,GRD%JM-1
        GRD%BY(:,1+IX) = GRD%BY(:,1+IX)/GRD%DY(:,1+IX)
    ENDDO
#ifdef SHARED_MEMORY
!$OMP END DO
!$OMP END PARALLEL
#endif

CALL MYFRTPROF_WALL('GET_BYG: BUOYANCY GRADIENT',1)
END SUBROUTINE GET_BYG
