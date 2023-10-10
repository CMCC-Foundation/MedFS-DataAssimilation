SUBROUTINE CNV_CTV

!-----------------------------------------------------------------------
!                                                                      !
! CONVERT FROM CONTROL TO V                                            !
!                                                                      !
! VERSION 1: S.DOBRICIC 2006                                           !
!-----------------------------------------------------------------------


 USE SET_KND
 USE GRD_STR
 USE CTL_STR
 USE EOF_STR
 USE MYFRTPROF, ONLY : MYFRTPROF_WALL
 USE RUN, ONLY : LLVARMDT, NCONF
 USE HYBRID
 USE EOFS3D

 IMPLICIT NONE

 INTEGER(I4)   :: I,J,K, KK

 CALL MYFRTPROF_WALL('CNV_CTV: CONTROL VECTOR TO V',0)

 KK = 0

  IF( NCONF .NE. 202 ) THEN
   DO K=1,ROS%NEOF
    DO J=1,GRD%JM
     DO I=1,GRD%IM
       KK = KK+1
       GRD%RO(I,J,K) = CTL%X_C(KK)
     ENDDO
    ENDDO
   ENDDO
  ENDIF

  IF(LLVARMDT .OR. NCONF .EQ. 202) THEN
   DO J=1,GRD%JM
    DO I=1,GRD%IM
      KK = KK+1
      GRD%CMDT(I,J) = CTL%X_C(KK)
    ENDDO
   ENDDO
  ENDIF

  IF( LL_HYBRID ) THEN
    KK = CVH_S
    DO K=1,ROSH%NEOF
     DO J=1,GRD%JM
      DO I=1,GRD%IM
       GRD%ROH(I,J,K) = CTL%X_C(KK)
       KK = KK+1
      ENDDO
     ENDDO
    ENDDO
  ENDIF

  IF( LL_EOFS3D ) ROE(:) = CTL%X_C(CVE_S:CVE_E)

 CALL MYFRTPROF_WALL('CNV_CTV: CONTROL VECTOR TO V',1)
END SUBROUTINE CNV_CTV
