SUBROUTINE CNV_CTV_AD

!-----------------------------------------------------------------------
!                                                                      !
! CONVERT FROM CONTROL TO V - ADJOINT                                  !
!                                                                      !
! VERSION 1: S.DOBRICIC 2006                                           !
!-----------------------------------------------------------------------


 USE GRD_STR
 USE CTL_STR
 USE EOF_STR
 USE MYFRTPROF, ONLY : MYFRTPROF_WALL
 USE RUN
 USE WEAKLY
 USE HYBRID
 USE EOFS3D

 IMPLICIT NONE

 INTEGER(I4)     :: I,J,K, KK

 CALL MYFRTPROF_WALL('CNV_CTV_AD: CONVERT CONTROL VECTOR',0)

   KK = 0

   IF(NCONF .NE. 202) THEN
    DO K=1,ROS%NEOF
     DO J=1,GRD%JM
      DO I=1,GRD%IM
       KK = KK+1
       CTL%G_C(KK) = GRD%RO_AD(I,J,K)
      ENDDO
     ENDDO
    ENDDO
   ENDIF

   IF(LLVARMDT.OR.NCONF.EQ.202) THEN
    DO J=1,GRD%JM
     DO I=1,GRD%IM
       KK = KK+1
       CTL%G_C(KK) = GRD%CMDT_AD(I,J)
     ENDDO
    ENDDO
   ENDIF

   IF( LL_HYBRID ) THEN
    KK = CVH_S
    DO K=1,ROSH%NEOF
     DO J=1,GRD%JM
      DO I=1,GRD%IM
       CTL%G_C(KK) = GRD%ROH_AD(I,J,K)
       KK = KK+1
      ENDDO
     ENDDO
    ENDDO
   ENDIF

   IF( LL_EOFS3D ) CTL%G_C(CVE_S:CVE_E) = ROE_AD(:)

   IF( LL_HYBRID .AND. LL_ALPHA_CTL ) CTL%G_C(CVA_S:CVA_E) = ALPHA_AD(:)

   IF(LL_WEAKLY) CTL%G_C(CVWC4_S:CVWC4_E) = ROW_AD

 CALL MYFRTPROF_WALL('CNV_CTV_AD: CONVERT CONTROL VECTOR',1)
END SUBROUTINE CNV_CTV_AD
