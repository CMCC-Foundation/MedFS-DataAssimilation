SUBROUTINE UNBIAS_SLA_AT

 USE SET_KND
 USE GRD_STR
 USE OBS_STR
 USE IOUNITS
 USE RUN, ONLY : LL_UNBIAS_SLA_AT, ZZ_SLA_BIAS
 USE MPIREL

 IMPLICIT NONE

 INTEGER(I4)   ::  K, KK
 REAL(R8) :: ZBIAS
 INTEGER(I4)   ::  i_start, i, iter,flag
 REAL(R8)      ::  sum_res, nb_obs, dxx, dyy, dsm, n_min

#include "obs_events.h"

 K=0
 KK=0
 ZBIAS=0._R8
 
 IF(SLA%NO.LE.0) RETURN
 !DO iter = 1,3
 
 !bias
 SLA%BIA(:) = 0.0
 ! Maximum distance between two consecutive observations (in km)
 dsm = 100. 
 ! Minimum number of observations in a section
 n_min = 4.0
 i_start = 1
 DO K=2,SLA%NO

    ! Compute distance between two consecutive observations 
    dxx = 6371.*3.14/180. * (SLA%LON(K)-SLA%LON(K-1)) * &
           cos(SLA%LAT(K)*3.14/180.)
    dyy = 6371.*3.14/180. * (SLA%LAT(K)-SLA%LAT(K-1))

    ! If the distance is bigger than dsm, we compute the bias
    ! for the section before
    IF((sqrt(dxx**2+dyy**2).GT.dsm) .AND. K.GT.i_start) THEN
       sum_res = 0.0
       nb_obs = 0.0
       DO i=i_start,K-1
          IF(SLA%FLC(i).eq.1)then
             sum_res = sum_res + SLA%RES(i)
             nb_obs  = nb_obs + 1.0
          ENDIF
       ENDDO
       sum_res = sum_res/MAX(nb_obs,1.0_R8)
       IF(nb_obs .GT. n_min) THEN
          DO i=i_start,K-1
             SLA%RES(i) = SLA%RES(i) - sum_res
             SLA%BIA(i) = sum_res
          ENDDO
       ELSE
          DO i=i_start,K-1
             SLA%RES(i) = SLA%RES(i) - sum_res
             SLA%BIA(i) = sum_res
             SLA%FLC(i) = 0
             SLA%EVE(i) = keve_TFEW
          ENDDO 
       ENDIF
       ! reset the starting point for the next section
       i_start = K
    ! If we arrived at the end of the observations,
    ! we compute the last bias
    ELSE IF((K.EQ.SLA%NO) .AND. (K.GE.i_start)) THEN
       sum_res = 0.0
       nb_obs = 0.0
       DO i=i_start,K
          IF(SLA%FLC(i).eq.1)then
             sum_res = sum_res + SLA%RES(i)
             nb_obs  = nb_obs + 1.0
          ENDIF
       ENDDO
       sum_res = sum_res/MAX(nb_obs,1.0_R8)
       IF(nb_obs .GT. n_min) THEN
          DO i=i_start,K
             SLA%RES(i) = SLA%RES(i) - sum_res
             SLA%BIA(i) = sum_res
          ENDDO
       ELSE
          DO i=i_start,K
             SLA%RES(i) = SLA%RES(i) - sum_res
             SLA%BIA(i) = sum_res
             SLA%FLC(i) = 0
             SLA%EVE(i) = keve_TFEW
          ENDDO
       ENDIF
    ENDIF
 ENDDO
 
 ! ENDDO ! iter
 
 IF( LL_UNBIAS_SLA_AT ) THEN
    WRITE(IOUNLOG,*) ' THE ALONG TRACK BIAS IS REMOVED'
 ELSE
    WRITE(IOUNLOG,*) ' THE BIAS IS NOT REMOVED'
    ZZ_SLA_BIAS = 0._R8
 ENDIF

 CALL FLUSH(IOUNLOG)

END SUBROUTINE UNBIAS_SLA_AT
