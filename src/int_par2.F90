SUBROUTINE INT_PAR2

!-----------------------------------------------------------------------
!                                                                      !
! CALCULATE INTERPOLATION PARAMETERS AND FIRST GUESS DEPARTURES        !
!                                                                      !
! VERSION 1: S.DOBRICIC 2006                                           !
! VERSION 2: A.S 2009 ADAPTED FOR AT-ONCE READING OF OBS
!-----------------------------------------------------------------------

 USE READFG
 USE IOUNITS
 USE GRD_STR
 USE RUN, ONLY : NCONF

 IMPLICIT NONE

 INTEGER(I4) :: JL, IAUX

!=== SLA PART

  !... PREPARE INTERPOLATION
  CALL INT_PAR_SLA

!=== INS PART

  !... PREPARE INTERPOLATION
  CALL INT_PAR_INS

!=== SST PART

  !... PREPARE INTERPOLATION
  CALL INT_PAR_SST

!=== SSS PART

  !... PREPARE INTERPOLATION
  CALL INT_PAR_SSS

END SUBROUTINE INT_PAR2
