SUBROUTINE SUCOV

!-----------------------------------------------------------------------
!                                                                      !
! DEFINE FILTER CONSTANTS, EOFS, ETC.                                  !
!                                                                      !
! VERSION 1: S.DOBRICIC 2006                                           !
!-----------------------------------------------------------------------

  USE SET_KND
  USE GRD_STR
  USE EOF_STR
  USE IOUNITS, ONLY : IOUNERR, IOUNLOG, IOUNOUT,IOUNEOF
  USE MYFRTPROF , ONLY : MYFRTPROF_WALL

  IMPLICIT NONE

  CALL MYFRTPROF_WALL('SUCOV: SET-UP B COVARIANCES',0)

  !.. RECURSIVE FILTER
  CALL SURF

  !... EMPIRICAL ORTHOGONAL FUNCTIONS
  CALL SUEOF

  CALL MYFRTPROF_WALL('SUCOV: SET-UP B COVARIANCES',1)
END
