PROGRAM MASTER

!-----------------------------------------------------------------------
!                                                                      !
! THE MAIN DRIVER OF THE VARIATIONAL ANALYSIS                          !
!                                                                      !
! VERSION 1: S.DOBRICIC 2006                                           !
!   HORIZONTAL COVARIANCE WITH RECURSIVE FILTERS, VERTICAL WITH EOFS,  !
!   ASSIMILATION OF SATELLITE OBSERVATIONS OF SLA, IN SITU OBSERVATIONS!
!   BY XBT AND ARGO FLOATS                                             !
!                                                                      !
! VERSION 2: S.DOBRICIC 2006                                           !
!   MULTIGRID METHOD.                                                  !
!                                                                      !
! VERSION 2: S.DOBRICIC 2006                                           !
!   INTERNAL BOUNDARIES FOR HORIZONTAL                                 !
!                                                                      !
!-----------------------------------------------------------------------

 USE SET_KND
 USE MYFRTPROF, ONLY : MYFRTPROF_PRINT,MYFRTPROF_WALL
 USE IOUNITS, ONLY : IOUNLOG

 IMPLICIT NONE

 CALL MYFRTPROF_WALL('MASTER: PROGRAM',0)

!... CALL THE SETUP CHAIN
 CALL SETUP0

!... CALL THE VARIATIONAL JOB DRIVER ROUTINE
 CALL VARJOB

!... PREPARING TO EXIT
 CALL TERMIN0

 CALL MYFRTPROF_WALL('MASTER: PROGRAM',1)
 CALL MYFRTPROF_PRINT(-1,IOUNLOG)

!... EXIT
 CALL TERMIN1

END PROGRAM MASTER
