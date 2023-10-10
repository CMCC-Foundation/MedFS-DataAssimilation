SUBROUTINE EIGENV(NCV,NCV2,B,EVC,EVA)

  USE SET_KND
  USE EOF_STR

  IMPLICIT NONE

  INTEGER(I4),INTENT(IN) :: NCV,NCV2
  REAL(R8), INTENT(IN)  :: B(ROS%KMT,ROS%KMT)
  REAL(R8), INTENT(OUT) :: EVC(ROS%KMT,ROS%NEOF)
  REAL(R8), INTENT(OUT) :: EVA(ROS%NEOF)

  INTEGER(I4)         :: JJ

! ARPACK SECTION

  INTEGER(I4)         :: NPNT, LWORKL,  INFO, IDO, IERR, NCONV,NEV
  INTEGER(I4)         :: IPARAM(11), IPNTR(11),I_ONE
  REAL(R8)            :: TOL, SIGMA, R_ONE, R_ZER
  LOGICAL             :: RVEC
  REAL(R8)            :: D_V(ROS%KMT,NCV), WORKL(NCV2)
  REAL(R8)            :: WORKD(3*ROS%KMT), D_S(NCV,2), RESID(ROS%KMT)
  LOGICAL             :: SELECT(NCV)
  CHARACTER           :: BMAT*1, WHICH*2

! PARTIAL EIGENVALUE DECOMPOSITION

!ARPACK

      BMAT   = 'I'
      WHICH  = 'LM'
      NPNT   = ROS%KMT
      NEV    =  ROS%NEOF
      LWORKL = NCV*(NCV+8)
      TOL    = 1.E-6_R8
      INFO   = 0
      IDO    = 0

      IPARAM(1) = 1
      IPARAM(3) = INT(ROS%NEOF*2.5)
      IPARAM(7) = 1

      R_ONE = 1._R8
      R_ZER = 0._R8
      I_ONE = 1

 10   CONTINUE

         CALL DSAUPD ( IDO, BMAT, NPNT, WHICH, NEV, TOL, RESID,                &
                       NCV, D_V, NPNT, IPARAM, IPNTR, WORKD, WORKL,            &
                       LWORKL, INFO )
         IF (IDO .EQ. -1 .OR. IDO .EQ. 1) THEN
           CALL DGEMV ( 'N', ROS%KMT, ROS%KMT, R_ONE, B, ROS%KMT,             &
                        WORKD(IPNTR(1)), I_ONE, R_ZER, WORKD(IPNTR(2)), I_ONE)
            GO TO 10
         END IF

      IF ( INFO .LT. 0 ) THEN
         PRINT *, ' '
         PRINT *, ' ERROR WITH _SAUPD, INFO = ', INFO
         PRINT *, ' CHECK DOCUMENTATION IN _SAUPD '
         PRINT *, ' '
      ELSE
         RVEC = .TRUE.
         CALL DSEUPD ( RVEC, 'ALL', SELECT, D_S, D_V, NPNT, SIGMA,             &
              BMAT, NPNT, WHICH, NEV, TOL, RESID, NCV, D_V, NPNT,              &
              IPARAM, IPNTR, WORKD, WORKL, LWORKL, IERR )
         IF ( IERR .NE. 0) THEN
             PRINT *, ' '
             PRINT *, ' ERROR WITH _SEUPD, INFO = ', IERR
             PRINT *, ' CHECK THE DOCUMENTATION OF _SEUPD. '
             PRINT *, ' '
         ELSE
             NCONV =  IPARAM(5)
             
             DO JJ=1,NCONV
                EVC(:,JJ) = D_V(:,NCONV-JJ+1)
                EVA(JJ)   = MAX(1.E-32_R8,D_S(NCONV-JJ+1,1))
                EVA(JJ)   = SQRT( EVA(JJ) )
             ENDDO
             IF(NCONV.LT.ROS%NEOF)THEN
                EVC(:,NCONV+1:ROS%NEOF) = 0.0
                EVA(NCONV+1:ROS%NEOF)   = 0.0
             ENDIF

         ENDIF
      ENDIF

! ----------------------------------------------------------------

END SUBROUTINE EIGENV
