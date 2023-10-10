SUBROUTINE matinv2(a,n,IERR)

  USE SET_KND

  IMPLICIT NONE

  INTEGER(I4), INTENT(IN) :: n                   ! dimension of matrix
  REAL(R8), INTENT(INOUT) :: a(n,n)              ! matrix to be inverted
  INTEGER(I4), INTENT(OUT) :: IERR

  !work variables
  INTEGER(I4) :: i, icol, irow, j, k             ! loop variables
  INTEGER(I4) :: indxc(n), indxr(n), ipiv(n)     ! tracking arrays

  REAL(R8) :: big, pivinv, dum                   ! pivotting vals and tolerance
  REAL(R8), PARAMETER :: tol = 1.0e-6_R8
  REAL(R8) :: dumr(1,n), dumc(n,1)               ! temporary row/column copies

  ipiv(:) = 0
  IERR = 0

  DO k=1, n                 !loop over columns
    big = 0.0_R8            !big is largest element in ith row
    DO i=1, n
      IF (ipiv(i) /= 1) THEN
        DO j = 1, n           !search over columns on ith row
          IF (ipiv(j) == 0) THEN
            IF (ABS(a(i,j)) >= big) THEN
              big = ABS(a(i,j))
              irow = i
              icol = j       !so pivotal element identified
            END IF
          ELSE IF (ipiv(j) > 1) THEN
            PRINT*, j, k, 'matinv A: Singular matrix'
            IERR = 1
            RETURN
          END IF
        END DO
      END IF
    END DO
    ipiv(icol) = ipiv(icol) + 1
    !pivot element identified at A(I,J) so now interchange rows
    !keep track of changes using indxr and indxc
    IF (irow /= icol) THEN
      DO j=1, n
        dumr(1,j) = a(irow,j)
        a(irow,j) = a(icol,j)
        a(icol,j) = dumr(1,j)
      END DO
    END IF

    indxr(k) = irow      !where pivot is now located
    indxc(k) = icol      !where pivot was located

    IF (ABS(a(icol,icol)) <= tol ) THEN
      PRINT*, k, 'matinv B: Close to singular matrix'
      IERR = 2
      RETURN
    END IF

    !divide pivotal column by pivot element
    pivinv = 1.0_R8 / a(icol,icol)
    a(icol,icol) = 1.0_R8
    DO j=1, n
      a(icol,j) = a(icol,j) * pivinv
    END DO

    !now elimate elements either side of A(icol,icol) by GE
    DO i=1, n
      IF (i /= icol) THEN
        dum = a(i,icol)                      !big is dummy variable
        a(i,icol) = 0.                       !just to make sure
        DO j=1, n
          a(i,j) = a(i,j) - a(icol,j) * dum
        END DO
      END IF
    END DO
  END DO

  !Now need to rearrange columns after pivoting
  DO k=n,1,-1
    IF (indxr(k) /= indxc(k)) THEN
      DO i=1, n
        dumc(i,1) = a(i,indxr(k))
        a(i,indxr(k)) = a(i,indxc(k))
        a(i,indxc(k)) = dumc(i,1)
      END DO
    END IF
  END DO

END SUBROUTINE matinv2
