              SUBROUTINE INIT_RANDOM_SEED()
                USE ISO_FORTRAN_ENV, ONLY: INT64
                IMPLICIT NONE
                INTEGER, ALLOCATABLE :: SEED(:)
                INTEGER :: I, N, UN, ISTAT, DT(8), PID
                INTEGER :: T

                CALL RANDOM_SEED(SIZE = N)
                ALLOCATE(SEED(N))
                ! FIRST TRY IF THE OS PROVIDES A RANDOM NUMBER GENERATOR
                   ! FALLBACK TO XOR:ING THE CURRENT TIME AND PID. THE
                   ! PID IS
                   ! USEFUL IN CASE ONE LAUNCHES MULTIPLE INSTANCES OF
                   ! THE SAME
                   ! PROGRAM IN PARALLEL.
                   CALL SYSTEM_CLOCK(T)
                   IF (T == 0) THEN
                      CALL DATE_AND_TIME(VALUES=DT)
                      T = (DT(1) - 1970) * 365 * 24 * 60 * 60 * 1000 &
                           + DT(2) * 31 * 24 * 60 * 60 * 1000 &
                           + DT(3) * 24 * 60 * 60 * 1000 &
                           + DT(5) * 60 * 60 * 1000 &
                           + DT(6) * 60 * 1000 + DT(7) * 1000 &
                           + DT(8)
                   END IF
                   DO I = 1, N
                      SEED(I) = LCG(T)
                   END DO
                CALL RANDOM_SEED(PUT=SEED)

CONTAINS

                FUNCTION LCG(S)
                  INTEGER :: LCG
                  INTEGER :: S
                  IF (S == 0) THEN
                     S = 104729
                  ELSE
                     S = MOD(S, 4294967296)
                  END IF
                  S = MOD(S * 279470273, 4294967291)
                  LCG = INT(MOD(S, INT(HUGE(0))), KIND(0))
                END FUNCTION LCG

              END SUBROUTINE INIT_RANDOM_SEED

