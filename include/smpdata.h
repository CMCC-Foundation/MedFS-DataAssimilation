!... Include file for shared memory parallelization configuration
!    Definitions also apply to recursive filter routines
!    if T and S are treated indipendently (LLTSAPART=.TRUE.) 
!    or in case of SR

! Index of variables in global arrays

INTEGER(KIND=i4),PARAMETER :: psv3d=4, &
                            & psv2d=1

INTEGER(KIND=i4),PARAMETER :: psvt=1, &
                            & psvs=2, &
                            & psvu=3, &
                            & psvv=4, &
                            & psve=5
