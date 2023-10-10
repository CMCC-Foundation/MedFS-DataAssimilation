!/*
!========== SCCS configuration management ==========================
!
!	SCCS	File name	: fcc_CDFP_PublicCommon.h
!	SCCS	Version		: 1.1
!	SCCS	Storage date	: 04/01/23
!
!==================================================================
!
!=========== RCS configuration management =========================
!
!	RCS	File name	: $RCSfile$
!	RCS	Version		: $Revision$
!	RCS	Storage date	: $Date$
!
!================================================================== 
!
!*****************************************************************************
!* MODULE : CDFP_PublicCommon
!*
!* ROLE :
!*    Fortran interface to some common data types and constants.
!*    See CDFP_PublicCommon.h for more information.
!*    Function definitions of fortran interface are in CDFP_PublicRead*.c
!*
!*
!*****************************************************************************
!
! Modifications:
!===============
!
! 2004/01/16: Ph. Poilbarbe, SU.DM0245
!		Creation
!  
!*/


!/*
!==============================================================================
!===            E X P O R T E D    T Y P E S    A N D    D A T A            ===
!==============================================================================
!*/

! Default/Missing value

REAL(8) ,PARAMETER :: dc_CDFP_DefReal8 = 18446744073709551616.0D+00
INTEGER(4),PARAMETER :: dc_CDFP_DefInteger32 = 2147483647

