C/*
C========== SCCS configuration management ==========================
C
C	SCCS	File name	: fcc_CDFP_PublicReadGrid.h
C	SCCS	Version		: 1.2
C	SCCS	Storage date	: 04/01/23
C
C==================================================================
C
C=========== RCS configuration management =========================
C
C	RCS	File name	: $RCSfile$
C	RCS	Version		: $Revision$
C	RCS	Storage date	: $Date$
C
C================================================================== 
C
C*****************************************************************************
C* MODULE : CDFP_PublicReadGrid
C*
C* ROLE :
C*    Fortran interface to some data types and constants.
C*    See CDFP_PublicReadGrid.h for more information.
C*    Function definitions of fortran interface are in CDFP_PublicReadGrid.c
C*
C*
C*****************************************************************************
C
C Modifications:
C===============
C
C 2001/10/17: Ph. Poilbarbe
C		Creation
C  
C*/


C/*
C==============================================================================
C===            E X P O R T E D    T Y P E S    A N D    D A T A            ===
C==============================================================================
C*/

C/* Semantic of data in grid */
C/* ATTENTION: Keep numbers since they are used for Fortran Interface */
      INTEGER*4	E_CDFP_GridDots
      INTEGER*4 E_CDFP_GridBoxes
      INTEGER*4	E_CDFP_GridDotsMercator
      INTEGER*4 E_CDFP_GridBoxesMercator
      PARAMETER (E_CDFP_GridDots		= 0)
      PARAMETER (E_CDFP_GridBoxes		= 1)
      PARAMETER (E_CDFP_GridDotsMercator	= 2)
      PARAMETER (E_CDFP_GridBoxesMercator	= 3)


