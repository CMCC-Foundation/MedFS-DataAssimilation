/*
========== SCCS configuration management ==========================

	SCCS	File name	: CDFP_PublicCommon.h
	SCCS	Version		: 1.2
	SCCS	Storage date	: 04/04/05

==================================================================

=========== RCS configuration management =========================

	RCS	File name	: $RCSfile$
	RCS	Version		: $Revision$
	RCS	Storage date	: $Date$

================================================================== 

*****************************************************************************
* MODULE : CDFP_PublicCommon
*
* ROLE :
*    Common declarations for public routines used to read published data.
*    I M P O R TA N T	 N O T E: This file has an equivalent for
*				  FORTRAN calls called fcc_CDFP_PublicCommon.h
*				  Any modification of fortran exported
*				  data/type must be dupplicated
*
* 
* EXPORTED DATA/TYPES:
*     Tnv_CDFP_GridType 	      Semantic for the position and meaning of
*				      grid data
*
* IMPORTED DATA:
*     < none >
*
* INTERNAL DATA:
*     < none >
*
* PUBLIC FUNCTIONS:
*     < none >
*
*
* INTERNAL FUNCTIONS:
*     < none >
*
* EXCEPTIONS:
*    < none >
*
* MODULES OF SUPPORT LAYER USED:
*     < none > Parts of support layer needed is included herein in order
*	       to publish independant source files
*
*
*
*****************************************************************************

 Modifications:
===============

 2004/01/16: Ph. Poilbarbe. SU.DM0245
		Creation
 2004/04/05: Ph. Poilbarbe
		Adapted for MS-Windows


  
*/
#ifndef CDFP_PublicCommon_h
#define CDFP_PublicCommon_h


/* Identifies file version for SCCS configuration management */
#ident "@(#)CDFP_PublicCommon.h	1.2 04/04/05"

/* Identifies file version for CVS/RCS configuration management */
#ident "$Id$"

/*
** For constant M_PI
*/
#include <math.h>




/*
==============================================================================
===                P U B L I C    T Y P E S    A N D    D A T A            ===
==============================================================================
*/


/*
**
** General types used to make size clearer and easily portable
**
*/
typedef int			CDFP_Integer32;
typedef double			CDFP_Real8;
typedef char			*CDFP_String;

#define VERIFY_SIZES	assert(sizeof(CDFP_Integer32) == 4);	\
			assert(sizeof(CDFP_Real8) == 8)


/*
**
** Constants used too
**
*/
#define	 dc_CDFP_DefReal8       ((CDFP_Real8)18446744073709551616.0) /* 2^64 fully representable */
#define	 jc_CDFP_DefInteger32   ((CDFP_Integer32)0x7FFFFFFF)
#define  dc_CDFP_Pi		((CDFP_Real8)(M_PI))
#define  dc_CDFP_Epsilon         ((CDFP_Real8) 1.0E-50)








/*
==============================================================================
===                 F U N C T I O N    P R O T O T Y P E S                 ===
==============================================================================
*/

/*
** Under Win32, C functions called by fortran may be declared with __stdcall
** attribute
*/
#ifdef CDFP_WIN32_STDCALL
#define STDCALL __stdcall
#else
#define STDCALL
#endif


/*
** Some fortran compilers add one underscore to function name while others add
** two or none at all. This macro make the name of a C function that would be
** called from fortran
*/

#ifdef CDFP_NO_UNDERSCORE
#define FTN_NAME(X) STDCALL X
#else
#ifdef CDFP_TWO_UNDERSCORES
#define FTN_NAME(X)	STDCALL X##__
#else
#define FTN_NAME(X)	STDCALL X##_
#endif
#endif



#endif
