/*
========== SCCS configuration management ==========================

	SCCS	File name	: CDFP_PublicReadGrid.h
	SCCS	Version		: 1.2
	SCCS	Storage date	: 04/01/23

==================================================================

=========== RCS configuration management =========================

	RCS	File name	: $RCSfile$
	RCS	Version		: $Revision$
	RCS	Storage date	: $Date$

==================================================================

*****************************************************************************
* MODULE : CDFP_PublicReadGrid  				         V1.1
*
* ROLE :
*    Interface to manipulate NetCdf files containing grids
*    I M P O R TA N T	 N O T E: This file has an equivalent for
*				  FORTRAN calls called fcc_CDFP_PublicRead.h
*				  Any modification of fortran exported
*				  data/type must be dupplicated
*
*    NOTE: When reading grids, data equal to the _FillValue attribute are
*	   always returned as dc_CDFP_DefReal8 whatever the _FillValue
*	   and the internal data types are. This allows consistant and
*	   unique way of testing missing data.
*
* PREREQUISITES:
*     When linking with a main program, the netcdf library V3.5 or
*     newer must be specified (option -lnetcdf), otherwise, functions
*     whose names begins with 'nc_' will be marqued as 'unresolved
*     externals'.
*
* EXPORTED DATA/TYPES:
*     Tnv_CDFP_GridType 	      Semantic for the position and meaning of
*				      grid data
*
* PUBLIC FUNCTIONS:
*     - Reads a grid from a netcdf grid file	      CDFP_ReadGrid
*     - Fortran interface (fcc=Fortran Calls C)       fcc_CDFP_ReadGrid
*
*****************************************************************************

 Modifications:
===============

 2001/10/10: F. Mertz
		Creation

 2001/10/17: Ph. Poilbarbe
		Modified to conform to standard names and format.
		Program converted to function. Fortran interface added

 2003/02/28: Ph. Poilbarbe
		Units atribute misspelled

 2004/01/16: Ph. Poilbarbe. SU.DM0245
		Common declarations for CDFP_PublicReadGrid and
		CDFP_PublicReadAlongTrackData grouped in common included file.
		Comments revamped.
		Fortran interface name adapted to fortran compilers which add
		two underscores instead of one (macro FTN_NAME)
		Verification of sizes of basic data (macro VERIFY_SIZES)


*/
#ifndef CDFP_PublicReadGrid_h
#define CDFP_PublicReadGrid_h


/* Identifies file version for SCCS configuration management */
#ident "@(#)CDFP_PublicReadGrid.h	1.2 04/01/23"

/* Identifies file version for CVS/RCS configuration management */
#ident "$Id$"


#include "CDFP_PublicCommon.h"




/*
==============================================================================
===                P U B L I C    T Y P E S    A N D    D A T A            ===
==============================================================================
*/

/*
**
** Identification of grid type.
**
** Change the semantic of sizes and the positions of data
**
** For a DOTS grid, each value is the sample value at the exact position
** of the dot.
**
** For a BOXES grid, each value sample is the value for the south-west
** corner or the area and has the same value for the whole surface until
** another box is encoutered. For example, a box cannot have a 90deg-north
** since there is no surface after this point.
**
** For MERCATOR grids, the latitude step is the one at the equator even
** if the equator is not in the grid.
**
**
*/
typedef enum {
  E_CDFP_GridDots		= 0,
  E_CDFP_GridBoxes		= 1,
  E_CDFP_GridDotsMercator	= 2,
  E_CDFP_GridBoxesMercator	= 3
} Tnv_CDFP_GridType;






/*
==============================================================================
===                 F U N C T I O N    P R O T O T Y P E S                 ===
==============================================================================
*/







/*
*****************************************************************************
*
* PUBLIC FUNCTION: CDFP_ReadGrid
*
*
* DESCRIPTION:
*    Read a grid from grid netcdf file
*
* RETURNED VALUE:
*    Error status: 0: No error Other: error
*
* PARAMETERS:
*   Name	      Unit    Description
*   ------------------------------------------------------------------------
*   Itv_FileName      /       Name of the NetCdf grid file
*   Ijv_GridNumber    /       Identifer of grid to reead
*   Onv_GridType      /       Grid type (Dots/Boxes/Mercator dots/Mercator boxes)
*   Otv_Unit	      /       Unit of the grid (text)
*   Ijv_UnitBufferSize
*		      /       Size of Otv_Unit buffer. If actual unit text
*			      is longer than this size, the end of the string
*			      is lost.
*   Odq_GridData      Otv_Unit
*			      Pointer to grid data (NbLon*NbLat*Depth data
*			      of type CDFP_Real8).
*			      If NULL no data is returned nor allocated
*   Odq_Longitudes    Degrees Pointer to longitude values (NbLon).
*			      If NULL no latitude is returned nor allocated
*   Odq_Latitudes     Degrees Pointer to latitude values (NbLat).
*			      If NULL no latitude is returned nor allocated
*   Ojv_NbLongitudes  /       Number of longitude values.
*   Ojv_NbLatitudes   /       Number of latitude values.
*   Ojv_GridDepth     /       Number layers of the grid (third dimension)
*
* NOTES:
*   If Odw_XXX is null, nothing is returned in the corresponding parameter.
*   if not, the pointer is filled with the address of allocated data. freeing
*   allocated data is the responsability of the caller.
*   In all cases NbLongitudes, NbLatitudes and GridDepth are returned
*   So using NULL values can be used ti tune what data you want or
*   get only dimensions.
*
*   Data are stored in Odw_GridData with latitude index varying first.
*   In a one dimonsionnal table, the pos function giving the index
*   of one element relatively to the Lat, Lon and depth indexes (all
*   starting at 0) is: Pos(Lat,Lon,Depth)=(Lon*NbLat+Lat)*Depth + Depth
*   It corresponds to a 3 dimensionnal declaration:
*     FORTRAN: Dimension GRID(Depth, NbLat, NbLon)
*     C      : Type GRID[NbLon][NbLat][Depth];
*     Where: NbLat, NbLon & Depth are the number of elements in latitude,
*	      longitude and depth respectively
*
* INTERNAL DATA ACCESSED (Side effect)
*    < none >
*
*****************************************************************************
*/
CDFP_Integer32 CDFP_ReadGrid
			(CDFP_String		Itv_FileName,
			 CDFP_Integer32		Ijv_GridNumber,
			 Tnv_CDFP_GridType	*Onv_GridType,
			 CDFP_String		Otv_Unit,
			 CDFP_Integer32		Ijv_UnitBufferSize,
			 CDFP_Real8		**Odq_GridData,
			 CDFP_Real8		**Odq_Longitudes,
			 CDFP_Real8		**Odq_Latitudes,
			 CDFP_Integer32		*Ojv_NbLongitudes,
			 CDFP_Integer32		*Ojv_NbLatitudes,
			 CDFP_Integer32		*Ojv_GridDepth)
;




/*
==============================================================================
==============================================================================
==============================================================================
==============================================================================


==============================================================================
===                   F O R T R A N    I N T E R F A C E                   ===
==============================================================================
*/









/*
*****************************************************************************
*
* PUBLIC FUNCTION: fcc_CDFP_ReadGrid
*
*
* DESCRIPTION:
*    Read a grid from grid netcdf file
*
* RETURNED VALUE:
*    < none >
*
* PARAMETERS:
*   Name	      Unit    Description
*   ------------------------------------------------------------------------
*   Itv_FileName      /       Name of the NetCdf grid file
*   Ijv_GridNumber    /       Identifer of grid to reead
*   Ojv_GridType      /       Grid type (Dots/Boxes/Mercator dots/Mercator boxes)
*   Otv_Unit	      /       Unit of the grid (text)
*   Ijv_LonDim        /        for the grid. See discussion below
*   Ijv_LatDim        /        Where xxxDim are maximum dimensions allowed
*   Ijv_DepthDim      /        REAL*8 GridData(DepthDim, LatDim, LonDim)
*   Odw_GridData      Otv_Unit This array is FORTRAN declared
*   Odw_Longitudes    Degrees  Array declared with LatDim dimension
*   Odw_Latitudes     Degrees  Array declared with LomDim dimension
*			      If NULL no latitude is returned nor allocated
*   Ojv_NbLongitudes  /       Number of longitude values.
*   Ojv_NbLatitudes   /       Number of latitude values.
*   Ojv_GridDepth     /       Number layers of the grid (third dimension)
*   Ojv_Status        /       Return satus: 0=Ok, Other=Not ok
*
* FORTRAN EQUIVALENT DECLARATION:
*  SUBROUTINE fcc_CDFP_ReadGrid
*  &		  (Itv_FileName,
*  &		   Ijv_GridNumber,
*  &		   Ojv_GridType,
*  &		   Otv_Unit,
*  &		   Ijv_LonDim,
*  &		   Ijv_LatDim,
*  &		   Ijv_DepthDim,
*  &		   Odw_GridData,
*  &		   Odw_Latitudes,
*  &		   Odw_Longitudes,
*  &		   Ojv_NbLongitudes,
*  &		   Ojv_NbLatitudes,
*  &		   Ojv_GridDepth
*  &		   Ojv_Status)
*    CHARACTER*(*) Itv_FileName
*    INTEGER*4     Ijv_GridNumber
*    INTEGER*4     Ojv_GridType     ! See fcc_CDFP_PublicReadGrid.h
*    CHARACTER*(*) Otv_Unit
*    INTEGER*4     Ijv_LonDim
*    INTEGER*4     Ijv_LatDim
*    INTEGER*4     Ijv_DepthDim
*    REAL*8	   Odw_GridData(Ijv_DepthDim, Ijv_LatDim, Ijv_LonDim)
*    REAL*8	   Odw_Longitudes(Ijv_LonDim)
*    REAL*8	   Odw_Latitudes(Ijv_LatDim)
*    INTEGER*4     Ojv_NbLongitudes
*    INTEGER*4     Ojv_NbLatitudes
*    INTEGER*4     Ojv_GridDepth
*    INTEGER*4     Ojv_Status
*
* NOTES:
*   If Ijv_LatDim is <=0, it is like setting GridData and Latitude to
*     NULL in C (nothing is returned in them)
*   If Ijv_LonDim is <=0, it is like setting GridData and Longitude to
*     NULL in C (nothing is returned in them)
*   If Ijv_DepthDim is <=0, it is like setting GridData to NULL in C
*     (nothing is returned in it)
*   if not, the arrays are filled with the data according to the declared
*   dimensions (Ijv_xxxDim) and the actual ones (Ojv_XXX).
*
*   If declared array dimensions are bigger than those of the grid the
*   'holes' are undefined (in fact they are unchanged)
*
*   Data are stored in Odw_GridData with latitude index varying first.
*   It corresponds to a 3 dimensionnal declaration:
*     FORTRAN: Dimension GRID(DepthDim, LatDim, LonDim)
*     Where: LonDim, latDim & DepthDim are the declared number of elements
*	     in longitude, latitude and depth respectively
*
* INTERNAL DATA ACCESSED (Side effect)
*    < none >
*
*****************************************************************************
*/
void FTN_NAME(fcc_cdfp_readgrid)
			(char			*Itv_FileName,
			 CDFP_Integer32		*Ijv_GridNumber,
			 CDFP_Integer32		*Ojv_GridType,
			 char			*Otv_Unit,
			 CDFP_Integer32		*Ijv_LonDim,
			 CDFP_Integer32		*Ijv_LatDim,
			 CDFP_Integer32		*Ijv_DepthDim,
			 CDFP_Real8		*Odw_GridData,
			 CDFP_Real8		*Odw_Longitudes,
			 CDFP_Real8		*Odw_Latitudes,
			 CDFP_Integer32		*Ojv_NbLongitudes,
			 CDFP_Integer32		*Ojv_NbLatitudes,
			 CDFP_Integer32		*Ojv_GridDepth,
			 CDFP_Integer32		*Ojv_Status,

/*******************Automatically added by fortran do not use them as actual parameters */
			 CDFP_Integer32		Ijv_FileNameSize,
			 CDFP_Integer32		Ijv_UnitSize)
;


#endif
