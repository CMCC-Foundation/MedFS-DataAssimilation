/*
========== SCCS configuration management ==========================

	SCCS	File name	: CDFP_PublicReadAlongTrackProduct.h
	SCCS	Version		: 1.1
	SCCS	Storage date	: 04/06/16

==================================================================

=========== RCS configuration management =========================

	RCS	File name	: $RCSfile$
	RCS	Version		: $Revision$
	RCS	Storage date	: $Date$

==================================================================

*****************************************************************************
* MODULE : CDFP_PublicReadAlongTrackProduct				 V2.0
*
* ROLE :
*    Interface to manipulate NetCdf files containing along track product (ATP)
*    I M P O R TA N T	 N O T E: This file has an equivalent for
*				  FORTRAN calls called
*				  fcc_CDFP_PublicReadAlongTrackProduct.h
*				  Any modification of fortran exported
*				  data/type must be dupplicated
*
*    NOTE: When reading ATP, data equal to the _FillValue attribute are
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
*     Thv_CDFP_ATPFile  	      Opaque type used for identifying a file
*
* PUBLIC FUNCTIONS:
*     - Opens an ATP file			CDFP_ATPOpen
*     - Close an opened file			CDFP_ATPClose
*     - Get file limits				CDFP_ATPGetMax
*     - Get list of all tracks			CDFP_ATPReadTrackList
*     - Get all cycles of a track		CDFP_ATPReadCycleList
*     - Get the list of all user variables	CDFP_ATPReadVariableList
*     - Get the number of data for a track	CDFP_ATPGetTrackNbData
*     - Get all values for a couple Cycle/Track	CDFP_ATPGetTrackDataFor1Cycle
*     - Get values for a couple Data/Track	CDFP_ATPGet1DataForCycles
*     - Retreive the netcdf file handle		CDFP_ATPGetNetCDFHandle
*     - Fortran interface (fcc=Fortran Calls C)	fcc_CDFP_ATPOpen
*						fcc_CDFP_ATPClose
*						fcc_CDFP_ATPGetMax
*						fcc_CDFP_ATPReadTrackList
*						fcc_CDFP_ATPReadCycleList
*						fcc_CDFP_ATPReadVariableList
*						fcc_CDFP_ATPGetTrackNbData
*						fcc_CDFP_ATPGetTrackDataFor1Cycle
*						fcc_CDFP_ATPGet1DataForCycles
*						fcc_CDFP_ATPGetNetCDFHandle
*
*****************************************************************************

 Modifications:
===============

 2004/01/16: Ph. Poilbarbe, SU.DM0245.
		Created
 2004/01/30: Ph. Poilbarbe.
		Check of file type was omitted.
 2004/06/16: Ph. Poilbarbe, SU.DM0354.
		Along-Track data (ATD) renamed to Along-Track Product (ATP)
		to avoid confusion with ADT files.

*/
#ifndef CDFP_PublicReadAlongTrackProduct_h
#define CDFP_PublicReadAlongTrackProduct_h


/* Identifies file version for SCCS configuration management */
#ident "@(#) CDFP_PublicReadAlongTrackProduct.h	1.1 04/06/16"

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
** Opaque type used to identify an opened file.
**
*/
typedef CDFP_Integer32	Thv_CDFP_ATPFile;








/*
==============================================================================
===                 F U N C T I O N    P R O T O T Y P E S                 ===
==============================================================================
*/










/*
******************************************************************************
**
** PUBLIC FUNCTION: CDFP_ATPOpen
**
**
** DESCRIPTION:
**    Opens a residual file for read.
**
** RETURNED VALUE:
**    File handle, or 0 if an error has occured.
**
** PARAMETERS:
**   Name	      Unit    Description
**   ------------------------------------------------------------------------
**   Itv_FileName      /       Name of the NetCdf reaiduals file
**
**
******************************************************************************
*/
Thv_CDFP_ATPFile CDFP_ATPOpen
			(CDFP_String		Itv_FileName)
;







/*
******************************************************************************
**
** PUBLIC FUNCTION: CDFP_ATPClose
**
**
** DESCRIPTION:
**    Closes a previously opened file
**
** RETURNED VALUE:
**    < none >
**
** PARAMETERS:
**   Name	      Unit    Description
**   ------------------------------------------------------------------------
**   IOhv_File		/	File handle as returned by CDFP_ATPOpen
**				If it is 0 or does not correspond to an
**				actually opened file, nothing is done.
**				Upon return, the value is 0.
**
**
******************************************************************************
*/
void CDFP_ATPClose
			(Thv_CDFP_ATPFile	*IOhv_File)
;







/*
******************************************************************************
**
** PUBLIC FUNCTION: CDFP_ATPGetMax
**
**
** DESCRIPTION:
**    Returns the maximum number of tracks and cycles by trace allowed by the
**    current file.
**
** RETURNED VALUE:
**    Error status: 0 if Ok, Other if not.
**
** PARAMETERS:
**   Name	      Unit    Description
**   ------------------------------------------------------------------------
**   Ihv_File		/	File handle as returned by CDFP_ATPOpen
**   Ojv_MaxTracks	/	Max number of track the file can contain
**   Ojv_MaxCycles	/	Max number of cycles each track can contain
**   Ojv_MaxData	/	Max actual number of data for all tracks
**   Ojv_NbVars		/	Number of user variables in the file
**
******************************************************************************
*/
CDFP_Integer32 CDFP_ATPGetMax
			(Thv_CDFP_ATPFile	Ihv_File,
			 CDFP_Integer32		*Ojv_MaxTracks,
			 CDFP_Integer32		*Ojv_MaxCycles,
			 CDFP_Integer32		*Ojv_MaxData,
			 CDFP_Integer32		*Ojv_NbVars)
;






/*
******************************************************************************
**
** PUBLIC FUNCTION: CDFP_ATPReadTrackList
**
**
** DESCRIPTION:
**    Returns the list of tracks contained in the file.
**
** RETURNED VALUE:
**    Error status: 0 if Ok, Other if not.
**
** PARAMETERS:
**   Name	      Unit    Description
**   ------------------------------------------------------------------------
**   Ihv_File		/	File handle as returned by CDFP_ATPOpen
**   Ojw_Tracks		/	Track numbers of tracks in file
**   Ojv_NbTracks	/	Actual number of tracks in file and in Ojw_Tracks
**				(0 <= Ojv_NbTracks <= Max number of tracks)
**
******************************************************************************
*/
CDFP_Integer32 CDFP_ATPReadTrackList
			(Thv_CDFP_ATPFile	Ihv_File,
			 CDFP_Integer32		*Ojw_Tracks,
			 CDFP_Integer32		*Ojv_NbTracks)
;







/*
******************************************************************************
**
** PUBLIC FUNCTION: CDFP_ATPReadCycleList
**
**
** DESCRIPTION:
**    Returns the list of cycles contained in a track.
**
** RETURNED VALUE:
**    Error status: 0 if Ok, Other if not.
**
** PARAMETERS:
**   Name	      Unit    Description
**   ------------------------------------------------------------------------
**   Ihv_File		/	File handle as returned by CDFP_ATPOpen
**   Ijv_Track		/	Track number about which we want cycle list
**   Ojw_Cycles		/	Cycles numbers associated with the track
**   Ojv_NbCycles	/	Actual number of cycles in the track
**				(0 <= Ojv_NbCycles <= Max cycles for tracks)
**
******************************************************************************
*/
CDFP_Integer32 CDFP_ATPReadCycleList
			(Thv_CDFP_ATPFile	Ihv_File,
			 CDFP_Integer32		Ijv_Track,
			 CDFP_Integer32		*Ojw_Cycles,
			 CDFP_Integer32		*Ojv_NbCycles)
;







/*
******************************************************************************
**
** PUBLIC FUNCTION: CDFP_ATPReadVariableList
**
**
** DESCRIPTION:
**    Returns the list of user variables contained in the file.
**
** RETURNED VALUE:
**    Error status: 0 if Ok, Other if not.
**
** PARAMETERS:
**   Name	      Unit    Description
**   ------------------------------------------------------------------------
**   Ihv_File		/	File handle as returned by CDFP_ATPOpen
**   Ojv_NbVars		/	Number of variables (as returned by
**				CDFP_ATPGetMax)
**   Otw_VarNames	/	Names of user variables
**				It is an array of pointers to strings
**				allocated by this routine. It is the
**				responsability of the caller to free the
**				returned array. Only one free is needed
**				to free the array and all the strings.
**				On entry the pointer passed must be set to
**				null and on exit it will be set to the
**				allocated area.
**
******************************************************************************
*/
CDFP_Integer32 CDFP_ATPReadVariableList
			(Thv_CDFP_ATPFile	Ihv_File,
			 CDFP_Integer32		*Ojv_NbVars,
			 CDFP_String		*Otw_VarNames[])
;







/*
******************************************************************************
**
** PUBLIC FUNCTION: CDFP_ATPGetTrackNbData
**
**
** DESCRIPTION:
**    Returns the number of data associated with a track
**
** RETURNED VALUE:
**    Error status: 0 if Ok, Other if not.
**
** PARAMETERS:
**   Name	      Unit    Description
**   ------------------------------------------------------------------------
**   Ihv_File		/	File handle as returned by CDFP_ATPOpen
**   Ijv_Track		/	Track number about which we want cycle list
**   Ojv_NbData		/	Number of data in track
**
******************************************************************************
*/
CDFP_Integer32 CDFP_ATPGetTrackNbData
			(Thv_CDFP_ATPFile	Ihv_File,
			 CDFP_Integer32		Ijv_Track,
			 CDFP_Integer32		*Ojv_NbData)
;







/*
******************************************************************************
**
** PUBLIC FUNCTION: CDFP_ATPGetTrackDataFor1Cycle
**
**
** DESCRIPTION:
**    Returns all the data associated with a (Track, Cycle) couple
**
** RETURNED VALUE:
**    Error status: 0 if Ok, Other if not.
**
** PARAMETERS:
**   Name	      Unit    Description
**   ------------------------------------------------------------------------
**   Ihv_File		/	File handle as returned by CDFP_ATPOpen
**   Itv_VarName	/	User variable name
**   Ijv_Track		/	Track number
**   Ijv_Cycle		/	Cycle number
**   Ojv_NbData		/	Number of data in track (as returned by
**				CDFP_ATPGetTrackNbData).
**				If NULL nothing is returned.
**   Odw_Data		/	Data for the couple (Cycle/Trace).
**				If NULL nothing is returned.
**   Odw_Latitudes	degree	Latitudes of the data point
**				If NULL nothing is returned.
**   Odw_Longitudes	degree	Longitudes of the data point
**				If NULL nothing is returned.
**   Odw_Dates		CJD	Dates of the data point.
**				If NULL nothing is returned.
**
** NOTE: All arrays (Odw_xxx) MUST be big enough to contain Ojv_NbData values.
**====== If they are set to NULL nothing is returned then allowing user to
**	 get only useful informations.
**
** CJD: Days since 1950-01-01 00:00:00.000000 UTC with precision of one
**		microsecond
**
******************************************************************************
*/
CDFP_Integer32 CDFP_ATPGetTrackDataFor1Cycle
			(Thv_CDFP_ATPFile	Ihv_File,
			 CDFP_String		Itv_VarName,
			 CDFP_Integer32		Ijv_Track,
			 CDFP_Integer32		Ijv_Cycle,
			 CDFP_Integer32		*Ojv_NbData,
			 CDFP_Real8		*Odw_Data,
			 CDFP_Real8		*Odw_Latitudes,
			 CDFP_Real8		*Odw_Longitudes,
			 CDFP_Real8		*Odw_Dates)
;







/*
******************************************************************************
**
** PUBLIC FUNCTION: CDFP_ATPGet1DataForCycles
**
**
** DESCRIPTION:
**    Returns the values associated with specified cycles of one data of a
**    track
**
** RETURNED VALUE:
**    Error status: 0 if Ok, Other if not.
**
** PARAMETERS:
**   Name	      Unit    Description
**   ------------------------------------------------------------------------
**   Ihv_File		/	File handle as returned by CDFP_ATPOpen
**   Itv_VarName	/	User variable name
**   Ijv_Track		/	Track number
**   Ijv_DataNumber	/	Index of data in track (starting at 0)
**   Ijv_NbCycles	/	Number of cycles asked.
**   Ijw_Cycles		/	Cycle list (Ijv_NbCycles)
**   Odv_Latitude	degree	Latitudes of the data point (one value)
**				If NULL nothing is returned.
**   Odv_Longitude	degree	Longitudes of the data point (one value)
**				If NULL nothing is returned.
**   Odw_Data		/	Values for each cycle of the Ijv_DataNumber-th
**				data of track Ijv_Track.
**				If NULL nothing is returned.
**   Odw_Dates		CJD	Dates of the values.
**				If NULL nothing is returned.
**
** NOTE: All arrays (Odw_xxx) MUST be big enough to contain Ojv_NbCycles values.
**====== If they are set to NULL nothing is returned then allowing user to
**	 get only useful informations.
**
** CJD: Days since 1950-01-01 00:00:00.000000 UTC with precision of one
**		microsecond
**
******************************************************************************
*/
CDFP_Integer32 CDFP_ATPGet1DataForCycles
			(Thv_CDFP_ATPFile	Ihv_File,
			 CDFP_String		Itv_VarName,
			 CDFP_Integer32		Ijv_Track,
			 CDFP_Integer32		Ijv_DataNumber,
			 CDFP_Integer32		Ijv_NbCycles,
			 CDFP_Integer32		*Ijw_Cycles,
			 CDFP_Real8		*Odv_Latitude,
			 CDFP_Real8		*Odv_Longitude,
			 CDFP_Real8		*Odw_Data,
			 CDFP_Real8		*Odw_Dates)
;







/*
******************************************************************************
**
** PUBLIC FUNCTION: CDFP_ATPGetNetCDFHandle
**
**
** DESCRIPTION:
**    Returns the handle as used by the NetCDF library (may be useful to
**    get some other informations like attributes).
**    For advanced users only (able to manipulate low level NetCDF interface)
**
** RETURNED VALUE:
**    NetCDF Handle or negative value if an error has occured
**
** PARAMETERS:
**   Name	      Unit    Description
**   ------------------------------------------------------------------------
**   Ihv_File		/	File handle as returned by CDFP_ATPOpen
**
******************************************************************************
*/
CDFP_Integer32 CDFP_ATPGetNetCDFHandle
			(Thv_CDFP_ATPFile	Ihv_File)
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
** NOTE: The prototypes for these functions in .h file are not necessary since
**       fortran do not need them. They are here as documentation.
*/









/*
******************************************************************************
**
** PUBLIC FUNCTION: fcc_CDFP_ATPOpen
**
**
** DESCRIPTION:
**    Opens a residual file for read.
**
** RETURNED VALUE:
**    File handle, or 0 if an error has occured.
**
** PARAMETERS:
**   Name	      Unit    Description
**   ------------------------------------------------------------------------
**   Itv_FileName      /       Name of the NetCdf reaiduals file
**   Ohv_FileHandle    /       File handle or 0 if an error has occured
**
** FORTRAN EQUIVALENT DECLARATION:
**  SUBROUTINE fcc_CDFP_ATPOpen
**  &              (Itv_FileName,
**  &               Ohv_FileHandle)
**    CHARACTER*(*) Itv_FileName
**    INTEGER*4     Ohv_FileHandle
**
******************************************************************************
*/
void FTN_NAME(fcc_cdfp_atpopen)
			(CDFP_String		Itv_FileName,
			 CDFP_Integer32		*Ohv_FileHandle,


/*******************Automatically added by fortran do not use them as actual parameters */
			 CDFP_Integer32		Ijv_FileNameSize)
;







/*
******************************************************************************
**
** PUBLIC FUNCTION: fcc_CDFP_ATPClose
**
**
** DESCRIPTION:
**    Closes a previously opened file
**
** RETURNED VALUE:
**    < none >
**
** PARAMETERS:
**   Name	      Unit    Description
**   ------------------------------------------------------------------------
**   IOhv_File		/	File handle as returned by CDFP_ATPOpen
**				If it is 0 or does not correspond to an
**				actually opened file, nothing is done.
**				Upon return, the value is 0.
**
** FORTRAN EQUIVALENT DECLARATION:
**  SUBROUTINE fcc_CDFP_ATPClose
**  &              (IOhv_File)
**    INTEGER*4     IOhv_File
**
******************************************************************************
*/
void FTN_NAME(fcc_cdfp_atpclose)
			(CDFP_Integer32	*IOhv_File)
;







/*
******************************************************************************
**
** PUBLIC FUNCTION: fcc_CDFP_ATPGetMax
**
**
** DESCRIPTION:
**    Returns the maximum number of tracks and cycles by trace allowed by the
**    current file.
**
** RETURNED VALUE:
**    <none>
**
** PARAMETERS:
**   Name	      Unit    Description
**   ------------------------------------------------------------------------
**   Ihv_File		/	File handle as returned by CDFP_ATPOpen
**   Ojv_MaxTracks	/	Max number of track the file can contain
**   Ojv_MaxCycles	/	Max number of cycles each track can contain
**   Ojv_MaxData	/	Max actual number of data for all tracks
**   Ojv_NbVars		/	Number of user variables in the file
**   Ojv_Status		/	0 if Ok other if not
**
** FORTRAN EQUIVALENT DECLARATION:
**  SUBROUTINE fcc_CDFP_ATPGetMax
**  &              (Ihv_File,
**  &		    Ojv_MaxTracks,
**  &		    Ojv_MaxCycles,
**  &		    Ojv_MaxData,
**  &		    Ojv_NbVars,
**  &		    Ojv_Status)
**    INTEGER*4     Ihv_File
**    INTEGER*4     Ojv_MaxTracks
**    INTEGER*4     Ojv_MaxCycles
**    INTEGER*4     Ojv_MaxData
**    INTEGER*4     Ojv_NbVars
**    INTEGER*4     Ojv_Status
**
******************************************************************************
*/
void FTN_NAME(fcc_cdfp_atpgetmax)
			(CDFP_Integer32		*Ihv_File,
			 CDFP_Integer32		*Ojv_MaxTracks,
			 CDFP_Integer32		*Ojv_MaxCycles,
			 CDFP_Integer32		*Ojv_MaxData,
			 CDFP_Integer32		*Ojv_NbVars,
			 CDFP_Integer32		*Ojv_Status)
;






/*
******************************************************************************
**
** PUBLIC FUNCTION: fcc_CDFP_ATPReadTrackList
**
**
** DESCRIPTION:
**    Returns the list of tracks contained in the file.
**
** RETURNED VALUE:
**    <none>
**
** PARAMETERS:
**   Name	      Unit    Description
**   ------------------------------------------------------------------------
**   Ihv_File		/	File handle as returned by CDFP_ATPOpen
**   Ojw_Tracks		/	Track numbers of tracks in file
**   Ojv_NbTracks	/	Actual number of tracks in file and in Ojw_Tracks
**				(0 <= Ojv_NbTracks <= Max number of tracks)
**   Ojv_Status		/	0 if Ok other if not
**
** FORTRAN EQUIVALENT DECLARATION:
**  SUBROUTINE fcc_CDFP_ATPReadTrackList
**  &              (Ihv_File,
**  &		    Ojw_Tracks,
**  &		    Ojv_NbTracks,
**  &		    Ojv_Status)
**    INTEGER*4     Ihv_File
**    INTEGER*4     Ojw_Tracks(*)
**    INTEGER*4     Ojv_NbTracks
**    INTEGER*4     Ojv_Status
**
******************************************************************************
*/
void FTN_NAME(fcc_cdfp_atpreadtracklist)
			(CDFP_Integer32		*Ihv_File,
			 CDFP_Integer32		*Ojw_Tracks,
			 CDFP_Integer32		*Ojv_NbTracks,
			 CDFP_Integer32		*Ojv_Status)
;







/*
******************************************************************************
**
** PUBLIC FUNCTION: fcc_CDFP_ATPReadCycleList
**
**
** DESCRIPTION:
**    Returns the list of cycles contained in a track.
**
** RETURNED VALUE:
**    <none>
**
** PARAMETERS:
**   Name	      Unit    Description
**   ------------------------------------------------------------------------
**   Ihv_File		/	File handle as returned by CDFP_ATPOpen
**   Ijv_Track		/	Track number about which we want cycle list
**   Ojw_Cycles		/	Cycles numbers associated with the track
**   Ojv_NbCycles	/	Actual number of cycles in the track
**				(0 <= Ojv_NbCycles <= Max cycles for tracks)
**   Ojv_Status		/	0 if Ok other if not
**
** FORTRAN EQUIVALENT DECLARATION:
**  SUBROUTINE fcc_CDFP_ATPReadCycleList
**  &              (Ihv_File,
**  &		    Ijv_Track,
**  &		    Ojw_Cycles,
**  &		    Ojv_NbCycles,
**  &		    Ojv_Status)
**    INTEGER*4     Ihv_File
**    INTEGER*4     Ijv_Track
**    INTEGER*4     Ojw_Cycles(*)
**    INTEGER*4     Ojv_NbCycles
**    INTEGER*4     Ojv_Status
**
******************************************************************************
*/
void FTN_NAME(fcc_cdfp_atpreadcyclelist)
			(CDFP_Integer32		*Ihv_File,
			 CDFP_Integer32		*Ijv_Track,
			 CDFP_Integer32		*Ojw_Cycles,
			 CDFP_Integer32		*Ojv_NbCycles,
			 CDFP_Integer32		*Ojv_Status)
;







/*
******************************************************************************
**
** PUBLIC FUNCTION: fcc_CDFP_ATPReadVariableList
**
**
** DESCRIPTION:
**    Returns the list of user variables contained in the file.
**
** RETURNED VALUE:
**    Error status: 0 if Ok, Other if not.
**
** PARAMETERS:
**   Name	      Unit    Description
**   ------------------------------------------------------------------------
**   Ihv_File		/	File handle as returned by CDFP_ATPOpen
**   Ojv_NbVars		/	Number of variables (as returned by
**				CDFP_ATPGetMax)
**   Otw_VarNames	/	Names of user variables
**   Ijv_NbMaxVars	/	Number of names Otw_VarNames cn contain
**				If there is more than this number of var in the
**				file, Otw_VarNames if fullfilled and last names
**				are lost. Ojv_NbVars can then be tested to see
**				if we are in this case
**   Ojv_Status		/	0 if Ok other if not
**
** FORTRAN EQUIVALENT DECLARATION:
**  SUBROUTINE fcc_CDFP_ATPReadVariableList
**  &              (Ihv_File,
**  &		    Ojv_NbVars,
**  &		    Otw_VarNames,
**  &		    Ijv_NbMaxVars
**  &		    Ojv_Status)
**    INTEGER*4     Ihv_File
**    INTEGER*4     Ojv_NbVars
**    CHARACTER*(*) Ojw_VarNames(Ijv_NbMaxVars)
**    INTEGER*4     Ijv_NbMaxVars
**    INTEGER*4     Ojv_Status
**
******************************************************************************
*/
void FTN_NAME(fcc_cdfp_atpreadvariablelist)
			(CDFP_Integer32		*Ihv_File,
			 CDFP_Integer32		*Ojv_NbVars,
			 CDFP_String		Otw_VarNames,
			 CDFP_Integer32		*Ijv_NbMaxVars,
			 CDFP_Integer32		*Ojv_Status,


/*******************Automatically added by fortran do not use them as actual parameters */
			 CDFP_Integer32		Ijv_VarNamesSize)
;







/*
******************************************************************************
**
** PUBLIC FUNCTION: fcc_CDFP_ATPGetTrackNbData
**
**
** DESCRIPTION:
**    Returns the number of data associated with a track
**
** RETURNED VALUE:
**    <none>
**
** PARAMETERS:
**   Name	      Unit    Description
**   ------------------------------------------------------------------------
**   Ihv_File		/	File handle as returned by CDFP_ATPOpen
**   Ijv_Track		/	Track number about which we want cycle list
**   Ojv_NbData		/	Number of data in track
**   Ojv_Status		/	0 if Ok other if not
**
** FORTRAN EQUIVALENT DECLARATION:
**  SUBROUTINE fcc_CDFP_ATPGetTrackNbData
**  &              (Ihv_File,
**  &		    Ijv_Track,
**  &		    Ojv_NbData,
**  &		    Ojv_Status)
**    INTEGER*4     Ihv_File
**    INTEGER*4     Ijv_Track
**    INTEGER*4     Ojv_NbData
**    INTEGER*4     Ojv_Status
**
******************************************************************************
*/
void FTN_NAME(fcc_cdfp_atpgettracknbdata)
			(CDFP_Integer32		*Ihv_File,
			 CDFP_Integer32		*Ijv_Track,
			 CDFP_Integer32		*Ojv_NbData,
			 CDFP_Integer32		*Ojv_Status)
;







/*
******************************************************************************
**
** PUBLIC FUNCTION: fcc_CDFP_ATPGetTrackDataFor1Cycle
**
**
** DESCRIPTION:
**    Returns all the data associated with a (Track, Cycle) couple
**
** RETURNED VALUE:
**    <none>
**
** PARAMETERS:
**   Name	      Unit    Description
**   ------------------------------------------------------------------------
**   Ihv_File		/	File handle as returned by CDFP_ATPOpen
**   Itv_VarName	/	User variable name
**   Ijv_Track		/	Track number
**   Ijv_Cycle		/	Cycle number
**   Ojv_NbData		/	Number of data in track (as returned by
**				fcc_CDFP_ATPGetTrackNbData).
**   Odw_Data		/	Data for the couple (Cycle/Trace).
**   Odw_Latitudes	degree	Latitudes of the data point
**   Odw_Longitudes	degree	Longitudes of the data point
**   Odw_Dates		CJD	Dates of the data point.
**   Ojv_Status		/	0 if Ok other if not
**
** FORTRAN EQUIVALENT DECLARATION:
**  SUBROUTINE fcc_CDFP_ATPGetTrackDataFor1Cycle
**  &              (Ihv_File,
**  &		    Itv_VarName,
**  &		    Ijv_Track,
**  &		    Ijv_Cycle,
**  &		    Ojv_NbData,
**  &		    Odw_Data,
**  &		    Odw_Latitudes,
**  &		    Odw_Longitudes,
**  &		    Odw_Dates,
**  &		    Ojv_Status)
**    INTEGER*4     Ihv_File
**    CHARACTER*(*) Itv_VarName
**    INTEGER*4     Ijv_Track
**    INTEGER*4     Ijv_Cycle
**    INTEGER*4     Ojv_NbData
**    REAL*8        Odw_Data(*)
**    REAL*8        Odw_Latitudes(*)
**    REAL*8        Odw_Longitudes(*)
**    REAL*8        Odw_Dates(*)
**    INTEGER*4     Ojv_Status
**
** NOTE: All arrays (Odw_xxx) MUST be big enough to contain Ojv_NbData values.
**======
**
** CJD: Days since 1950-01-01 00:00:00.000000 UTC with precision of one
**		microsecond
**
******************************************************************************
*/
void FTN_NAME(fcc_cdfp_atpgettrackdatafor1cycle)
			(CDFP_Integer32		*Ihv_File,
			 CDFP_String		Itv_VarName,
			 CDFP_Integer32		*Ijv_Track,
			 CDFP_Integer32		*Ijv_Cycle,
			 CDFP_Integer32		*Ojv_NbData,
			 CDFP_Real8		*Odw_Data,
			 CDFP_Real8		*Odw_Latitudes,
			 CDFP_Real8		*Odw_Longitudes,
			 CDFP_Real8		*Odw_Dates,
			 CDFP_Integer32		*Ojv_Status,


/*******************Automatically added by fortran do not use them as actual parameters */
			 CDFP_Integer32		Ijv_VarNameSize)
;







/*
******************************************************************************
**
** PUBLIC FUNCTION: fcc_CDFP_ATPGet1DataForCycles
**
**
** DESCRIPTION:
**    Returns the values associated with specified cycles of one data of a
**    track
**
** RETURNED VALUE:
**    <none>
**
** PARAMETERS:
**   Name	      Unit    Description
**   ------------------------------------------------------------------------
**   Ihv_File		/	File handle as returned by CDFP_ATPOpen
**   Itv_VarName	/	User variable name
**   Ijv_Track		/	Track number
**   Ijv_DataNumber	/	Index of data in track (starting at 0)
**   Ijv_NbCycles	/	Number of cycles asked.
**   Ijw_Cycles		/	Cycle list (Ijv_NbCycles)
**   Odv_Latitude	degree	Latitudes of the data point (one value)
**   Odv_Longitude	degree	Longitudes of the data point (one value)
**   Odw_Data		/	Values for each cycle of the Ijv_DataNumber-th
**				data of track Ijv_Track.
**   Odw_Dates		CJD	Dates of the values.
**   Ojv_Status		/	0 if Ok other if not
**
** FORTRAN EQUIVALENT DECLARATION:
**  SUBROUTINE fcc_CDFP_ATPGet1DataForCycles
**  &              (Ihv_File,
**  &		    Itv_VarName,
**  &		    Ijv_Track,
**  &		    Ijv_DataNumber,
**  &		    Ijv_NbCycles,
**  &		    Ijw_Cycles,
**  &		    Odv_Latitude,
**  &		    Odv_Longitude,
**  &		    Odw_Data,
**  &		    Odw_Dates,
**  &		    Ojv_Status)
**    INTEGER*4     Ihv_File
**    CHARACTER*(*) Itv_VarName
**    INTEGER*4     Ijv_Track
**    INTEGER*4     Ijv_DataNumber
**    INTEGER*4     Ijv_NbCycles
**    INTEGER*4     Ijw_Cycles(Ijv_NbCycles)
**    REAL*8        Odv_Latitude
**    REAL*8        Odv_Longitude
**    REAL*8        Odw_Data(Ijv_NbCycles)
**    REAL*8        Odw_Dates(Ijv_NbCycles)
**    INTEGER*4     Ojv_Status
**
** NOTE: All arrays (Odw_xxx) MUST be big enough to contain Ojv_NbCycles values.
**====== If they are set to NULL nothing is returned then allowing user to
**	 get only useful informations.
**
** CJD: Days since 1950-01-01 00:00:00.000000 UTC with precision of one
**		microsecond
**
******************************************************************************
*/
void FTN_NAME(fcc_cdfp_atpget1dataforcycles)
			(CDFP_Integer32		*Ihv_File,
			 CDFP_String		Itv_VarName,
			 CDFP_Integer32		*Ijv_Track,
			 CDFP_Integer32		*Ijv_DataNumber,
			 CDFP_Integer32		*Ijv_NbCycles,
			 CDFP_Integer32		*Ijw_Cycles,
			 CDFP_Real8		*Odv_Latitude,
			 CDFP_Real8		*Odv_Longitude,
			 CDFP_Real8		*Odw_Data,
			 CDFP_Real8		*Odw_Dates,
			 CDFP_Integer32		*Ojv_Status,


/*******************Automatically added by fortran do not use them as actual parameters */
			 CDFP_Integer32		Ijv_VarNameSize)
;







#endif
