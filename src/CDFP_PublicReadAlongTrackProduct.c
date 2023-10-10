/*
========== SCCS configuration management ==========================

	SCCS	File name	: CDFP_PublicReadAlongTrackProduct.c
	SCCS	Version		: 1.2
	SCCS	Storage date	: 06/06/09

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
*				  fcc_CDFP_PublicReadAlongTrackData.h
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
 2006/06/09: Ph. Poilbarbe. SU.FA0502. CDFP_ATPGetTrackDataFor1Cycle, misspelling
		of variable name when computing Odw_Dates (Odw_Data was used instead
		which may produce an invalid date or a core dump).
		Thanks to Kaoru Ichikawa indentifying and reporting the bug.

*/

/* Identifies file version for SCCS configuration management */
#ident "@(#) CDFP_PublicReadAlongTrackProduct.c	1.2 06/06/09"

/* Identifies file version for CVS/RCS configuration management */
#ident "$Id$"

#include <assert.h>
#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <netcdf.h>

#include "CDFP_PublicReadAlongTrackProduct.h"





/*
==============================================================================
===            I N T E R N A L    T Y P E S    A N D    D A T A            ===
==============================================================================
*/


/*
** Macros used to do some repetitive checks
*/
#define CHECK_NULL(x)	if ((x) == NULL) goto Error
#define IS_OPENED(x)	(((x) > 0) &&				\
			 ((x) <= jc_MaxATPFilesOpened) &&	\
			 (sb_Files[x] != NULL))
#define CHECK_OPENED(x)	if (! IS_OPENED(x))			\
			{					\
			  ErrorMsg("Try tu use file #%d while it is not opened");	\
			  goto Error;				\
			}


/*
** Max function
*/
#define MAX(X,Y)	((X) > (Y) ? (X) : (Y))

/*
**
** Attribute names known by NetCdf packages
**
*/
#define tc_UnitAttr		"units"
#define tc_ScaleFactorAttr	"scale_factor"
#define tc_AddOffsetAttr	"add_offset"
#define tc_FillValueAttr	"_FillValue"

/*
**
** ATP specific names
**
*/
#define tc_TracksDim		"Tracks"
#define tc_CyclesDim		"Cycles"
#define tc_FileTypeAttr		"FileType"
#define tc_DeltaTVar		"DeltaT"
#define tc_TracksVar		"Tracks"
#define tc_NbPointsVar		"NbPoints"
#define tc_CyclesVar		"Cycles"
#define tc_LongitudesVar	"Longitudes"
#define tc_LatitudesVar		"Latitudes"
#define tc_BeginDatesVar	"BeginDates"
#define tc_DataIndexesVar	"DataIndexes"
#define tc_GlobalCycleListVar	"GlobalCyclesList"

/*
**
** Maximum number of opened files at one time
**
*/
#define jc_MaxATPFilesOpened	 20

/*
**
** File kind (actually one)
**
*/
typedef enum {
  E_ATPWithMeanProfile
} Tnv_FileType;


/*
**
** Opened file informations cached for performances
**
*/
typedef struct {
	CDFP_String	tv_FileName;
	CDFP_Integer32	jv_File;
	Tnv_FileType	nv_FileType;
	CDFP_Integer32	jv_MaxTracks;
	CDFP_Integer32	jv_MaxCycles;
        CDFP_Integer32	jv_MaxDataInTrack;	/* Computed */
	CDFP_Integer32	jv_NbUserVars;		/* Computed */
	CDFP_String	*sq_UserVarNames;
	CDFP_Real8	dv_DeltaT;
	CDFP_Integer32	*jw_Tracks;	/* Tracks actually in the file (jv_MaxTracks) */
	CDFP_Integer32	*jw_NbPoints;	/* Number of points in each track (jv_MaxTracks) */
	CDFP_Integer32	*jw_Cycles;	/* Cycles in each track (jv_MaxTracks*jv_MaxCycles) */
	CDFP_Real8	*dw_BeginDates;	/* Cycles in each track (jv_MaxTracks*jv_MaxCycles) */
	CDFP_Real8	*dw_Data;	/* Longitudes of measurement (Max(jv_MaxDataInTrack,jv_MaxTracks)) */
} Tsv_FileDescription;

/*
**
**
** Notice that index 0 is never used (File id 0 is always a closed file
** in order to allow vars initialized to 0).
**
*/


/*
**
** List of opened/closed files
** Notice that index 0 is never used (File id 0 is always a closed file
** in order to allow vars initialized to 0).
**
*/
static Tsv_FileDescription	*sb_Files[jc_MaxATPFilesOpened+1];


/*
**
** Current function name
**
*/
static CDFP_String		tv_FunctionName;



/*
==============================================================================
==============================================================================
==============================================================================
==============================================================================


==============================================================================
===                  I N T E R N A L    F U N C T I O N S                  ===
==============================================================================
*/






/*
******************************************************************************
**
** PRIVATE FUNCTION: QSortInteger32Callback
**
**
** DESCRIPTION:
**    Callback to QSort used to sort 32 bits integers in ascending order
**
**
** RETURNED VALUE:
**	Boolean: <0, 0 or >0 according to qsort(3)
**
** PARAMETERS:
**   Name	      Unit    Description
**   ------------------------------------------------------------------------
**   Ijv_X		/	First number to compare
**   Ijv_Y		/	Second number to compare
**
******************************************************************************
*/
static int QSortInteger32Callback
			(const void	*Ijv_X,
			 const void	*Ijv_Y)
{
  return *((CDFP_Integer32 *)Ijv_X) - *((CDFP_Integer32 *)Ijv_Y);
}






/*
******************************************************************************
**
** PRIVATE FUNCTION: IsUserVarName
**
**
** DESCRIPTION:
**    Indicates if a variable is a user varname or not
**
**
** RETURNED VALUE:
**	Boolean: 0 (False)    if not a user var name
**		 Other (true) if it is a user var name
**
** PARAMETERS:
**   Name	      Unit    Description
**   ------------------------------------------------------------------------
**   Itv_VarName       /       Var name to check
**
******************************************************************************
*/
static CDFP_Integer32 IsUserVarName
			(const CDFP_String	VarName)
{
  return (strcmp(VarName, tc_DeltaTVar)			!= 0)	&&
         (strcmp(VarName, tc_TracksVar)			!= 0)	&&
         (strcmp(VarName, tc_NbPointsVar)		!= 0)	&&
         (strcmp(VarName, tc_CyclesVar)			!= 0)	&&
         (strcmp(VarName, tc_BeginDatesVar)		!= 0)	&&
	 (strcmp(VarName, tc_DataIndexesVar)		!= 0)	&&
         (strcmp(VarName, tc_LatitudesVar)		!= 0)	&&
         (strcmp(VarName, tc_LongitudesVar)		!= 0)	&&
	 (strcmp(VarName, tc_GlobalCycleListVar)	!= 0);
}





/*
******************************************************************************
**
** INTERNAL FUNCTION: ErrorMsg
**
**
** DESCRIPTION:
**    Shows an error message
**
** RETURNED VALUE:
**    <none>
**
** PARAMETERS:
**   Name              Unit		Description
**   ------------------------------------------------------------------------
**   Itv_Format		/		Format of message in a printf manner
**   ...		/		Parameters of message
**
**
******************************************************************************
*/
static void ErrorMsg
			(const CDFP_String	Itv_Message,
			 ...)
{
  va_list		ArgP;

  errno	= 0;
  fprintf(stderr,"\n************* ERROR in %s **********\n", tv_FunctionName);

  va_start(ArgP, Itv_Message);
  vfprintf(stderr, Itv_Message, ArgP);
  va_end(ArgP);

  fprintf(stderr,"\n");
}





/*
******************************************************************************
**
** INTERNAL FUNCTION: NetCdfErrorCheck
**
**
** DESCRIPTION:
**    Shows an error message if a NetCdf error occurs
**
** RETURNED VALUE:
**    True if error, False if not
**
** PARAMETERS:
**   Name              Unit		Description
**   ------------------------------------------------------------------------
**   Ijv_CdfError	/		NetCdfError Code
**   Itv_Format		/		Format of message in a printf manner
**   ...		/		Parameters of message
**
**
******************************************************************************
*/
static CDFP_Integer32 NetCdfErrorCheck
			(const CDFP_Integer32	Ijv_NetCdfError,
			 const CDFP_String	Itv_Message,
			 ...)
{
  va_list		ArgP;

  if (Ijv_NetCdfError == NC_NOERR)
    return 0;

  fprintf(stderr,"\n************* ERROR **********\n");

  va_start(ArgP, Itv_Message);
  vfprintf(stderr, Itv_Message, ArgP);
  va_end(ArgP);

  fprintf(stderr,"\n");

  fprintf(stderr, "NetCdf error message: %s\n", nc_strerror(Ijv_NetCdfError));
  return 1;
}




/*
******************************************************************************
**
** INTERNAL FUNCTION: New
**
**
** DESCRIPTION:
**    Allocates memory and writes an error message if it cannot
**
** RETURNED VALUE:
**    Pointer to allocated memory or NULL if an error as occured
**
** PARAMETERS:
**   Name		Unit	Description
**   ------------------------------------------------------------------------
**   Ijv_Size		Bytes	Number of bytes to allocate.
**   Itv_Name		/	Name of file involved or NULL
**
**
******************************************************************************
*/
static void *New
		(CDFP_Integer32	Ijv_Size,
		 CDFP_String	Itv_Name)
{
  void *Tmp;
  Tmp	= calloc(1, Ijv_Size);
  if (Tmp == NULL)
  {
    ErrorMsg("Not enough memory to allocate %d bytes%s%s%s",
		Ijv_Size,
		(Itv_Name == NULL ? "" : " (For file "),
		(Itv_Name == NULL ? "" : Itv_Name),
		(Itv_Name == NULL ? "" : ")")
	);
  }
  return Tmp;
}




/*
******************************************************************************
**
** INTERNAL FUNCTION: FtnStrDup
**
**
** DESCRIPTION:
**    Converts a fortran string to a C string (NUL terminated)
**
** RETURNED VALUE:
**    C string (allocated, must be freed)
**
** PARAMETERS:
**   Name	      Unit	      Description
**   ------------------------------------------------------------------------
**   Iaw_FString       / 	      Fortran string
**   Ijv_Length        / 	      Fortran string length
**
**
******************************************************************************
*/
static CDFP_String FtnStrDup
			(char		*Iaw_FString,
			 CDFP_Integer32	 Ijv_Length)
{

  CDFP_String		 Result;
  char			*Find;
  CDFP_Integer32	 Length;

  Length	= (Ijv_Length > 0 ? Ijv_Length : 0);

  Result	= calloc(1, Length+1);

  if (Result != NULL)
  {
    /* Copy the whole fortran string */
    memcpy(Result, Iaw_FString, Length);


    /* Find last non blank character */
    Find = Result+Length-1;
    while ((Find >= Result) && (*Find == ' '))
    {
      Find--;
    }

    /* Mark the end of the string */
    *(Find+1) = '\0';
  }
  else
    ErrorMsg("Not enougn memory to allocate %d characters to convert fortran string to C",
		Length+1);
  return Result;
}




/*
******************************************************************************
**
** INTERNAL FUNCTION: Str2Ftn
**
**
** DESCRIPTION:
**    Converts a C string (NUL terminated) to a fortran string
**
** RETURNED VALUE:
**    < none >
**
** PARAMETERS:
**   Name	      Unit	      Description
**   ------------------------------------------------------------------------
**   Iaw_FString       / 	      Fortran string
**   Ijv_Length        / 	      Fortran string length
**   Itv_String        / 	      C String to convert
**
**
******************************************************************************
*/
static void Str2Ftn
                (char		*Oaw_FString,
                 CDFP_Integer32	Ijv_Length,
                 CDFP_String	Itv_String)
{
  CDFP_Integer32     Length;

  memset(Oaw_FString, ' ', Ijv_Length);
  if ((Itv_String != NULL) && (*Itv_String != '\0'))
  {
    Length      = strlen(Itv_String);
    if (Length>Ijv_Length)
      Length    = Ijv_Length;

    memcpy(Oaw_FString, Itv_String, Length);
  }

  return;
}





/*
******************************************************************************
**
** INTERNAL FUNCTION: GetFloatAttribute
**
**
** DESCRIPTION:
**    Retreive the value of a numérical attribute as a float. Returns a default
**    value if attribute does not exist
**
** RETURNED VALUE:
**    True if error, False if not
**
** PARAMETERS:
**   Name              Unit		Description
**   ------------------------------------------------------------------------
**   Ijv_NetCdfHandle	/		Netcdf file handle
**   Ijv_VarId		/		Netcdf var id of var to retreive
**   Itv_Attribute	/		Attribute name
**   Odv_Value		/		Value retreived
**   Obv_Exists		/		Indicates if attribute exists
**   Itv_VarName	/		Variable name (for error messages)
**
**
******************************************************************************
*/
static CDFP_Integer32 GetFloatAttribute
			(const CDFP_Integer32	Ijv_NetCdfHandle,
			 const CDFP_Integer32	Ijv_VarId,
			 const CDFP_String	Itv_Attribute,
			 CDFP_Real8		*Odv_Value,
			 CDFP_Integer32		*Obv_Exists,
			 const CDFP_String	Itv_VarName)
{
  CDFP_Integer32	CdfStatus;
  size_t		Length;
  CDFP_Integer32	Exists	= 1;

  CdfStatus	= nc_inq_attlen(Ijv_NetCdfHandle, Ijv_VarId, Itv_Attribute, &Length);

  if (CdfStatus == NC_ENOTATT)
  {
    *Odv_Value	= dc_CDFP_DefReal8;
    Exists	= 0;
  }
  else
  {
    if (NetCdfErrorCheck(CdfStatus, "Getting attribute %s:%s",
			 Itv_VarName,
			 Itv_Attribute))
      goto Error;

    if (Length != 1)
    {
      ErrorMsg("Attribute %s:%s length %d is invalid, only 1 is valid",
			 Itv_VarName,
			 Itv_Attribute,
			 Length);
      goto Error;
    }
    CdfStatus	= nc_get_att_double(Ijv_NetCdfHandle,
				    Ijv_VarId,
				    Itv_Attribute,
				    Odv_Value);
    if (NetCdfErrorCheck(CdfStatus, "Getting attribute %s:%s",
			 Itv_VarName,
			 Itv_Attribute))
      goto Error;
  }

  *Obv_Exists	= Exists;
  return 0;

/*********************************************Error treatment */
Error:
  *Obv_Exists	= 0;
  *Odv_Value	= dc_CDFP_DefReal8;
  return 1;
}





/*
******************************************************************************
**
** INTERNAL FUNCTION: FreeFile
**
**
** DESCRIPTION:
**    Frees a file descriptor created by NewFile
**
** RETURNED VALUE:
**	<none>
**
** PARAMETERS:
**   Name		Unit		Description
**   ------------------------------------------------------------------------
**   Ijv_FileId		/		Index of opened file
**
******************************************************************************
*/
static void FreeFile
			(const CDFP_Integer32	Ijv_FileId)
{
  if (IS_OPENED(Ijv_FileId))
  {
    if (sb_Files[Ijv_FileId]->jv_File != 0)
      nc_close(sb_Files[Ijv_FileId]->jv_File);
    free(sb_Files[Ijv_FileId]->tv_FileName);
    free(sb_Files[Ijv_FileId]->jw_Tracks);
    free(sb_Files[Ijv_FileId]->jw_NbPoints);
    free(sb_Files[Ijv_FileId]->jw_Cycles);
    free(sb_Files[Ijv_FileId]->dw_BeginDates);
    free(sb_Files[Ijv_FileId]->sq_UserVarNames);
    free(sb_Files[Ijv_FileId]);
    sb_Files[Ijv_FileId]	= NULL;
  }
}





/*
******************************************************************************
**
** INTERNAL FUNCTION: NewFile
**
**
** DESCRIPTION:
**    Allocates a new file handle at spécified position
**
** RETURNED VALUE:
**    True if error, False if not
**
** PARAMETERS:
**   Name              Unit		Description
**   ------------------------------------------------------------------------
**   Ijv_FileId		/		Index of file to be opened
**   Itv_FileName	/		Name of the file to be opened
**
******************************************************************************
*/
static CDFP_Integer32 NewFile
			(const CDFP_Integer32	Ijv_FileId,
			 const CDFP_String	Itv_FileName)
{
  assert(sb_Files[Ijv_FileId] == NULL);

  sb_Files[Ijv_FileId]	= New(sizeof(*sb_Files[Ijv_FileId]), Itv_FileName);
  CHECK_NULL(sb_Files[Ijv_FileId]);

  sb_Files[Ijv_FileId]->tv_FileName	= strdup(Itv_FileName);
  if (sb_Files[Ijv_FileId]->tv_FileName == NULL)
  {
    ErrorMsg("Not enough memory to store file name '%s'", Itv_FileName);
    goto Error;
  }

  return 0;

/*********************************************Error treatment */
Error:
  FreeFile(Ijv_FileId);
  return 1;
}





/*
******************************************************************************
**
** INTERNAL FUNCTION: ReadIntVar
**
**
** DESCRIPTION:
**	Returns the value of an integer variable. Attributes scale_factor
**	and add_offset are not applyed since this is for control variables
**	and they should not have such attributes.
**
** RETURNED VALUE:
**    True if error, False if not
**
** PARAMETERS:
**   Name              Unit		Description
**   ------------------------------------------------------------------------
**   Ijv_NetCdfHandle	/		Netcdf file handle
**   Ojw_Values		/		Values read
**   Itv_VarName	/		Netcdf var name
**   Ijv_StartIndex1	/		Starting index along first dimension
**					if <0 a 0 dim var is expected
**   Ijv_Count1		/		Number of values along first dimension
**   Ijv_StartIndex2	/		Starting index along second dimension
**					if <0 1 dim var is expected
**   Ijv_Count2		/		Number of values along second dimension
**
**
******************************************************************************
*/
static CDFP_Integer32 ReadIntVar
			(const CDFP_Integer32	Ijv_CdfFileHandle,
			 CDFP_Integer32		*Ojw_Values,
			 const CDFP_String	Itv_VarName,
			 CDFP_Integer32		Ijv_StartIndex1,
			 CDFP_Integer32		Ijv_Count1,
			 CDFP_Integer32		Ijv_StartIndex2,
			 CDFP_Integer32		Ijv_Count2)
{

  size_t		Starts[2];
  size_t		Counts[2];
  CDFP_Integer32	CdfStatus;
  CDFP_Integer32	VarId;		/* NetCdf Id of a variable */
  nc_type		Type;		/* Physical storage of a var */
  CDFP_Integer32	Exists;		/* Indicates if an attribute exists */
  CDFP_Integer32	ExpectedDims;
  CDFP_Integer32	NbDims;
  CDFP_Real8		FillValueFlt;
  CDFP_Integer32	FillValue;	/* FillValueInt converted into int */

  CDFP_Integer32	Index;


  ExpectedDims	= (Ijv_StartIndex1 >= 0 ? (Ijv_StartIndex2 >= 0 ? 2 : 1) : 0);

  CdfStatus	= nc_inq_varid(Ijv_CdfFileHandle, Itv_VarName, &VarId);
  if (NetCdfErrorCheck(CdfStatus, "Getting %s variable id", Itv_VarName))
    goto Error;

  CdfStatus	= nc_inq_vartype(Ijv_CdfFileHandle, VarId, &Type);
  if (NetCdfErrorCheck(CdfStatus, "Getting %s variable type", Itv_VarName))
    goto Error;

  if (Type != NC_INT)
  {
    ErrorMsg("Variable %s is not an integer variable", Itv_VarName);
    goto Error;
  }

  CdfStatus	= nc_inq_varndims(Ijv_CdfFileHandle, VarId, &NbDims);
  if (NetCdfErrorCheck(CdfStatus, "Getting %s variable dimensions", Itv_VarName))
    goto Error;

  if (NbDims != ExpectedDims)
  {
    ErrorMsg("Variable %s has %d dimension while %d are expected",
		NbDims,
		ExpectedDims);
    goto Error;
  }

  Starts[0]	= Ijv_StartIndex1;
  Starts[1]	= Ijv_StartIndex2;
  Counts[0]	= (NbDims >= 1 ? Ijv_Count1 : 1); /* 1 are sentinnels */
  Counts[1]	= (NbDims >= 2 ? Ijv_Count2 : 1);

  /*
  ** Read the values
  */
  if (ExpectedDims == 0)
    CdfStatus	= nc_get_var_int(Ijv_CdfFileHandle, VarId, Ojw_Values);
  else
    CdfStatus	= nc_get_vara_int(Ijv_CdfFileHandle, VarId, Starts, Counts, Ojw_Values);
  if (NetCdfErrorCheck(CdfStatus, "Getting %s variable values", Itv_VarName))
    goto Error;


  /*
  ** Take into account the _FillValue attribute and replace it by a consistent
  ** one (jc_CDFP_DefInteger32)
  */
  if (GetFloatAttribute(Ijv_CdfFileHandle,
			VarId,
			tc_FillValueAttr,
			&FillValueFlt,
			&Exists,
			Itv_VarName))
    goto Error;

  if (! Exists)
    FillValue	= NC_FILL_INT; /* Default netcdf fill value */
  else
    FillValue	= FillValueFlt;

  for (Index=0; Index<(Ijv_Count1*(Ijv_StartIndex2 >= 0 ? Ijv_Count2 : 1)); Index++)
  {
    if (Ojw_Values[Index] == FillValue)
      Ojw_Values[Index]	= jc_CDFP_DefInteger32;
  }

  return 0;

/*********************************************Error treatment */
Error:
  return 1;
}





/*
******************************************************************************
**
** INTERNAL FUNCTION: ReadFloatVar
**
**
** DESCRIPTION:
**	Returns the value of an integer variable. Attributes _FillValue,
**	add_offset and scale_factor are applyed.
**
** RETURNED VALUE:
**    True if error, False if not
**
** PARAMETERS:
**   Name              Unit		Description
**   ------------------------------------------------------------------------
**   Ijv_NetCdfHandle	/		Netcdf file handle
**   Odw_Values		/		Values read
**   Itv_VarName	/		Netcdf var name
**   Ijv_StartIndex1	/		Starting index along first dimension
**					if <0 a 0 dim var is expected
**   Ijv_Count1		/		Number of values along first dimension
**   Ijv_StartIndex2	/		Starting index along second dimension
**					if <0 1 dim var is expected
**   Ijv_Count2		/		Number of values along second dimension
**
**
******************************************************************************
*/
static CDFP_Integer32 ReadFloatVar
			(const CDFP_Integer32	Ijv_CdfFileHandle,
			 CDFP_Real8		*Odw_Values,
			 const CDFP_String	Itv_VarName,
			 CDFP_Integer32		Ijv_StartIndex1,
			 CDFP_Integer32		Ijv_Count1,
			 CDFP_Integer32		Ijv_StartIndex2,
			 CDFP_Integer32		Ijv_Count2)
{

  size_t		Starts[2];
  size_t		Counts[2];
  CDFP_Integer32	CdfStatus;
  CDFP_Integer32	VarId;		/* NetCdf Id of a variable */
  nc_type		Type;		/* Physical storage of a var */
  CDFP_Integer32	Exists;		/* Indicates if an attribute exists */
  CDFP_Integer32	ExpectedDims;
  CDFP_Integer32	NbDims;
  CDFP_Real8		FillValue;
  CDFP_Real8		ScaleFactor;
  CDFP_Real8		AddOffset;

  CDFP_Integer32	NbData;
  CDFP_Real8		*Current;


  ExpectedDims	= (Ijv_StartIndex1 >= 0 ? (Ijv_StartIndex2 >= 0 ? 2 : 1) : 0);

  CdfStatus	= nc_inq_varid(Ijv_CdfFileHandle, Itv_VarName, &VarId);
  if (NetCdfErrorCheck(CdfStatus, "Getting %s variable id", Itv_VarName))
    goto Error;

  CdfStatus	= nc_inq_varndims(Ijv_CdfFileHandle, VarId, &NbDims);
  if (NetCdfErrorCheck(CdfStatus, "Getting %s variable dimensions", Itv_VarName))
    goto Error;

  if (NbDims != ExpectedDims)
  {
    ErrorMsg("Variable %s has %d dimension while %d are expected",
		NbDims,
		ExpectedDims);
    goto Error;
  }

  Starts[0]	= Ijv_StartIndex1;
  Starts[1]	= Ijv_StartIndex2;
  Counts[0]	= (NbDims >= 1 ? Ijv_Count1 : 1); /* 1 are sentinnels */
  Counts[1]	= (NbDims >= 2 ? Ijv_Count2 : 1);

  /*
  ** Read the values
  */
  if (ExpectedDims == 0)
    CdfStatus	= nc_get_var_double(Ijv_CdfFileHandle, VarId, Odw_Values);
  else
    CdfStatus	= nc_get_vara_double(Ijv_CdfFileHandle, VarId, Starts, Counts, Odw_Values);
  if (NetCdfErrorCheck(CdfStatus, "Getting %s variable values", Itv_VarName))
    goto Error;


  /*
  ** Take into account the _FillValue attribute and replace it by a consistent
  ** one (jc_CDFP_DefInteger32)
  */
  if (GetFloatAttribute(Ijv_CdfFileHandle,
			VarId,
			tc_FillValueAttr,
			&FillValue,
			&Exists,
			Itv_VarName))
    goto Error;

  if (! Exists)
  {/* Find default netcdf fill value */
    CdfStatus	= nc_inq_vartype(Ijv_CdfFileHandle, VarId, &Type);
    if (NetCdfErrorCheck(CdfStatus, "Getting %s variable type", Itv_VarName))
      goto Error;
    switch (Type)
    {
      case NC_BYTE:	FillValue	= NC_FILL_BYTE;		break;
      case NC_SHORT:	FillValue	= NC_FILL_SHORT;	break;
      case NC_INT:	FillValue	= NC_FILL_INT;		break;
      case NC_FLOAT:	FillValue	= NC_FILL_FLOAT;	break;
      case NC_DOUBLE:	FillValue	= NC_FILL_DOUBLE;	break;
      default:
	ErrorMsg("Invalid or unknown type %d for variable %s", Type, Itv_VarName);
        goto Error;
    }
  }


  /*
  **
  ** Get conversion factor attributes
  **
  */
  if (GetFloatAttribute(Ijv_CdfFileHandle,
			VarId,
			tc_ScaleFactorAttr,
			&ScaleFactor,
			&Exists,
			Itv_VarName))
    goto Error;

  if (GetFloatAttribute(Ijv_CdfFileHandle,
			VarId,
			tc_AddOffsetAttr,
			&AddOffset,
			&Exists,
			Itv_VarName))
    goto Error;

  NbData	= Counts[0]*Counts[1];
  Current	= Odw_Values;
  while (NbData-- >0)
  {
    if (*Current == FillValue)
      *Current = dc_CDFP_DefReal8;
    else
    {
      if (ScaleFactor != dc_CDFP_DefReal8)
	*Current	*= ScaleFactor;

      if (AddOffset != dc_CDFP_DefReal8)
	*Current	+= AddOffset;
    }

    Current++;
  }

  return 0;

/*********************************************Error treatment */
Error:
  return 1;
}





/*
******************************************************************************
**
** INTERNAL FUNCTION: FindTrack
**
**
** DESCRIPTION:
**	Returns the index of a track known by its number. Returns also
**	the position of first data of the track.
**
** RETURNED VALUE:
**    True if error, False if not
**
** PARAMETERS:
**   Name              Unit		Description
**   ------------------------------------------------------------------------
**   Isv_File		/		File descriptor
**   Ijv_Track		/		Number (id) of track to find
**   Ojv_Index		/		Index of track in track list
**   Ojv_DataIndex	/		Index of first data of track in
**					data arrays. If NULL nothing is
**					returned.
**
******************************************************************************
*/
static CDFP_Integer32 FindTrack
			(Tsv_FileDescription	*Isv_File,
			 CDFP_Integer32		Ijv_Track,
			 CDFP_Integer32		*Ojv_Index,
			 CDFP_Integer32		*Ojv_DataIndex)
{
  CDFP_Integer32	Index;
  CDFP_Integer32	DataIndex;

  if ((Ijv_Track < 0) || (Ijv_Track == jc_CDFP_DefInteger32))
  {
    ErrorMsg("Invalid track number %d for %s", Ijv_Track, Isv_File->tv_FileName);
    goto Error;
  }


  DataIndex	= 0;
  for (Index=0; Index<Isv_File->jv_MaxTracks; Index++)
  {
    if (Isv_File->jw_Tracks[Index] == Ijv_Track)
    {
      *Ojv_Index	= Index;
      if (Ojv_DataIndex != NULL)
        *Ojv_DataIndex	= DataIndex;
      return 0;
    }
    if (Isv_File->jw_NbPoints[Index] != jc_CDFP_DefInteger32)
      DataIndex	+= Isv_File->jw_NbPoints[Index];
  }

  ErrorMsg("Track #%d not in file %s", Ijv_Track, Isv_File->tv_FileName);

/*********************************************Error treatment */
Error:
  return 1;
}





/*
******************************************************************************
**
** INTERNAL FUNCTION: FindCycle
**
**
** DESCRIPTION:
**	Returns the index within a track of a cycle known by its number.
**
** RETURNED VALUE:
**    True if error, False if not
**
** PARAMETERS:
**   Name              Unit		Description
**   ------------------------------------------------------------------------
**   Isv_File		/		File descriptor
**   Ijv_TrackIndex	/		Index (position) of track
**   Ijv_Cycle		/		Cycle number to find
**   Ojv_Index		/		Index of cycle in track
**
******************************************************************************
*/
static CDFP_Integer32 FindCycle
			(Tsv_FileDescription	*Isv_File,
			 CDFP_Integer32		Ijv_TrackIndex,
			 CDFP_Integer32		Ijv_Cycle,
			 CDFP_Integer32		*Ojv_Index)
{
  CDFP_Integer32	Index;
  CDFP_Integer32	*CycleList;

  assert(Ijv_TrackIndex >= 0 && Ijv_TrackIndex < Isv_File->jv_MaxTracks);

  if ((Ijv_Cycle < 0) || (Ijv_Cycle == jc_CDFP_DefInteger32))
  {
    ErrorMsg("Invalid cycle number %d for %s", Ijv_Cycle, Isv_File->tv_FileName);
    goto Error;
  }

  CycleList	= Isv_File->jw_Cycles +
		  (Ijv_TrackIndex * Isv_File->jv_MaxCycles);

  for (Index=0; Index < Isv_File->jv_MaxCycles; Index++)
  {
    if (CycleList[Index] == Ijv_Cycle)
    {
      *Ojv_Index	= Index;
      return 0;
    }
  }

  ErrorMsg("Cycle/Track #%d/%d not in file %s",
		Ijv_Cycle,
		Isv_File->jw_Tracks[Ijv_TrackIndex],
		Isv_File->tv_FileName);

/*********************************************Error treatment */
Error:
  return 1;
}




/*
==============================================================================
==============================================================================
==============================================================================
==============================================================================


==============================================================================
===                    P U B L I C    F U N C T I O N S                    ===
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
{
  CDFP_Integer32	FileId;
  Tsv_FileDescription	*File;
  CDFP_Integer32	CdfFileHandle	= 0;
  CDFP_Integer32	CdfStatus;
  CDFP_Integer32	NbCdfVars;

  CDFP_Integer32	DimId;		/* NetCdf Id of a dimension */
  CDFP_Integer32	Index;
  nc_type		Type;		/* attribute/var type */
  size_t		TmpVal;
  char			Buffer[50];	/* To read strings */



  VERIFY_SIZES;

  tv_FunctionName	= "CDFP_ATPOpen";
  /*
  ** Find an unused entry in opened file list
  */
  FileId	= 1;
  while (FileId <= jc_MaxATPFilesOpened)
  {
    if (sb_Files[FileId] == NULL)
      break;
    FileId++;
  }
  if (FileId > jc_MaxATPFilesOpened)
  {
    ErrorMsg("Too much opened files (Max %d at one time)", jc_MaxATPFilesOpened);
    goto Error;
  }

  if (NewFile(FileId, Itv_FileName))
    goto Error;

  File	= sb_Files[FileId];

  CdfStatus	= nc_open(Itv_FileName, NC_NOWRITE, &CdfFileHandle);

  if (NetCdfErrorCheck(CdfStatus, "Opening file %s", Itv_FileName))
    goto Error;

  /*
  **
  ** Get file type
  **
  */
  CdfStatus = nc_inq_att(CdfFileHandle,
			 NC_GLOBAL,
			 tc_FileTypeAttr,
			 &Type,
			 &TmpVal);

  if (NetCdfErrorCheck(CdfStatus, "Getting attribute %s info", tc_FileTypeAttr))
    goto Error;

  if ((Type != NC_CHAR) || (TmpVal >= sizeof(Buffer)))
  {
    ErrorMsg("Attribute %s is not a string or is invalid", tc_FileTypeAttr);
    goto Error;
  }

  memset(Buffer, 0, sizeof(Buffer));
  CdfStatus	= nc_get_att_text(CdfFileHandle, NC_GLOBAL, tc_FileTypeAttr, Buffer);
  if (NetCdfErrorCheck(CdfStatus, "Getting attribute %s", tc_FileTypeAttr))
    goto Error;

  if ((strcmp(Buffer, "ALONG_TRACK_PRODUCT")  == 0) ||
      (strcmp(Buffer, "ALONG_TRACK_DATA")     == 0) ||	/* very Old name */
      (strcmp(Buffer, "SEA_LEVEL_ANOMALY")    == 0))	/* very Old name */
  {
    File->nv_FileType	= E_ATPWithMeanProfile;
  }
  else
  {
    ErrorMsg("File type is '%s' while '%s' is expected, invalid file %s",
		Buffer,
		"ALONG_TRACK_PRODUCT",
		Itv_FileName);
    goto Error;
  }


  /*
  ** Get the maximum number of tracks the file can contain
  */
  CdfStatus	= nc_inq_dimid(CdfFileHandle, tc_TracksDim, &DimId);
  if (NetCdfErrorCheck(CdfStatus, "Getting dimension %s id", tc_TracksDim))
    goto Error;

  CdfStatus	= nc_inq_dimlen(CdfFileHandle, DimId, &TmpVal);
  if (NetCdfErrorCheck(CdfStatus, "Getting dimension %s", tc_TracksDim))
    goto Error;
  File->jv_MaxTracks	= TmpVal;

  /*
  ** Get the maximum number of cycles per track the file can contain
  */
  CdfStatus	= nc_inq_dimid(CdfFileHandle, tc_CyclesDim, &DimId);
  if (NetCdfErrorCheck(CdfStatus, "Getting dimension %s id", tc_CyclesDim))
    goto Error;

  CdfStatus	= nc_inq_dimlen(CdfFileHandle, DimId, &TmpVal);
  if (NetCdfErrorCheck(CdfStatus, "Getting dimension %s", tc_CyclesDim))
    goto Error;
  File->jv_MaxCycles	= TmpVal;

  /*
  ** Allocates some data (the ones for which the sizes are known)
  */
  File->jw_Tracks	= New(sizeof(*File->jw_Tracks)*File->jv_MaxTracks, Itv_FileName);
  CHECK_NULL(File->jw_Tracks);
  File->jw_NbPoints	= New(sizeof(*File->jw_NbPoints)*File->jv_MaxTracks, Itv_FileName);
  CHECK_NULL(File->jw_NbPoints);
  File->jw_Cycles	= New(sizeof(*File->jw_Cycles)*
				File->jv_MaxTracks*
				File->jv_MaxCycles,
			      Itv_FileName);
  CHECK_NULL(File->jw_Cycles);
  File->dw_BeginDates	= New(sizeof(*File->dw_BeginDates)*
				File->jv_MaxTracks*
				File->jv_MaxCycles,
			      Itv_FileName);
  CHECK_NULL(File->dw_BeginDates);

  /*
  ** Read the corresponding data
  */
  if (ReadIntVar(CdfFileHandle, File->jw_Tracks, tc_TracksVar, 0, File->jv_MaxTracks, -1, -1))
    goto Error;

  if (ReadIntVar(CdfFileHandle, File->jw_NbPoints, tc_NbPointsVar, 0, File->jv_MaxTracks, -1, -1))
    goto Error;

  if (ReadIntVar(CdfFileHandle, File->jw_Cycles, tc_CyclesVar, 0, File->jv_MaxTracks, 0, File->jv_MaxCycles))
    goto Error;

  if (ReadFloatVar(CdfFileHandle, File->dw_BeginDates, tc_BeginDatesVar, 0, File->jv_MaxTracks, 0, File->jv_MaxCycles))
    goto Error;

  if (ReadFloatVar(CdfFileHandle, &File->dv_DeltaT, tc_DeltaTVar, -1, -1, -1, -1))
    goto Error;

  /*
  ** Computes the actual max number of points in a track
  */
  File->jv_MaxDataInTrack	= 1;
  for (Index=0; Index<File->jv_MaxTracks; Index++)
  {
    if ((File->jw_NbPoints[Index] != jc_CDFP_DefInteger32) &&
	(File->jw_NbPoints[Index] > File->jv_MaxDataInTrack))
      File->jv_MaxDataInTrack	= File->jw_NbPoints[Index];
  }

  /*
  ** Check user variables and allocates memory to store the names
  */
  CdfStatus	= nc_inq_nvars(CdfFileHandle, &NbCdfVars);
  if (NetCdfErrorCheck(CdfStatus, "Getting number of variables"))
    goto Error;


    /*
    ** Array of pointer to strings and strings themselves allocated at one
    ** time and then freed using a unique call to free
    */
  File->sq_UserVarNames	= New((sizeof(*File->sq_UserVarNames)+NC_MAX_NAME+1)*NbCdfVars,
			      Itv_FileName);
  CHECK_NULL(File->sq_UserVarNames);

  File->jv_NbUserVars	= 0;
  for (Index=0; Index<NbCdfVars; Index++)
  {
    File->sq_UserVarNames[File->jv_NbUserVars]	= ((char *)File->sq_UserVarNames) +
						  (sizeof(*File->sq_UserVarNames)*NbCdfVars) +
						  (File->jv_NbUserVars*(NC_MAX_NAME+1));
    CdfStatus	= nc_inq_varname(CdfFileHandle,
				 Index,
				 File->sq_UserVarNames[File->jv_NbUserVars]);
    if (NetCdfErrorCheck(CdfStatus, "Getting name of variable #%d", Index))
      goto Error;

    if (IsUserVarName(File->sq_UserVarNames[File->jv_NbUserVars]))
      File->jv_NbUserVars++;
  }


/*
**
** End function
**
*/
  sb_Files[FileId]->jv_File	= CdfFileHandle;
  return FileId;


/*********************************************Error treatment */
Error:
/*
**
** Free allocated resources
**
*/
  if (CdfFileHandle != 0)
    nc_close(CdfFileHandle);
  FreeFile(FileId);
  return 0;
}







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
{
  tv_FunctionName	= "CDFP_ATPClose";

  FreeFile(*IOhv_File);
  *IOhv_File	= 0;
}







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
{
  tv_FunctionName	= "CDFP_ATPGetMax";

  CHECK_OPENED(Ihv_File);

  *Ojv_MaxTracks	= sb_Files[Ihv_File]->jv_MaxTracks;
  *Ojv_MaxCycles	= sb_Files[Ihv_File]->jv_MaxCycles;
  *Ojv_MaxData		= sb_Files[Ihv_File]->jv_MaxDataInTrack;
  *Ojv_NbVars		= sb_Files[Ihv_File]->jv_NbUserVars;

  return 0;


/*********************************************Error treatment */
Error:
  return 1;
}






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
{
  Tsv_FileDescription	*File;
  CDFP_Integer32	Filled;
  CDFP_Integer32	Index;

  tv_FunctionName	= "CDFP_ATPReadTrackList";
  CHECK_OPENED(Ihv_File);

  File		= sb_Files[Ihv_File];
  Filled	= 0;
  for (Index=0; Index < File->jv_MaxTracks; Index++)
  {
    if (File->jw_Tracks[Index] != jc_CDFP_DefInteger32)
      Ojw_Tracks[Filled++]	= File->jw_Tracks[Index];
  }

  *Ojv_NbTracks	= Filled;
  /*
  ** Make sure track numbers are in ascending number (there is not specific
  ** order in file)
  */
  qsort(Ojw_Tracks, Filled, sizeof(*Ojw_Tracks), QSortInteger32Callback);

  return 0;


/*********************************************Error treatment */
Error:
  return 1;
}







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
{
  Tsv_FileDescription	*File;
  CDFP_Integer32	Filled;
  CDFP_Integer32	Index;
  CDFP_Integer32	TrackIndex;
  CDFP_Integer32	*TrackCycles;

  tv_FunctionName	= "CDFP_ATPReadCycleList";
  CHECK_OPENED(Ihv_File);

  File		= sb_Files[Ihv_File];
  if (FindTrack(File, Ijv_Track, &TrackIndex, NULL))
    goto Error;

  Filled	= 0;
  TrackCycles	= File->jw_Cycles + TrackIndex*File->jv_MaxCycles;
  for (Index=0; Index < File->jv_MaxCycles; Index++)
  {
    if (TrackCycles[Index] != jc_CDFP_DefInteger32)
      Ojw_Cycles[Filled++]	= TrackCycles[Index];
  }


  *Ojv_NbCycles	= Filled;
  /*
  ** Make sure cycle numbers are in ascending number (there is not specific
  ** order in file)
  */
  qsort(Ojw_Cycles, Filled, sizeof(*Ojw_Cycles), QSortInteger32Callback);

  return 0;


/*********************************************Error treatment */
Error:
  return 1;
}







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
{
  Tsv_FileDescription	*File;
  CDFP_String		*Result		= NULL;
  CDFP_Integer32	Index;

  tv_FunctionName	= "CDFP_ATPReadVariableList";
  CHECK_OPENED(Ihv_File);
  if (*Otw_VarNames != NULL)
  {
    ErrorMsg("PROGRAM ERROR: Pointer to be returned is not NULL upon entry. Not freed ?");
    goto Error;
  }

  File		= sb_Files[Ihv_File];
  Result	= New(File->jv_NbUserVars*(sizeof(*Result)+(NC_MAX_NAME+1)),
		      File->tv_FileName);
  for (Index=0; Index < File->jv_NbUserVars; Index++)
  {
    Result[Index]	= ((char *)Result)+
			   (sizeof(*Result)*File->jv_NbUserVars) +
			   (Index * (NC_MAX_NAME+1));
    strncpy(Result[Index], File->sq_UserVarNames[Index], NC_MAX_NAME);
  }
  *Ojv_NbVars	= File->jv_NbUserVars;
  *Otw_VarNames	= Result;

  return 0;


/*********************************************Error treatment */
Error:
  free(Result);
  return 1;
}







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
{
  Tsv_FileDescription	*File;
  CDFP_Integer32	TrackIndex;

  tv_FunctionName	= "CDFP_ATPGetTrackNbData";

  CHECK_OPENED(Ihv_File);

  File		= sb_Files[Ihv_File];
  if (FindTrack(File, Ijv_Track, &TrackIndex, NULL))
    goto Error;

  *Ojv_NbData	= File->jw_NbPoints[TrackIndex];

  return 0;


/*********************************************Error treatment */
Error:
  return 1;
}







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
{
  Tsv_FileDescription	*File;
  CDFP_Integer32	TrackIndex;
  CDFP_Integer32	CycleIndex;
  CDFP_Integer32	DataPosition;
  CDFP_Integer32	NbData;
  CDFP_Real8		BeginDate;
  CDFP_Integer32	Index;


  tv_FunctionName	= "CDFP_ATPGetTrackDataFor1Cycle";
  CHECK_OPENED(Ihv_File);

  File		= sb_Files[Ihv_File];
  if (FindTrack(File, Ijv_Track, &TrackIndex, &DataPosition))
    goto Error;

  if (FindCycle(File, TrackIndex, Ijv_Cycle, &CycleIndex))
    goto Error;

  NbData	= File->jw_NbPoints[TrackIndex];

  if (Odw_Data != NULL)
  {
    if (ReadFloatVar(File->jv_File,
		     Odw_Data,
		     Itv_VarName,
		     DataPosition,
		     NbData,
		     CycleIndex,
		     1))
      goto Error;
  }

  if (Odw_Latitudes != NULL)
  {
    if (ReadFloatVar(File->jv_File,
		     Odw_Latitudes,
		     tc_LatitudesVar,
		     DataPosition,
		     NbData,
		     -1,
		     -1))
      goto Error;
  }

  if (Odw_Longitudes != NULL)
  {
    if (ReadFloatVar(File->jv_File,
		     Odw_Longitudes,
		     tc_LongitudesVar,
		     DataPosition,
		     NbData,
		     -1,
		     -1))
      goto Error;
  }

  if (Odw_Dates != NULL)
  {
    if (ReadFloatVar(File->jv_File,
		     Odw_Dates,
		     tc_DataIndexesVar,
		     DataPosition,
		     NbData,
		     -1,
		     -1))
      goto Error;
    BeginDate	= File->dw_BeginDates[TrackIndex*File->jv_MaxCycles + CycleIndex];
    for (Index=0; Index<NbData; Index++)
    {
      Odw_Dates[Index]	= BeginDate + (Odw_Dates[Index]*File->dv_DeltaT)/86400.0;
    }
  }

  if (Ojv_NbData != NULL)
  {
    *Ojv_NbData	= NbData;
  }
  return 0;


/*********************************************Error treatment */
Error:
  return 1;
}







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
{
  Tsv_FileDescription	*File;
  CDFP_Integer32	TrackIndex;
  CDFP_Integer32	CycleIndex;
  CDFP_Integer32	DataPosition;
  CDFP_Integer32	NbData;
  CDFP_Real8		*BeginDates;
  CDFP_Real8		*Values		= NULL;
  CDFP_Real8		DeltaDate;
  CDFP_Integer32	Index;


  tv_FunctionName	= "CDFP_ATPGet1DataForCycles";
  CHECK_OPENED(Ihv_File);

  File		= sb_Files[Ihv_File];
  if (FindTrack(File, Ijv_Track, &TrackIndex, &DataPosition))
    goto Error;


  NbData	= File->jw_NbPoints[TrackIndex];

  if ((Ijv_DataNumber < 0) || (Ijv_DataNumber >= NbData))
  {
    ErrorMsg("Invalid data number #%d for track #%d in file %s (Expected [0, %d]",
		Ijv_DataNumber,
		Ijv_Track,
		File->tv_FileName,
		NbData-1);
    goto Error;
  }


  if (ReadFloatVar(File->jv_File,
		   &DeltaDate,
		   tc_DataIndexesVar,
		   DataPosition+Ijv_DataNumber,
		   1,
		   -1,
		   -1))
    goto Error;

  DeltaDate	= (DeltaDate*File->dv_DeltaT)/86400.0;

  if (Odw_Data != NULL)
  {
    Values	= New(sizeof(*Values)*File->jv_MaxCycles, File->tv_FileName);
    CHECK_NULL(Values);

    if (ReadFloatVar(File->jv_File,
		     Values,
		     Itv_VarName,
		     DataPosition+Ijv_DataNumber,
		     1,
		     0,
		     File->jv_MaxCycles))
      goto Error;
  }


  BeginDates	= File->dw_BeginDates + TrackIndex*File->jv_MaxCycles;

  for (Index=0; Index<Ijv_NbCycles; Index++)
  {
    if (FindCycle(File, TrackIndex, Ijw_Cycles[Index], &CycleIndex))
      goto Error;

    if (Odw_Data != NULL)
    {
      Odw_Data[Index]	= Values[CycleIndex];
    }

    if (Odw_Dates != NULL)
    {
      Odw_Dates[Index]	= BeginDates[CycleIndex] + DeltaDate;
    }
  }

  if (Odv_Latitude != NULL)
  {
    if (ReadFloatVar(File->jv_File,
		     Odv_Latitude,
		     tc_LatitudesVar,
		     DataPosition+Ijv_DataNumber,
		     1,
		     -1,
		     -1))
      goto Error;
  }

  if (Odv_Longitude != NULL)
  {
    if (ReadFloatVar(File->jv_File,
		     Odv_Longitude,
		     tc_LongitudesVar,
		     DataPosition+Ijv_DataNumber,
		     1,
		     -1,
		     -1))
      goto Error;
  }

  free(Values);
  return 0;


/*********************************************Error treatment */
Error:
  free(Values);
  return 1;
}







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
{
  tv_FunctionName	= "CDFP_ATPGetNetCDFHandle";

  CHECK_OPENED(Ihv_File);

  return sb_Files[Ihv_File]->jv_File;

/*********************************************Error treatment */
Error:
  return -1;
}








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
{
  CDFP_String		FileName	= NULL;

  tv_FunctionName	= "fcc_CDFP_ATPOpen";
  FileName		= FtnStrDup(Itv_FileName, Ijv_FileNameSize);
  CHECK_NULL(FileName);

  *Ohv_FileHandle	= CDFP_ATPOpen(FileName);

  goto EndFunction;


/*********************************************Error treatment */
Error:
  *Ohv_FileHandle	= 0;

/*********************************************Common part */
EndFunction:
/*
**
** Free allocated resources
**
*/
  free(FileName);

  return;
}







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
{
  CDFP_ATPClose(IOhv_File);
}







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
{

  if (CDFP_ATPGetMax(*Ihv_File,
		     Ojv_MaxTracks,
		     Ojv_MaxCycles,
		     Ojv_MaxData,
		     Ojv_NbVars))
    goto Error;

  *Ojv_Status	= 0;
  goto EndFunction;


/*********************************************Error treatment */
Error:
  *Ojv_Status	= 1;

/*********************************************Common part */
EndFunction:
  return;
}






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
{

  if (CDFP_ATPReadTrackList(*Ihv_File,
			    Ojw_Tracks,
			    Ojv_NbTracks))
    goto Error;

  *Ojv_Status	= 0;
  goto EndFunction;


/*********************************************Error treatment */
Error:
  *Ojv_Status	= 1;

/*********************************************Common part */
EndFunction:
  return;
}







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
{

  if (CDFP_ATPReadCycleList(*Ihv_File,
			    *Ijv_Track,
			    Ojw_Cycles,
			    Ojv_NbCycles))
    goto Error;

  *Ojv_Status	= 0;
  goto EndFunction;


/*********************************************Error treatment */
Error:
  *Ojv_Status	= 1;

/*********************************************Common part */
EndFunction:
  return;
}







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
{
  CDFP_Integer32	Index;
  CDFP_String		*VarList	= NULL;

  if (CDFP_ATPReadVariableList(*Ihv_File,
			       Ojv_NbVars,
			       &VarList))
    goto Error;

  for (Index=0;
       Index < (*Ojv_NbVars < *Ijv_NbMaxVars ? *Ojv_NbVars : *Ijv_NbMaxVars);
       Index++)
  {
    Str2Ftn(Otw_VarNames + (Index*Ijv_VarNamesSize),
	    Ijv_VarNamesSize,
	    VarList[Index]);
  }

  *Ojv_Status	= 0;
  goto EndFunction;


/*********************************************Error treatment */
Error:
  *Ojv_Status	= 1;

/*********************************************Common part */
EndFunction:
/*
**
** Free allocated resources
**
*/
  free(VarList);
  return;
}







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
{

  if (CDFP_ATPGetTrackNbData(*Ihv_File, *Ijv_Track, Ojv_NbData))
    goto Error;

  *Ojv_Status	= 0;
  goto EndFunction;


/*********************************************Error treatment */
Error:
  *Ojv_Status	= 1;

/*********************************************Common part */
EndFunction:
  return;
}







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
/* void FTN_NAME(fcc_cdfp_atpgettrackdatafor1cycle) */
void FTN_NAME(fcc_cdfp_atpgtdf1c)
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
{
  CDFP_String		VarName;

/* AS  tv_FunctionName	= "fcc_CDFP_ATPGetTrackDataFor1Cycle"; */
  tv_FunctionName	= "fcc_CDFP_ATPGTDF1C";


  VarName		= FtnStrDup(Itv_VarName, Ijv_VarNameSize);
  CHECK_NULL(VarName);

  if (CDFP_ATPGetTrackDataFor1Cycle(*Ihv_File,
				    VarName,
				    *Ijv_Track,
				    *Ijv_Cycle,
				    Ojv_NbData,
				    Odw_Data,
				    Odw_Latitudes,
				    Odw_Longitudes,
				    Odw_Dates))
    goto Error;


  *Ojv_Status	= 0;
  goto EndFunction;


/*********************************************Error treatment */
Error:
  *Ojv_Status	= 1;

/*********************************************Common part */
EndFunction:
/*
**
** Free allocated resources
**
*/
  free(VarName);
  return;
}







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
{
  CDFP_String		VarName;

  tv_FunctionName	= "fcc_CDFP_ATPGet1DataForCycles";


  VarName		= FtnStrDup(Itv_VarName, Ijv_VarNameSize);
  CHECK_NULL(VarName);

  if (CDFP_ATPGet1DataForCycles(*Ihv_File,
				VarName,
				*Ijv_Track,
				*Ijv_DataNumber,
				*Ijv_NbCycles,
				Ijw_Cycles,
				Odv_Latitude,
				Odv_Longitude,
				Odw_Data,
				Odw_Dates))
    goto Error;


  *Ojv_Status	= 0;
  goto EndFunction;


/*********************************************Error treatment */
Error:
  *Ojv_Status	= 1;

/*********************************************Common part */
EndFunction:
/*
**
** Free allocated resources
**
*/
  free(VarName);
  return;
}


