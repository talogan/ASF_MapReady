/*==============================================================================
Filename:	ADEOSfname.c

Description:
	This module contains the function(s) used for determining if a
filename is a ADEOS type file name - if it conforms to ADEOS's file 
naming convention.

External Functions:
	parse_ADEOS_filename
	
Static Functions:
	None
	
External Variables Defined:
	ADEOS_File_Id_Table

File Scope Static Variables:
	None
	
Notes:

==============================================================================*/

static char SccsFile[] = "%M%" ;
static char SccsRevision[] = "%R%" ;
static char SccsDate[] = "%G%";
static char SccsLastChanger[] = "%W%";
static char SccsState[] = "%I%";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <syslog.h>
#include "faifdefs.h"
#include "ADEOS.h"


#ifdef __STDC__
int parse_ADEOS_filename(char *) ;
#else
int parse_ADEOS_filename() ;
#endif

extern void *util_do_malloc() ;
 
File_Identifier ADEOS_File_Id_Table[] =
{
   { ADEOS_REQR, ADEOS_REQR_HSTR }, /* REQR */
   { ADEOS_RDRD, ADEOS_RDRD_HSTR }, /* RDRD */
   { ADEOS_ELMD, ADEOS_ELMD_HSTR }, /* ELMD */
   { ADEOS_ELMP, ADEOS_ELMP_HSTR }, /* ELMP */
   { ADEOS_OPL1, ADEOS_OPL1_HSTR }, /* OPL1 */
   { ADEOS_RPLN, ADEOS_RPLN_HSTR }, /* RPLN */
   { ADEOS_ORST, ADEOS_ORST_HSTR }, /* ORST */
   { ADEOS_STAD, ADEOS_STAD_HSTR }, /* STAD */
   { ADEOS_TMDF, ADEOS_TMDF_HSTR }, /* TMDF */
   { ADEOS_STGS, ADEOS_STGS_HSTR }, /* STGS */
   { ADEOS_REAC, ADEOS_REAC_HSTR }, /* REAC */
   { ADEOS_SRRD, ADEOS_SRRD_HSTR }, /* Shipment Report */
   { SENTINEL,   NULL            }
} ;




/*==============================================================================
Function:	int parse_ADEOS_filename(filename)

Description:
	Parse and validate the input string filename as a ADEOS
filename.  This is an initial means of determining the file type of a
specific flight agency file.  If filename conforms to the ADEOS file
naming convention, then it is further validated as an ADEOS file.  The
table ADEOS_File_Id_Table is consulted to validate file id information
that are obtained from the filename.  The return status of this
function is the file id code for a valid file type or ERROR if the
filename is invalid.

Parameters:
	filename - Filename part of a file's filepath specification.
Note that the names used must be valid Unix filenames.

Returns:	valid ADEOS file type or ERROR
Creator:	Norbert Piega	
Creation Date:	10/04/1995
Notes:		
==============================================================================*/
#ifdef __STDC__
int
parse_ADEOS_filename(char *filename)
#else
int
parse_ADEOS_filename(filename)
   char *filename ;
#endif
{
   int i, type_matched = ERROR ;
   char file_id[ADEOS_ID_STRLEN+1] ;
   char *start, *fname ;
   char logmsg[MAX_SYSLOG_MSGLEN+1] ;

   /* Precaution.  Skip past the last '/' in the file specification
   */
   if ((start = strrchr(filename, '/')) != NULL)
   {
      fname = (char *)util_do_malloc(sizeof(char)*(strlen(start+1)+1)) ;
      strcpy(fname, start+1) ;
   }
   else
      fname = filename ;

   if ((int)strlen(fname) < ADEOS_ID_STRLEN)
   {
      sprintf(logmsg, "WARNING, %s is an invalid ADEOS file name\n", fname) ;
      syslog(LOG_ERR, logmsg) ;
      return(ERROR) ;
   }

   /* Compare file id from filename to table of
   -- file types and associated ids
   */
   file_id[0] = toupper(fname[0]) ;
   file_id[1] = toupper(fname[1]) ;
   file_id[2] = toupper(fname[2]) ;
   file_id[3] = toupper(fname[3]) ;
   file_id[4] = '\0' ;
   i = 0 ;
   while (ADEOS_File_Id_Table[i].file_id_number != SENTINEL)
      if (strcmp(ADEOS_File_Id_Table[i].file_identifier, file_id) == 0)
      {
         type_matched = ADEOS_File_Id_Table[i].file_id_number ;
         break ; 
      }
      else
         i++ ;

   return(type_matched) ;

} /* parse_ADEOS_filename */


/* End of File */
