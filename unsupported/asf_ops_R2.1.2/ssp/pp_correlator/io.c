/* SccsId[]= @(#)io.c	2.41 3/24/98 */
static char sccsid_io[]= "@(#)PPio.c:2.41";

/* SccsId[ ]= @(#)io.c	1.1 2/5/92 */
#include "pp.h"

#define PERMS 0666
/* IO library:
 * written by quyen dinh nguyen
 * 11/12/91.   
 *
 * Modified by Shelby Yang to utilize info_handler to handle errors.
 * DO not mistake this file as being the same as other io.c file scattered
 * in other places.  
 */

/* To open a file and assign a channel to it. This must be
   done before any attempt is made to access the file. The 
   return value (initdk) is the file descriptor. The file can
   be closed with the closedk subroutine.
   
   Remember, always open files before you need to access to them
   and close them after you don't need them any more. In UNIX,
   there is a limit (20) of number files can be opened at once.

   Note that, if the file is not existing, the routine will create
   the new file  with PERM=0666.

   Calling sequence(from FORTRAN):
         fd = initdk(lun,filename)
   where:
         fd is the long int for file descriptor.

         lun is the dummy variable to be compatible with VMS calls.

         filename is the name of the file. Include directory paths 
         if necessary.
 */
int initdk(lun, filename)
int *lun; char *filename;
{  int i;
   int fd;
   for(i=0; i < strlen(filename); i++)
     if( *(filename+i) == ' ') *(filename+i) = '\0' ;
   if((fd=open(filename, O_RDWR|O_CREAT|O_TRUNC, PERMS)) < 0){
       if( (fd = open(filename, O_RDONLY, PERMS)) > 0)
           info_handler(0,NULL,"initdk: open file %s as READ ONLY",filename);
   }
   if( fd < 0 ) fd = creat(filename,PERMS);
   if(fd == -1)
     info_handler(ierr_2,filename,"initdk: cannot create file %s",filename);
   return(fd);
}

/* To write data into a previous opened file. This routine
   will wait until the write operations are completed.
  
   Calling sequence (from FORTRAN):
         nbytes = iowrit( chan, buff, bytes)
	 call iowrit(chan,buff,bytes)
   where:
         nbytes is the number bytes that transfered.
   
         chan is the file descriptor.

         buff is the buffer or array containing the data you
         wish to write.

         bytes is the number of bytes you wish to write.
*/ 
int iowrit(chan, buff, bytes)
int *chan, *bytes;
char *buff;
{  
   int nbytes;
   nbytes = write(*chan, buff, *bytes);

   if(nbytes != *bytes) 
     info_handler(ierr_3,raw_file,"iowrit: cannot write %d bytes of data",*bytes);

   return(nbytes);
}

/* To read data from a previously opened file. This routine will
   wait until after its operations are completed.

   Calling sequence (from FORTRAN):
       nbytes = ioread( chan, buff, bytes)
       call ioread( chan, buff, bytes)
   where:
       nbytes is the number bytes that transfered.
  
       chan is the file descriptor.
 
       buff is the buffer or array containning the data you wish
       to read.

       bytes is the number of bytes you wish to read.

 */
int ioread(chan, buff, bytes)
int *chan, *bytes ;
char *buff;
{  
   int nbytes;
   nbytes = read(*chan, buff, *bytes);

   if(nbytes != *bytes)
     info_handler(ierr_4,raw_file,"ioread: cannot read %d bytes of data",*bytes);
 
   return(nbytes);
}


/* To position the file pointer. This routine will call the lseek 
   to update the file pointer.

   Calling sequence (from FORTRAN):
      file_loc = ioseek(chan,loc_byte)
      call ioseek(chan,loc_byte)
   where:
        file_loc is the returned file location.

        chan is the file descriptor.

        loc_byte is byte location that requested to be set. This value
        must be greater or equal to zero for positioning the file at
        that location. If loc_byte is negative, the file pointer will
        move abs(loc_byte) from the current location.
        
*/

int ioseek(chan, loc_byte)
int *chan, *loc_byte;
{  
   int ibytes,nloc;
   ibytes = *loc_byte ;
   if(ibytes >= 0) nloc = lseek(*chan, ibytes, 0);
   else {
      ibytes = - ibytes;
      nloc = lseek(*chan, ibytes, 1);
   }
   return(nloc);
}



/* To close the file previously opened by initdk.

   Calling sequence (from FORTRAN):
      istatus = closedk( lun, chan)
      call closedk( lun, chan)
   where:
      istatus is the return value (0 is success, -1 is error)
 
      lun is the dummy variable to be compatible the VAX VMS call.

      chan is the file descriptor that you want to close.
 */

int closedk(lun,chan)
int *lun, *chan;
{
   return(close(*chan));
}



