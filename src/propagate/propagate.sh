#!/bin/sh
#
# NAME: propagate - automates orbital propagation
#
# SYNOPSIS: propagate <inVectFile> <outVectFile>
#
# DESCRIPTION:
#      	
#        Uses gen_oe and ASAP to propagate output state vectors (outVectFile)
#        using the input state vectors and desired time (inVectFile).
#        The format of the inVectFile is that created by parse_lzp_par (see
#        that man page).  The output is simply a file containing lines of
#        time, state vector
#        for each propagation step requested.
#
#        This is a simple script that was designed to overcome all of the
#        hardcoded file names used by gen_oe and asap.
#
# EXTERNAL ASSOCIATES:
#    NAME:               USAGE:
#    ---------------------------------------------------------------
#
# FILE REFERENCES:
#    NAME:               USAGE:
#    ---------------------------------------------------------------
#
# PROGRAM HISTORY:
#    VERS:   DATE:  AUTHOR:      PURPOSE:
#    ---------------------------------------------------------------
#    1.0	    Tom Logan
#
# HARDWARE/SOFTWARE LIMITATIONS:
#
# ALGORITHM DESCRIPTION:
#
# ALGORITHM REFERENCES:
#
# BUGS
#	
#****************************************************************************
#								            *
#   propagate - automates orbital propagation				    *
#   Copyright (C) 2001  ASF Advanced Product Development    	    	    *
#									    *
#   This program is free software; you can redistribute it and/or modify    *
#   it under the terms of the GNU General Public License as published by    *
#   the Free Software Foundation; either version 2 of the License, or       *
#   (at your option) any later version.					    *
#									    *
#   This program is distributed in the hope that it will be useful,	    *
#   but WITHOUT ANY WARRANTY; without even the implied warranty of    	    *
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the   	    *
#   GNU General Public License for more details.  (See the file LICENSE     *
#   included in the asf_tools/ directory).				    *
#									    *
#   You should have received a copy of the GNU General Public License       *
#   along with this program; if not, write to the Free Software		    *
#   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.               *
#									    *
#       ASF Advanced Product Development LAB Contacts:			    *
#	APD E-mail:	apd@asf.alaska.edu 				    *
# 									    *
#	Alaska SAR Facility			APD Web Site:	            *	
#	Geophysical Institute			www.asf.alaska.edu/apd	    *
#      	University of Alaska Fairbanks					    *
#	P.O. Box 757320							    *
#	Fairbanks, AK 99775-7320					    *
#								  	    *
#**************************************************************************/
#
if [ ! $# -eq "2" ]
then
        echo ""
        echo "Usage: $0 <inVectfile> <outVectfile>"
        echo ""
        echo "Propagates state vectors using gen_oe and ASAP."
        echo "Version 1.00, ASF SAR TOOLS" 
        echo ""
        exit 1
fi


if [ -r "./propagate_lock_file" ]
then
	cat <<EOF
propagate_lock_file found in '`pwd`'.

It is possible that ths lock file is stale, and that no competing
process is trying to create or use new temporary files in 
'`pwd`'.  
If this is the case, delete the lock file and try again.
EOF
	>&2
  	exit 1
else
	touch ./propagate_lock_file
fi


#echo "Generating ASAP input file"
gen_oe $1

#echo " "
#echo "Running ASAP"
#echo " "
asap
/bin/mv fort.7 $2
/bin/rm ./temp_file_formerly_known_as_5

#echo " "
#echo "Propagate:  Created output file $2"
#echo " "

rm ./propagate_lock_file

 
