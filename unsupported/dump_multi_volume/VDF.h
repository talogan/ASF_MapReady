/*
VDF.h:
	This header file contains structures for the VDR-- Volume
Descriptor Record, FPR-- File Pointer Records, and TDR-- Text Description
Records, which together make up a Volume Descriptor File (VDF).
*/
#include	<stdio.h>

typedef enum {
	volumeDescriptorRecord=192,
	filePointer=219,
	textRecord=18
} subType1Code;

typedef struct volume_descriptor_record {
	char    record_num[4];	            /* I*4 Record number */
	unsigned char	subType1,           /* First record sub-type code*/
		subType2, 	            /* 2nd record sub-type code (192)*/
		subType3,   	            /* 3rd record sub-type code (18)*/
		subType4;                   /* 4th record sub-type code (18)*/
	char 	record_length[4];	    /* Bytes in record, incl. header.*/
	char	ascii_ebcdic_flag[2],       /* ASCII or EBCDIC flag.*/
		blank[2],
		format_control_doc[12],     /* format control doc */
		format_control_vers[2],     /* version number of above. */
		format_control_revision[2], /* revision number of above. */
		facility_revision[12],      /* facitily release/revision code.*/
		tape_id[16],                /* Id of volume for this VDF */
		vol_id[16],                 /* Volume ID */
		set_id[16],                 /* Volume set ID */
		total_volumes[2],           /* I*2 # volumes in order */
		first_volume[2],            /* I*2 first tape number */
		last_volume[2],	            /* I*2 last tape number */
		this_volume[2],	            /* I*2 this volume number */
		crap_middle[60],	    /* 101 padding */
		numFPR[4],	            /* 161 # file ptr records in VDF */
		numRecords[4],	            /* 165 # records in VDF */
		crap_end[192];	            /* 169 irrelevent */
} volume_descriptor_record;

/* See comment below-- Divide by Zero when volume_descriptor_record
   is the wrong size.*/
const int vdrsizeFlag=1/(360==sizeof(volume_descriptor_record));


typedef struct file_pointer_record {
	char    record_num[4];	      /* I*4 Record number */
	unsigned char	subType1,     /* First record sub-type code*/
		subType2,             /* Second record sub-type code (192)*/
		subType3,             /* Third record sub-type code (18)*/
		subType4;             /* Fourth record sub-type code (18)*/
	char 	record_length[4];     /* # bytes in record, incl. header.*/
	char	ascii_ebcdic_flag[2], /* ASCII or EBCDIC flag.*/
		blank[2],
		ser_num[4],	      /* I*4 sequential number of the file */
		i_id[16],	      /* image ID string */
		f_type[28],	      /* File type string */
		f_code[4],	      /* File code: "SARL", "IMOP", or "SART" */
		d_type[28],	      /* Data type string */
		d_code[4],	      /* Data code:"MBAA","BINO","COMP","REAL"*/
		n_line[8],	      /* I*8 number of records */
		line_size[8],	      /* I*8 size of first record */
		max_rec_len[8],	      /* I*8 max record length */
		record_length_type[12],  /* rec_len type string*/
		record_length_code[4],   /* rec_len code: "FIXD" or "VARE" */
		phys_vol_start[2],    /* I*2 physical volume for file start */
		phys_vol_end[2],      /* I*2 physical volume for file end   */
		rec_phys_start[8],    /* I*8 First Record Number on Volume */
		rec_phys_end[8],      /* I*8 Last Record Number on Volume */
                crap_end[200];        /* No use defined */
} file_pointer_record;

/*This next lines are to make sure that we find out at compile time
  if the size of the above records are suddently no longer 360 bytes.
  If 360!=sizeof(fpr), then we'll get a divide-by-zero error at compile
  time, and we should change the above definition so we get exactly
  360 bytes.*/
const int fprsizeFlag=1/(360==sizeof(file_pointer_record));
