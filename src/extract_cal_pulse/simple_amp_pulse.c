/******************************************************************************
NAME: extract_cal_pulse

SYNOPSIS:  extract_cal_pulse <infile> <outfile>

DESCRIPTION:  Intent is to remove the calibration pulses present in Seasat data

EXTERNAL ASSOCIATES:
    NAME:               USAGE:
    ---------------------------------------------------------------

FILE REFERENCES:
    NAME:               USAGE:
    ---------------------------------------------------------------

PROGRAM HISTORY:
    VERS:   DATE:  AUTHOR:      PURPOSE:
    ---------------------------------------------------------------
     
HARDWARE/SOFTWARE LIMITATIONS:

ALGORITHM DESCRIPTION:
	
ALGORITHM REFERENCES:

BUGS:

******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define SAMPLES_PER_LINE   6840	/* decoded samples per output line */
#define O_SAMP_PER_LINE	   6840
#define NUM_LINES 	   6234

void byteswap(void *buf,int len);

main (int argc, char *argv[])
{
  int i, j, k, val;
  int curr_line;
  int last_dwp, this_dwp;
  int last_start;
  FILE *fpin_dat, *fpin_hdr;
  FILE *fpout_dat, *fpout_hdr;
  FILE *fptmp;
  float  buf[SAMPLES_PER_LINE];
  double tbuf[O_SAMP_PER_LINE], fbuf[SAMPLES_PER_LINE], fbuf2[SAMPLES_PER_LINE];
  double total[SAMPLES_PER_LINE];
  double tot=0.0;
  
  int    dwp_positions[256];
  double ave;
  int    ndwps;
  double pulse[10][O_SAMP_PER_LINE];

  double op[O_SAMP_PER_LINE];
  
  char indat[256], inhdr[256];
  char outdat[256], outhdr[256], outdis[256];
  char tmp[256];

  int ocnt;
  double ocal[20];
  char intmp[256];

  if (argc != 3) {
     printf("Usage: %s <in> <out>\n",argv[0]);
     printf("\tin         - input img and meta file base name\n");
     printf("\tout        - output img and meta file base name\n\n");
     exit(1);
  }

  strcpy(indat,argv[1]); strcat(indat,".img");
  strcpy(inhdr,argv[1]); strcat(inhdr,".meta");
  strcpy(outdat,argv[2]); strcat(outdat,".img");
  strcpy(outhdr,argv[2]); strcat(outhdr,".meta");

  printf("\t%s: opening input files...\n",argv[0]);
  fpin_dat = fopen(indat,"rb");
  if (fpin_dat == NULL) {printf("ERROR: Unable to open input data file %s\n",indat); exit(1);}  

  for (i=0; i<SAMPLES_PER_LINE; i++) tbuf[i] = 0.0; 
  for (i=0; i<SAMPLES_PER_LINE; i++) total[i] = 0.0; 

  printf("\treading file; creating average\n");

  /* Loop through the entire file creating averaged calibration pulse */
  for (j=0; j< NUM_LINES; j++) {
      fread(buf,SAMPLES_PER_LINE*sizeof(float),1,fpin_dat);
      byteswap(buf,SAMPLES_PER_LINE);
      for (i=0; i<SAMPLES_PER_LINE; i++) { 
	 tbuf[i] += (double) buf[i];
         // if (i<10) printf("buf = %f; tbuf[%i] = %lf ",buf[i],i,tbuf[i]);
      }  // printf("\n");
      if ((int)curr_line%1000==0) printf("\treading line %i\n",curr_line);
      curr_line++;
  }

  /* average the summed up pulse */
  double sum = 0.0;
  for (i=0; i<SAMPLES_PER_LINE; i++) { 
      tbuf[i] = tbuf[i]/(double) NUM_LINES; 
      sum += tbuf[i];
  }
  sum = sum / (double) SAMPLES_PER_LINE;
  for (i=0; i<SAMPLES_PER_LINE; i++) tbuf[i] = tbuf[i]-sum;
  
  /* write the averaged pulse to output file */
  sprintf(tmp,"%s.ave",argv[1],last_start); 
  fptmp = fopen(tmp,"w");
  if (fptmp == NULL) {printf("ERROR: Unable to open output data file %s\n",tmp); exit(1);}  
  for (i=0; i<SAMPLES_PER_LINE; i++) { 
    fprintf(fptmp,"%lf\n",tbuf[i]);
  }
  fclose(fptmp);

  printf("\n");
  fclose(fpin_dat);
  fpin_dat = fopen(indat,"rb");
  if (fpin_dat == NULL) {printf("ERROR: Unable to open input data file %s\n",indat); exit(1);}  
  
  printf("Opening output file\n");
  fpout_dat = fopen(outdat,"wb");
  if (fpout_dat == NULL) {printf("ERROR: Unable to open output data file %s\n",outdat); exit(1);}  

  /* Now that we have the pulse, apply it to the data file */
  for (j=0 ; j<NUM_LINES; j++) {
    fread(buf,SAMPLES_PER_LINE*sizeof(float),1,fpin_dat);
    byteswap(buf,SAMPLES_PER_LINE);
    for (i=0; i<SAMPLES_PER_LINE; i++) {
        fbuf[i] = (double) buf[i];
	if (tbuf[i]>=1.5) fbuf2[i] = fbuf[i] - tbuf[i]; 
	else fbuf2[i] = fbuf[i];
	buf[i] = (float) fbuf2[i]; 
    }
    byteswap(buf,SAMPLES_PER_LINE);
    fwrite(buf,SAMPLES_PER_LINE*sizeof(float),1,fpout_dat); 
    tot+=1.0;
  }


  sprintf(tmp,"cp %s %s\n",inhdr,outhdr);
  system(tmp);
  printf("Done correcting file, wrote %lf lines of output\n\n",tot);
  fclose(fpout_dat);

  exit(0);
}


void byteswap(void *buf,int len)
{
   unsigned char *bufptr = (unsigned char *) buf;
   unsigned char save1, save2;
   int i, cnt;

   for (i=0; i<len; i++)
    {
        cnt = i*sizeof(float);  
        save1 = bufptr[cnt];
        bufptr[cnt] = bufptr[cnt+3];
        bufptr[cnt+3] = save1;

        save2 = bufptr[cnt+1];
        bufptr[cnt+1] = bufptr[cnt+2];
        bufptr[cnt+2] = save2;
     }   
    return; 
}

