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

	
	loop through the entire file 
	    for each DWP
	    	read line
		for first 1000 lines
			oversample line
			sum oversampled lines
		read next DWP
	    calculate the average pulse (sum/total line)
	    calculate the mean for the pulse
    
	    {save the (pulse-mean) for examination: *.ov.pulse}
    
	    convert the average pulse to complex
	    FFT pulse into freq space
    
	    {save the magnitude of the complex freq: b.freq}
    
	    remove unwanted frequencies...  WHICH ARE WHAT???
    
	    {save the magnitude of the complex freq: b.freq2}
    
	    FFT pulse back in time space
	    calculate the mean for the pulse
    
	    {save the (pulse-mean) for examination: b.time2}
	    {save the results of time2 -ov.pulse}
    
	    reset total lines, last_dwp, last_start line
   
	
ALGORITHM REFERENCES:

BUGS:

******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

typedef struct {
        int       major_cnt;
	long int  major_sync_loc;
	int       lsd_year;
	int       station_code;
	long int  msec;
	int       day_of_year;
	int       clock_drift;  
	int       no_scan_indicator_bit;
	int       bits_per_sample;
	int       mfr_lock_bit;
	int       prf_rate_code;
	int       delay;
	int       scu_bit;
	int       sdf_bit;
	int       adc_bit;
	int       time_gate_bit;
	int       local_prf_bit;
	int       auto_prf_bit;
	int       prf_lock_bit;
	int       local_delay_bit;
}  SEASAT_header_ext;

typedef struct {
   float real;
   float imag;
} complexFloat;

#define SAMPLES_PER_LINE   13680	/* decoded samples per output line */
#define MAX_DWPS	   10000

int get_values(FILE *fp,SEASAT_header_ext *s);
void write_values(FILE *fp,SEASAT_header_ext *s);
void spectra(FILE *fp,int sl, int nl,double iqmean,int *ocnt,double *ocal);
void match_pulse(double *fbuf,double *ovpulse,double *fbuf2);
void filter(unsigned char buf[], unsigned char obuf[]);
void create_smoothed_dat(FILE *fpin_dat,char *filename);
int fix_dwps(int dwp_positions[],double pulse[][SAMPLES_PER_LINE],int ndwps);

int PULSE_LENGTH=400;          /* number of lines to sum to create the pulses */

main (int argc, char *argv[])
{

  char indat[256], inhdr[256];
  char outdat[256], outhdr[256], outdis[256];
  char tmp[256];

  int i, j, k, val;
  int curr_line;
  int last_dwp, this_dwp;
  int last_start;
  SEASAT_header_ext *hdr;
  FILE *fpin_dat, *fpin_hdr;
  FILE *fpout_dat;
  FILE *fptmp;
  unsigned char buf[SAMPLES_PER_LINE];
  unsigned char obuf[SAMPLES_PER_LINE];
  unsigned char ubuf[SAMPLES_PER_LINE]; 
  double tbuf[SAMPLES_PER_LINE], fbuf[SAMPLES_PER_LINE], fbuf2[SAMPLES_PER_LINE];
  double total[SAMPLES_PER_LINE];
  double tot=0.0;
  
  int    dwp_positions[MAX_DWPS];
  double ave;
  int    ndwps;
  double pulse[MAX_DWPS][SAMPLES_PER_LINE];
  int ocnt;
  double ocal[20];
  char intmp[256];
  double max[1000000];
  int loc[1000000];
  

  if (argc != 3 && argc != 4) {
     printf("Usage: %s <in> <out> [window]\n",argv[0]);
     printf("\tin         - input data and header file base name\n");
     printf("\tout        - output data and header file base name\n");
     printf("\twindow     - length of summation window for pulses (default 400)\n\n");
     exit(1);
  }

  strcpy(indat,argv[1]); strcat(indat,".dat");
  strcpy(inhdr,argv[1]); strcat(inhdr,".hdr");
  strcpy(outdat,argv[2]); strcat(outdat,".dat");
  strcpy(outhdr,argv[2]); strcat(outhdr,".hdr");

  if (argc == 4) PULSE_LENGTH=atoi(argv[3]);
  printf("Using pulses of length %i\n",PULSE_LENGTH);
  
  hdr = (SEASAT_header_ext *) malloc(sizeof(SEASAT_header_ext));

  printf("\t%s: opening input files...\n",argv[0]);
  fpin_dat = fopen(indat,"rb");
  if (fpin_dat == NULL) {printf("ERROR1: Unable to open input data file %s\n",indat); exit(1);}  
  fpin_hdr = fopen(inhdr,"r");
  if (fpin_hdr == NULL) {printf("ERROR2: Unable to open input header file %s\n",inhdr); exit(1);}  

  strcpy(intmp,argv[1]); strcat(intmp,".smooth");
  printf("Creating smoothed dat file called %s\n",intmp);
  create_smoothed_dat(fpin_dat,intmp);
  FILE *fpin_smooth = fopen(intmp,"r");
  if (fpin_smooth == NULL) {printf("ERROR3: Unable to open input file test.smooth\n"); exit(1);}  

  for (i=0; i<SAMPLES_PER_LINE; i++) { tbuf[i] = 0.0; total[i] = 0.0; }
  for (i=0; i<256000; i++) { max[i] = 0.0; loc[i] = 0; }
    
  /* Start reading the input file header - get the first DWP */
  val=get_values(fpin_hdr,hdr);
  if (val!=20) {printf("ERROR4: unable to read from header file\n"); exit(1);}
  last_dwp = hdr->delay;
  printf("\tfound initial dwp of %i\n",last_dwp);

  this_dwp = last_dwp;
  curr_line = 0;
  ndwps = 0;
  dwp_positions[ndwps] = 0;
  printf("Set dwp position %i to %i\n",ndwps,dwp_positions[ndwps]);
  ndwps++;
  
  last_start = 1;
  printf("\treading file; creating average pulses\n");

  /* Loop through the entire file creating averaged calibration pulses */
  while (val == 20) {

      /* read and sum lines until DWP shift is found in the hdr file */
      while (this_dwp == last_dwp) {
        fread(buf,SAMPLES_PER_LINE,1,fpin_dat);
        fread(ubuf,SAMPLES_PER_LINE,1,fpin_smooth);
	filter(ubuf,obuf);

        for (i=0; i<SAMPLES_PER_LINE; i++) {
	  tbuf[i] += (double) buf[i]; total[i]++;
	  if ((double)obuf[i]>max[curr_line]) { max[curr_line] = (double) obuf[i]; loc[curr_line] = i; }
	}
	
   	if ((int)curr_line%PULSE_LENGTH==0 && curr_line != 0) {
          printf("\tread to line %i.  DWP is %i",curr_line,this_dwp);
          dwp_positions[ndwps] = curr_line;
          printf("\tSet dwp position %i to %i\n",ndwps,dwp_positions[ndwps]);
	  
          /* average the summed up pulse */
          double sum = 0.0;
          for (i=0; i<SAMPLES_PER_LINE; i++) 
            { 
     	      if (total[i] != 0) tbuf[i] = tbuf[i]/total[i]; 
	      else tbuf[i] = 0.0;
              sum += tbuf[i];
            }
          sum = sum / (double) SAMPLES_PER_LINE;
	  
	  /* save averaged pulses as we go */
          for (i=0; i<SAMPLES_PER_LINE; i++) { 
            pulse[ndwps][i] = tbuf[i]-sum;  
	    tbuf[i] = 0.0;
	    total[i] = 0.0;
          }

  	  ndwps++;
        }
	
        curr_line++;
        val=get_values(fpin_hdr,hdr);
        if (val!=20) break;
        this_dwp = hdr->delay;
      }

      /* check if we found a new DWP or hit end of file */
      if (this_dwp != last_dwp) {
        printf("\tread to line %i.  New DWP is %i\n",curr_line,this_dwp);
        dwp_positions[ndwps] = curr_line;
        printf("Set dwp position %i to %i\n",ndwps,dwp_positions[ndwps]);
      } else if (val!=20) {
        printf("\tread to end of file\n");
	dwp_positions[ndwps] = curr_line;
        printf("\t\tSet dwp position %i to %i\n",ndwps,dwp_positions[ndwps]);
      } 

      /* average the summed up pulse */
      double sum = 0.0;
      for (i=0; i<SAMPLES_PER_LINE; i++) 
        { 
	  if (total[i] != 0) tbuf[i] = tbuf[i]/total[i]; 
	  else tbuf[i] = 0.0;
          sum += tbuf[i];
        }
      sum = sum / (double) SAMPLES_PER_LINE;

      /* save averaged pulses as we go */
      for (i=0; i<SAMPLES_PER_LINE; i++) { 
         pulse[ndwps][i] = tbuf[i]-sum;  
	 tbuf[i] = 0.0;
	 total[i] = 0.0;
      }
      
      ndwps++;

      last_dwp = this_dwp;
      last_start = curr_line;
  }
  
  /* write the maximums to file */
  sprintf(tmp,"%s%.7i.max",argv[1],last_start); 
  fptmp = fopen(tmp,"w");
  if (fptmp == NULL) {printf("ERROR7: Unable to open output data file %s\n",tmp); exit(1);}  
  for (i=0; i<curr_line; i++) { 
     fprintf(fptmp,"%i %lf\n",loc[i],max[i]);
  }
  fclose(fptmp);
  
  printf("Calling fix_dwps with ndwps = %i\n",ndwps);
  ndwps = fix_dwps(dwp_positions,pulse,ndwps);
  printf("Returned from fix_dwps with ndwps = %i\n",ndwps);
  
  for (k=0; k<ndwps; k++) {
    sprintf(tmp,"%s%.7i.pulse",argv[1],dwp_positions[k]); 
    fptmp = fopen(tmp,"w");
    if (fptmp == NULL) {printf("ERROR7: Unable to open output data file %s\n",tmp); exit(1);}  
    for (i=0; i<SAMPLES_PER_LINE; i++) { 
      fprintf(fptmp,"%lf\n",pulse[k][i]);
    }
    fclose(fptmp);
  }

  printf("\n");
  fclose(fpin_hdr);
  fseek(fpin_dat,0,SEEK_SET);

  printf("Opening output file\n");
  fpout_dat = fopen(outdat,"wb");
  if (fpout_dat == NULL) {printf("ERROR9: Unable to open output data file %s\n",outdat); exit(1);}  

  printf("Number of DWPs is %i\n",ndwps);

  /* Now that we have the pulses, apply them to the data file */
  for (k=0; k<ndwps-1; k++) {
    printf("%i : trying range of %i to %i\n",k,dwp_positions[k],dwp_positions[k+1]);
    for (j=dwp_positions[k]; j<dwp_positions[k+1]; j++) {
      fread(buf,SAMPLES_PER_LINE,1,fpin_dat);
      for (i=0; i<SAMPLES_PER_LINE; i++) {
         fbuf[i] = (double) buf[i];

// Restict code
//	 if (i>6380 && i<7940) fbuf2[i] = fbuf[i] - pulse[k+1][i];
//	 else                  fbuf2[i] = fbuf[i];
	 
	 fbuf2[i] = fbuf[i] - pulse[k+1][i];
	 buf[i] = (int) fbuf2[i]; 
      }
      fwrite(buf,SAMPLES_PER_LINE,1,fpout_dat); 
      tot+=1.0;
    }
  }

  sprintf(tmp,"cp %s %s\n",inhdr,outhdr);
  system(tmp);

  printf("Done correcting file, wrote %lf lines of output\n\n",tot);
  fclose(fpout_dat);

  exit(0);
}


int get_values(FILE *fp,SEASAT_header_ext *s) {
  int val;
  if (s==NULL) {printf("empty pointer passed to get_values\n"); exit(1);}
  val = fscanf(fp,"%i %li %i %i %i %li %i %i %i %i %i %i %i %i %i %i %i %i %i %i\n",
    &(s->major_cnt),&(s->major_sync_loc),&(s->station_code),&(s->lsd_year),
    &(s->day_of_year),&(s->msec),&(s->clock_drift),&(s->no_scan_indicator_bit),
    &(s->bits_per_sample),&(s->mfr_lock_bit),&(s->prf_rate_code),&(s->delay),
    &(s->scu_bit),&(s->sdf_bit),&(s->adc_bit),&(s->time_gate_bit),&(s->local_prf_bit),
    &(s->auto_prf_bit),&(s->prf_lock_bit),&(s->local_delay_bit));
  return(val);
}

void write_values(FILE *fp,SEASAT_header_ext *s) {
  if (s==NULL) {printf("empty pointer passed to write values\n"); exit(1);}
  if (fp==NULL) {printf("null file pointer passed to write values\n"); exit(1);}
  fprintf(fp,"%i %li %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i\n",
    s->major_cnt,s->major_sync_loc,s->station_code,s->lsd_year,
    s->day_of_year,s->msec,s->clock_drift,s->no_scan_indicator_bit,
    s->bits_per_sample,s->mfr_lock_bit,s->prf_rate_code,s->delay,
    s->scu_bit,s->sdf_bit,s->adc_bit,s->time_gate_bit,s->local_prf_bit,
    s->auto_prf_bit,s->prf_lock_bit,s->local_delay_bit);
}

#define KERNEL 121
#define KERNEL2 5

void filter(unsigned char buf[], unsigned char obuf[])
{
  int sum = 0;
  int tot = 0;
  int i, j;
  
  /* printf("In filter: kernel = %i  kernel/2 = %i  SAMPLES_PER_LINE-KERNEL/2 = %i\n",KERNEL,KERNEL/2,
  SAMPLES_PER_LINE-KERNEL/2); */
  
  /* copy in first unfiltered section */
  for (i=0; i<KERNEL2/2; i++) {
     obuf[i] = buf[i];
  }
  
  /* run the filter on the data */
  for (i=KERNEL2/2;i<SAMPLES_PER_LINE-KERNEL2/2;i++)
    {
      sum = 0;
      for (j=-KERNEL2/2; j<KERNEL2/2+1; j++) sum += (int) buf[i+j];
      obuf[i] = (unsigned char) ((double)sum / (double) KERNEL2);
    }
  
  /* copy in the last unfiltered section */
  for (i=SAMPLES_PER_LINE-KERNEL2/2 ;i<SAMPLES_PER_LINE; i++) obuf[i] = (unsigned char) buf[i];
}
	

void create_smoothed_dat(FILE *fpin_dat,char *filename)
{
  unsigned char buf[SAMPLES_PER_LINE];
  int **tbuf;
  int *tmp;
  int sum[SAMPLES_PER_LINE];
  unsigned char obuf[SAMPLES_PER_LINE];
  int curr_line;
  int i,j, line;
  int cnt;
  FILE *fpout;
  
  curr_line = 0;

  if( access( filename, F_OK ) != -1 ) {
    // file exists
    printf("Smoothed file already exists\n");

  } else {
    // file doesn't exist
    fpout = fopen(filename,"wb");

    /* write out first lines */
    for (i=0;i<KERNEL/2; i++) {
      cnt = fread(buf,1,SAMPLES_PER_LINE, fpin_dat);
      fwrite(buf,1,SAMPLES_PER_LINE, fpout);
      curr_line++;
    }
  
    while (cnt == SAMPLES_PER_LINE) {
      for (j=0;j<SAMPLES_PER_LINE; j++) sum[j] = 0;
      for (line = curr_line-KERNEL/2; line < curr_line+KERNEL/2+1; line++)
        {
          long int where = (long int) (line*SAMPLES_PER_LINE);
          fseek(fpin_dat,where,SEEK_SET);
  	  cnt = fread(buf,1,SAMPLES_PER_LINE, fpin_dat);
          for (j=0;j<SAMPLES_PER_LINE; j++) sum[j] += (int) buf[j];
        }
      for (j=0;j<SAMPLES_PER_LINE; j++) obuf[j] = (unsigned int) (sum[j]/KERNEL);
      fwrite(obuf,1,SAMPLES_PER_LINE, fpout);
      curr_line++;
    }
  
    /* dump out last lines */
    for (i=KERNEL/2+1; i<KERNEL; i++) {
      long int where = (long int) (curr_line*SAMPLES_PER_LINE);
      fseek(fpin_dat,where,SEEK_SET);
      cnt = fread(buf,1,SAMPLES_PER_LINE, fpin_dat);
      fwrite(buf,1,SAMPLES_PER_LINE, fpout);
      curr_line++;
    }
  
    fseek(fpin_dat,0,SEEK_SET);
    fclose(fpout);
  }

}

int fix_dwps(int dwp_positions[],double pulse[][SAMPLES_PER_LINE],int ndwps)
{
  int i,j;
  
  int ldwps[MAX_DWPS];
  double lpulse[MAX_DWPS][SAMPLES_PER_LINE];
  int cnt;
  int first_time;
  
  /* zero out the local values */
  cnt = 0;
  for (i=0; i<MAX_DWPS; i++) {
    ldwps[i] = 0;
    for (j=0; j< SAMPLES_PER_LINE; j++) lpulse[i][j] = 0;
  }
  
  /* copy in the first pulse */
  ldwps[cnt] = dwp_positions[0];
  for (j=0; j<SAMPLES_PER_LINE; j++) lpulse[cnt][j] = pulse[1][j];
  cnt++;
  first_time = 1;
  
  /* loop through all pulses */
  for (i=1; i< ndwps; i++)  {

    /* if this pulse covers less than a PULSE_LENGTH discard 
       first time, use the last one; 2nd time use the next one */
    if (dwp_positions[i]-dwp_positions[i-1]<PULSE_LENGTH) {

      printf("Fixing values at %i\n",dwp_positions[i]);
      
      ldwps[cnt] = dwp_positions[i];
   
      if (first_time == 1) {
        for (j=0; j<SAMPLES_PER_LINE; j++) lpulse[cnt][j] = pulse[i-1][j];
	first_time=0;
      } else {
        for (j=0; j<SAMPLES_PER_LINE; j++) lpulse[cnt][j] = pulse[i+1][j];
        first_time=1;
      }
      cnt++;
    } else {
      ldwps[cnt] = dwp_positions[i];
      for (j=0; j<SAMPLES_PER_LINE; j++) lpulse[cnt][j] = pulse[i][j];
      cnt++;
    }
  }
  
  /* copy into return variables */  
  for (i=0; i<cnt; i++) {
    dwp_positions[i] = ldwps[i];
    for (j=0; j<SAMPLES_PER_LINE; j++) pulse[i][j]=lpulse[i][j];
  }
  for (; i<MAX_DWPS; i++) {
    dwp_positions[i] = 0;
    for (j=0; j<SAMPLES_PER_LINE; j++) pulse[i][j]=0;
  }
    
  return(cnt);
}

