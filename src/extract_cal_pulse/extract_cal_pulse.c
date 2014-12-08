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
#include <fftw3.h>

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
#define EXPANSION_FACTOR   4
#define O_SAMP_PER_LINE	   (SAMPLES_PER_LINE*EXPANSION_FACTOR)
#define FFT_LEN            16384
#define O_SAMP_PER_FFT     (FFT_LEN*EXPANSION_FACTOR)
#define MAX_DWPS	   1000
#define SAVE_PULSES	   1		/* set to zero to turn off */

int get_values(FILE *fp,SEASAT_header_ext *s);
void spectra(FILE *fp,int sl, int nl,double iqmean,int *ocnt,double *ocal);
double *fft_oversamp(double *in, int len, int exp_factor, double *out);
void match_pulse(double *fbuf,double *ovpulse,double *fbuf2);
int fix_dwps(int dwp_positions[],double pulse[][O_SAMP_PER_LINE],int ndwps);
void baseband(double *tbuf, double *total, double ovpulse[][O_SAMP_PER_LINE], int ndwps);

int PULSE_LENGTH=4000;          /* number of lines to sum to create the pulses */
complexFloat *b;
fftwf_plan   b_longb, b_longf;
double tbuf[O_SAMP_PER_LINE], fbuf[SAMPLES_PER_LINE], fbuf2[SAMPLES_PER_LINE];
double total[O_SAMP_PER_LINE];

main (int argc, char *argv[])
{
  int i, j, k, val;
  int curr_line;
  int last_dwp, this_dwp;
  int last_start;
  SEASAT_header_ext *hdr;
  FILE *fpin_dat, *fpin_hdr;
  FILE *fpout_dat, *fpout_hdr;
  FILE *fptmp;
  unsigned char buf[SAMPLES_PER_LINE];
  
  int    dwp_positions[MAX_DWPS];
  double ave;
  int    ndwps;
  double ovpulse[MAX_DWPS][O_SAMP_PER_LINE];

  double op[O_SAMP_PER_LINE];
  
  char indat[256], inhdr[256];
  char outdat[256], outhdr[256], outdis[256];
  char tmp[256];

  int ocnt;
  double ocal[20];
  char intmp[256];

  if (argc != 3 && argc != 4) {
     printf("Usage: %s <in> <out>\n",argv[0]);
     printf("\tin         - input data and header file base name\n");
     printf("\tout        - output data and header file base name\n\n");
     exit(1);
  }

  strcpy(indat,argv[1]); strcat(indat,".dat");
  strcpy(inhdr,argv[1]); strcat(inhdr,".hdr");
  strcpy(outdat,argv[2]); strcat(outdat,".dat");
  strcpy(outhdr,argv[2]); strcat(outhdr,".hdr");

  if (argc == 4) PULSE_LENGTH=atoi(argv[3]);
  printf("Using pulses of length %i\n",PULSE_LENGTH);

  hdr = (SEASAT_header_ext *) malloc(sizeof(SEASAT_header_ext));

  printf("%s: opening input files...\n",argv[0]);
  fpin_dat = fopen(indat,"rb");
  if (fpin_dat == NULL) {printf("ERROR: Unable to open input data file %s\n",indat); exit(1);}  
  fpin_hdr = fopen(inhdr,"r");
  if (fpin_hdr == NULL) {printf("ERROR: Unable to open input header file %s\n",inhdr); exit(1);}  
  printf("%s: opening output file\n",argv[0]);
  fpout_dat = fopen(outdat,"wb");
  if (fpout_dat == NULL) {printf("ERROR: Unable to open output data file %s\n",outdat); exit(1);}  

  b = (complexFloat *) fftwf_malloc(sizeof(fftwf_complex)*O_SAMP_PER_FFT);
  b_longb = fftwf_plan_dft_1d(O_SAMP_PER_FFT, (fftwf_complex *) b, (fftwf_complex *)b, FFTW_BACKWARD, FFTW_MEASURE);
  b_longf = fftwf_plan_dft_1d(O_SAMP_PER_FFT, (fftwf_complex *) b, (fftwf_complex *)b, FFTW_FORWARD, FFTW_MEASURE);

  for (i=0; i<O_SAMP_PER_LINE; i++) { tbuf[i] = 0.0; total[i] = 0;}

  /* Start reading the input file header - get the first DWP */
  val=get_values(fpin_hdr,hdr);
  if (val!=20) { printf("ERROR: unable to read from header file\n"); exit(1);}
  last_dwp = hdr->delay;
  printf("\tfound initial dwp of %i\n",last_dwp);

  this_dwp = last_dwp;
  curr_line = 0;
  ndwps = 0;
  dwp_positions[ndwps] = 0;
  printf("Set dwp position %i to %i\n",ndwps,dwp_positions[ndwps]);
  ndwps++;
  
  last_start = 1;
  printf("\treading file; creating oversampled average\n");
  double sum = 0.0;

  /* Loop through the entire file creating averaged calibration pulses */
  while (val == 20) {
  
      /* read, oversample, and sum lines until DWP shift is found in the hdr file */
      while (this_dwp == last_dwp) {
        fread(buf,SAMPLES_PER_LINE,1,fpin_dat);
        for (i=0; i<SAMPLES_PER_LINE; i++) fbuf[i] = (double) buf[i];
        fft_oversamp(fbuf,SAMPLES_PER_LINE,EXPANSION_FACTOR,op);
        for (i=0; i<O_SAMP_PER_LINE; i++) /* if (op[i]>=4.0 && op[i]<28.0) */ { tbuf[i] += op[i]; total[i]++; }
	
   	if ((int)curr_line%PULSE_LENGTH==0 && curr_line != 0) {
          printf("\tread to line %i.  DWP is %i. ",curr_line,this_dwp);
          dwp_positions[ndwps] = curr_line;
          printf("Set dwp position %i to %i\n",ndwps,dwp_positions[ndwps]);
          baseband(tbuf,total,ovpulse,ndwps);

          /* write the pulses to file as we go */
	  sprintf(tmp,"pulse%.6i.txt",dwp_positions[ndwps]);
	  fptmp = fopen(tmp,"w");
          for (i=0; i<O_SAMP_PER_LINE; i++) { 
	    fprintf(fptmp,"%lf\n",ovpulse[ndwps][i]);
  	  }
	  fclose(fptmp);
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
        printf("Set dwp position %i to %i\n",ndwps,dwp_positions[ndwps]);
      }

      baseband(tbuf,total,ovpulse,ndwps);
      ndwps++;
      last_dwp = this_dwp;
      last_start = curr_line;
  }

  printf("\n");
  fclose(fpin_hdr);
  
  printf("Calling fix_dwp with ndwps = %i\n",ndwps);
  ndwps = fix_dwps(dwp_positions,ovpulse,ndwps);
  printf("Returned from fix_dwps with ndwps = %i\n",ndwps);

  /* calculate the spectra - just for fun (?) */
//  spectra(fpin_dat,dwp_positions[(ndwps-1)],10002,ave,&ocnt,ocal);
//  sprintf(intmp,"%s%.7i.spectra.out",argv[1],dwp_positions[(ndwps-1)]); 
//  rename("spectra.out",intmp);
//  sprintf(intmp,"%s%.7i.spectra.fixed",argv[1],dwp_positions[(ndwps-1)]); 
//  rename("spectra.fixed",intmp);

  fseek(fpin_dat,0,SEEK_SET);

  /* Now that we have the pulses, apply them to the data file */
  sum = 0.0; curr_line =0;
  for (k=0; k<ndwps-1; k++) {
    printf("%i : trying range of %i to %i\n",k,dwp_positions[k],dwp_positions[k+1]);
    for (j=dwp_positions[k]; j<dwp_positions[k+1]; j++) {
      fread(buf,SAMPLES_PER_LINE,1,fpin_dat);
      sum =0.0; 
      for (i=0; i<SAMPLES_PER_LINE; i++) { fbuf[i] = (double) buf[i]; sum += fbuf[i]; }
      sum = sum / SAMPLES_PER_LINE;
      for (i=0; i<SAMPLES_PER_LINE; i++) fbuf[i] -= sum;
      match_pulse(fbuf,ovpulse[k],fbuf2);
      for (i=0; i<SAMPLES_PER_LINE; i++) { buf[i] = (int) (fbuf2[i]+sum); }
      fwrite(buf,SAMPLES_PER_LINE,1,fpout_dat); 
      curr_line+=1.0;
    }
  }

  sprintf(tmp,"cp %s %s\n",inhdr,outhdr);
  system(tmp);
  printf("Done correcting file, wrote %i lines of output\n\n",curr_line);
  fclose(fpout_dat);
  fclose(fpin_dat);
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

void baseband(double *tbuf, double *total, double ovpulse[][O_SAMP_PER_LINE], int ndwps) 
{
    int i;
    double sum;
    
    /* average the summed up pulse */
    sum = 0.0; 
    for (i=0; i<O_SAMP_PER_LINE; i++) { 
      if (total[i] != 0) tbuf[i] = tbuf[i]/total[i]; 
      else tbuf[i] = 0.0;
      sum += tbuf[i];
    }
    sum = sum / (double) O_SAMP_PER_LINE;

    /* save the basebanded pulse */
    for (i=0; i<O_SAMP_PER_LINE; i++) { 
      ovpulse[ndwps][i] = tbuf[i]-sum;
      tbuf[i] = 0.0;
      total[i] = 0.0;
    }
}
	
/* need to deal with the actual DWP changes AND remove the high frequency sinusoid */
int fix_dwps(int dwp_positions[],double pulse[][O_SAMP_PER_LINE],int ndwps)
{
  int i,j,k;
  
  int ldwps[MAX_DWPS];
  double lpulse[MAX_DWPS][O_SAMP_PER_LINE];
  int cnt;
  int first_time;
  FILE *fptmp;
  double sum = 0.0;
  char tmp[256];
  
  /* zero out the local values */
  cnt = 0;
  for (i=0; i<MAX_DWPS; i++) {
    ldwps[i] = 0;
    for (j=0; j< O_SAMP_PER_LINE; j++) lpulse[i][j] = 0;
  }
  
  /* copy in the first pulse */
  ldwps[cnt] = dwp_positions[0];
  for (j=0; j<O_SAMP_PER_LINE; j++) lpulse[cnt][j] = pulse[1][j];
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
        for (j=0; j<O_SAMP_PER_LINE; j++) lpulse[cnt][j] = pulse[i-1][j];
	first_time=0;
      } else {
        for (j=0; j<O_SAMP_PER_LINE; j++) lpulse[cnt][j] = pulse[i+1][j];
        first_time=1;
      }
      cnt++;
    } else {
      ldwps[cnt] = dwp_positions[i];
      for (j=0; j<O_SAMP_PER_LINE; j++) lpulse[cnt][j] = pulse[i][j];
      cnt++;
    }
  }
  
  /* copy into return variables */  
  for (i=0; i<cnt; i++) {
    dwp_positions[i] = ldwps[i];
    for (j=0; j<O_SAMP_PER_LINE; j++) pulse[i][j]=lpulse[i][j];
  }
  for (; i<MAX_DWPS; i++) {
    dwp_positions[i] = 0;
    for (j=0; j<O_SAMP_PER_LINE; j++) pulse[i][j]=0;
  }
  
  for (i=0; i<cnt; i++) {

      /* Now, loop through all pulses and remove the high frequency sinusoid */
      for (k=0; k<O_SAMP_PER_LINE; k++) { b[k].real = pulse[i][k]+16.0; b[k].imag = 0.0; }
      for (k=O_SAMP_PER_LINE; k<O_SAMP_PER_FFT; k++) { b[k].real = 0.0; b[k].imag = 0.0; }
      
      fftwf_execute(b_longf);

      sprintf(tmp,"pulse%.6i.freq1",dwp_positions[i]);
      fptmp = fopen(tmp,"w");
      printf("%i:writing freq1 file %s\n",dwp_positions[i],tmp);
      for (k=0; k<O_SAMP_PER_LINE; k++) {
         tbuf[k] = sqrt(b[k].real*b[k].real+b[k].imag*b[k].imag)/(double)O_SAMP_PER_LINE;
         fprintf(fptmp,"%i %lf\n",k,tbuf[k]);
      }
      fclose(fptmp);

      double mean = 0.0, diff=0.0, sqdiff=0.0, sumsq = 0.0, stddev;
      double DEVS = 2.0;
      int notch_count = 0;

      /* get the mean and standard deviation of the spectra */
      for (k=0; k<FFT_LEN; k++) mean = mean + sqrt(b[k].real*b[k].real+b[k].imag*b[k].imag); 
      mean = mean / FFT_LEN;
      for (k=0; k<FFT_LEN; k++) {
        diff = sqrt(b[k].real*b[k].real+b[k].imag*b[k].imag)-mean;	        /* deviation from mean */
        sqdiff = diff*diff;	/* square of deviation from mean */
        sumsq = sumsq + sqdiff;	/* sum of squared deviations */
      }
      stddev = sqrt(sumsq/(FFT_LEN-1));

      printf("Mean is %lf\n",mean);
      printf("Standard Deviation is %lf\n",stddev);

      /* find out where all of the notch filters need to be */
      int notch_cnt = 0;
//      for (k=0; k<FFT_LEN; k++) {
//        if ((sqrt(b[k].real*b[k].real+b[k].imag*b[k].imag)-mean)> DEVS*stddev) {
//           b[k].real = mean;
//           b[k].imag = mean;
//           notch_cnt++;
//        }
//      }      
      printf("Notched out %i values\n\n",notch_cnt);
      
      sprintf(tmp,"pulse%.6i.freq2",dwp_positions[i]);
      fptmp = fopen(tmp,"w");
      printf("%i:writing freq2 file %s\n",dwp_positions[i],tmp);
      for (k=0; k<O_SAMP_PER_LINE; k++) {
         tbuf[k] = sqrt(b[k].real*b[k].real+b[k].imag*b[k].imag)/(double)O_SAMP_PER_LINE;
         fprintf(fptmp,"%i %lf\n",k,tbuf[k]);
      }
      fclose(fptmp);
    
      fftwf_execute(b_longb);
    
      /* average the summed up pulse */
      sum = 0.0;
      for (k=0; k<O_SAMP_PER_LINE; k++) 
        { 
          tbuf[k] = sqrt(b[k].real*b[k].real+b[k].imag*b[k].imag)/(double)O_SAMP_PER_LINE;
          sum += tbuf[k];
        }
      sum = sum / (double) O_SAMP_PER_LINE;
      
      sprintf(tmp,"pulse%.6i.time2",dwp_positions[i]);
      fptmp = fopen(tmp,"w");
      printf("%i:writing time2 file %s\n",dwp_positions[i],tmp);
      for (k=0; k<O_SAMP_PER_LINE; k++) {
         tbuf[k] -= sum;
         fprintf(fptmp,"%i %lf\n",k,tbuf[k]);
         pulse[i][k] = tbuf[k];  /* save averaged pulses as we go */
      }
      fclose(fptmp);

   }

   return(cnt);
}
