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

int get_values(FILE *fp,SEASAT_header_ext *s);
void write_values(FILE *fp,SEASAT_header_ext *s);
void spectra(FILE *fp,int sl, int nl,double iqmean,int *ocnt,double *ocal);
double *fft_oversamp(double *in, int len, int exp_factor, double *out);
void match_sinc(double *fbuf,double *fbuf2, double max, int best);

#define SAMPLES_PER_LINE   13680	/* decoded samples per output line */
#define EXPANSION_FACTOR   8
#define O_SAMP_PER_LINE	   (SAMPLES_PER_LINE*EXPANSION_FACTOR)
#define FFT_LEN            16384
#define O_SAMP_PER_FFT     (FFT_LEN*EXPANSION_FACTOR)

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
  double tbuf[O_SAMP_PER_LINE], fbuf[SAMPLES_PER_LINE], fbuf2[SAMPLES_PER_LINE];
  double total=0.0;
  
  int    dwp_positions[256];
  double ave;
  int    ndwps;
  double ovpulse[10][O_SAMP_PER_LINE];

  complexFloat *b;
  fftwf_plan   b_longb, b_longf;
  double op[O_SAMP_PER_LINE];
  
  char indat[256], inhdr[256];
  char outdat[256], outhdr[256], outdis[256];
  char tmp[256];

  int ocnt;
  int best;
  double ocal[20];
  char intmp[256];

  if (argc != 3 && argc != 4) {
     printf("Usage: %s <in> <out> [offset]\n",argv[0]);
     printf("\tin         - input data and header file base name\n");
     printf("\tout        - output data and header file base name\n");
     printf("\t[offset]   - offset of sinc in the data file\n\n");
     exit(1);
  }

  strcpy(indat,argv[1]); strcat(indat,".dat");
  strcpy(inhdr,argv[1]); strcat(inhdr,".hdr");
  strcpy(outdat,argv[2]); strcat(outdat,".dat");
  strcpy(outhdr,argv[2]); strcat(outhdr,".hdr");

  if (argc == 4) best = atoi(argv[3]);
  else best = 6384;

  hdr = (SEASAT_header_ext *) malloc(sizeof(SEASAT_header_ext));

  printf("\t%s: opening input files...\n",argv[0]);
  fpin_dat = fopen(indat,"rb");
  if (fpin_dat == NULL) {printf("ERROR: Unable to open input data file %s\n",indat); exit(1);}  
  fpin_hdr = fopen(inhdr,"r");
  if (fpin_hdr == NULL) {printf("ERROR: Unable to open input header file %s\n",inhdr); exit(1);}  

  b = (complexFloat *) fftwf_malloc(sizeof(fftwf_complex)*O_SAMP_PER_FFT);
  b_longb = fftwf_plan_dft_1d(O_SAMP_PER_FFT, (fftwf_complex *) b, (fftwf_complex *)b, FFTW_BACKWARD, FFTW_MEASURE);
  b_longf = fftwf_plan_dft_1d(O_SAMP_PER_FFT, (fftwf_complex *) b, (fftwf_complex *)b, FFTW_FORWARD, FFTW_MEASURE);

  for (i=0; i<O_SAMP_PER_LINE; i++) tbuf[i] = 0.0;

  /* Start reading the input file header - get the first DWP */
  val=get_values(fpin_hdr,hdr);
  if (val!=20) {printf("ERROR: unable to read from header file\n"); exit(1);}
  last_dwp = hdr->delay;
  printf("\tfound initial dwp of %i\n",last_dwp);

  this_dwp = last_dwp;
  curr_line = 0;
  ndwps = 0;
  dwp_positions[ndwps++] = 0;
  printf("\tSet dwp position %i to %i\n",ndwps-1,dwp_positions[ndwps-1]);
  
  last_start = 1;

  printf("\treading file; searching for DWP positions\n");
  double sum = 0.0;
  double max = 0.0;

  /* Loop through the entire file creating averaged calibration pulses */
  while (val == 20) {

      /* read, oversample, and sum lines until DWP shift is found in the hdr file */
      while (this_dwp == last_dwp) {
        fread(buf,SAMPLES_PER_LINE,1,fpin_dat);
        curr_line++;
        val=get_values(fpin_hdr,hdr);
        if (val!=20) break;
        this_dwp = hdr->delay;
	// if (curr_line%1000==1) printf("\tskipping line %i\n",curr_line);
      }

      /* check if we found a new DWP or hit end of file */
      if (this_dwp != last_dwp) {
        printf("\tread to line %i.  New DWP is %i\n",curr_line,this_dwp);
        dwp_positions[ndwps++] = curr_line;
        printf("Set dwp position %i to %i\n",ndwps-1,dwp_positions[ndwps-1]);
      } else if (val!=20) {
        printf("\tread to end of file\n");
	dwp_positions[ndwps] = curr_line;
        printf("Set dwp position %i to %i\n",ndwps,dwp_positions[ndwps]);
      }
  }
    
  printf("\n");
  fpin_dat = fopen(indat,"rb");
  if (fpin_dat == NULL) {printf("ERROR: Unable to open input data file %s\n",indat); exit(1);}  
  
  printf("Opening output file\n");
  fpout_dat = fopen(outdat,"wb");
  if (fpout_dat == NULL) {printf("ERROR: Unable to open output data file %s\n",outdat); exit(1);}  

  printf("Number of DWPs is %i\n",ndwps);

  /* Now that we have the DWPs, go through them */
  sum = 0.0;
  
  printf("\tUsing best fit of %i\n",best);
  for (k=0; k<ndwps; k++) {
    printf("%i : trying range of %i to %i\n",k,dwp_positions[k],dwp_positions[k+1]);

    for (j=dwp_positions[k]; j<dwp_positions[k+1]; j++) {
      fread(buf,SAMPLES_PER_LINE,1,fpin_dat);
      sum =0.0; max =0.0;
      for (i=0; i<SAMPLES_PER_LINE; i++) {
          fbuf[i] = (double) buf[i]; 
	  sum += fbuf[i]; 
      }
      sum = sum / SAMPLES_PER_LINE;
      for (i=0; i<SAMPLES_PER_LINE; i++) { fbuf[i] -= sum; if (fabs(fbuf[i])>max) max = fabs(fbuf[i]); }

      match_sinc(fbuf,fbuf2,max,best);

      for (i=0; i<SAMPLES_PER_LINE; i++) { 
        buf[i] = (int) (fbuf2[i]+sum); 
	if (buf[i]>32) { printf("ERROR: value at %i out of range %i (sum %lf)\n",i,buf[i],sum); buf[i] = 32; }
      }
      fwrite(buf,SAMPLES_PER_LINE,1,fpout_dat); 
      total+=1.0;
      if ((int)total%10000 == 1) printf("\tprocessing line %lf\n",total);
    }
  }

  sprintf(tmp,"cp %s %s\n",inhdr,outhdr);
  system(tmp);

  printf("Done correcting file, wrote %f lines of output\n\n",total);
  
  fclose(fpout_dat);
//  fclose(fpout_hdr);

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

	
