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
void match_pulse(double *fbuf,double *ovpulse,double *fbuf2);

#define SAMPLES_PER_LINE   13680	/* decoded samples per output line */
#define EXPANSION_FACTOR   4
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
  double total[O_SAMP_PER_LINE];
  
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

  int ocnt, this_line;
  double ocal[20];
  char intmp[256];

  if (argc != 3) {
     printf("Usage: %s <in> <out>\n",argv[0]);
     printf("\tin         - input data and header file base name\n");
     printf("\tout        - output data and header file base name\n\n");
     exit(1);
  }

  strcpy(indat,argv[1]); strcat(indat,".dat");
  strcpy(inhdr,argv[1]); strcat(inhdr,".hdr");
  strcpy(outdat,argv[2]); strcat(outdat,".dat");
  strcpy(outhdr,argv[2]); strcat(outhdr,".hdr");

  hdr = (SEASAT_header_ext *) malloc(sizeof(SEASAT_header_ext));

  printf("\t%s: opening input files...\n",argv[0]);
  fpin_dat = fopen(indat,"rb");
  if (fpin_dat == NULL) {printf("ERROR: Unable to open input data file %s\n",indat); exit(1);}  
  fpin_hdr = fopen(inhdr,"r");
  if (fpin_hdr == NULL) {printf("ERROR: Unable to open input header file %s\n",inhdr); exit(1);}  

  b = (complexFloat *) fftwf_malloc(sizeof(fftwf_complex)*O_SAMP_PER_FFT);
  b_longb = fftwf_plan_dft_1d(O_SAMP_PER_FFT, (fftwf_complex *) b, (fftwf_complex *)b, FFTW_BACKWARD, FFTW_MEASURE);
  b_longf = fftwf_plan_dft_1d(O_SAMP_PER_FFT, (fftwf_complex *) b, (fftwf_complex *)b, FFTW_FORWARD, FFTW_MEASURE);

  for (i=0; i<O_SAMP_PER_LINE; i++) { tbuf[i] = 0.0; total[i] = 0;}

  /* Start reading the input file header - get the first DWP */
  val=get_values(fpin_hdr,hdr);
  if (val!=20) {printf("ERROR: unable to read from header file\n"); exit(1);}
  last_dwp = hdr->delay;
  printf("\tfound initial dwp of %i\n",last_dwp);

  this_dwp = last_dwp;
  curr_line = 0;
  ndwps = 0;
  dwp_positions[ndwps++] = 0;
  printf("Set dwp position %i to %i\n",ndwps-1,dwp_positions[ndwps-1]);
  
  last_start = 1;

  printf("\treading file; creating oversampled average\n");
  double sum = 0.0;

  /* Loop through the entire file creating averaged calibration pulses */
  while (val == 20) {
      this_line = 0;
      /* read, oversample, and sum lines until DWP shift is found in the hdr file */
      while (this_dwp == last_dwp) {
        fread(buf,SAMPLES_PER_LINE,1,fpin_dat);
	if (this_line < 10000) {
          for (i=0; i<SAMPLES_PER_LINE; i++) fbuf[i] = (double) buf[i];
	  fft_oversamp(fbuf,SAMPLES_PER_LINE,EXPANSION_FACTOR,op);
          for (i=0; i<O_SAMP_PER_LINE; i++) /* if (op[i]>=4.0 && op[i]<28.0) */ { tbuf[i] += op[i]; total[i]++; }
	  this_line++;
   	  if ((int)this_line%100==0) printf("\toversampling line %i\n",curr_line);
	}
	
        curr_line++;
        val=get_values(fpin_hdr,hdr);
        if (val!=20) break;
        this_dwp = hdr->delay;
	if (curr_line%1000==1) printf("\tskipping line %i\n",curr_line);
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

      /* average the summed up pulse */
      sum = 0.0; 
      for (i=0; i<O_SAMP_PER_LINE; i++) { 
          if (total[i] != 0) tbuf[i] = tbuf[i]/total[i]; 
	  else tbuf[i] = 0.0;
          sum += tbuf[i];
        }
      sum = sum / (double) O_SAMP_PER_LINE;

      /* write the averaged normalized pulse to output file */
      sprintf(tmp,"%s%.7i.ov.pulse",argv[1],last_start); 
      fptmp = fopen(tmp,"w");
      if (fptmp == NULL) {printf("ERROR: Unable to open output data file %s\n",tmp); exit(1);}  
      for (i=0; i<O_SAMP_PER_LINE; i++) { 
          fprintf(fptmp,"%lf\n",(tbuf[i]-sum));
	  ovpulse[ndwps-1][i] = tbuf[i]-sum;  /* save averaged pulses as we go */
      }
      fclose(fptmp);

      /* calculate the spectra - just for fun (?) */
      spectra(fpin_dat,dwp_positions[(ndwps-1)],10002,ave,&ocnt,ocal);
      sprintf(intmp,"%s%.7i.spectra.out",argv[1],dwp_positions[(ndwps-1)]); 
      rename("spectra.out",intmp);
      sprintf(intmp,"%s%.7i.spectra.fixed",argv[1],dwp_positions[(ndwps-1)]); 
      rename("spectra.fixed",intmp);

      /* Since the spectra ALWAYS shows a peak at 0.25, we remove that here */
      for (k=0; k<O_SAMP_PER_LINE; k++) { b[k].real = tbuf[k]; b[k].imag = 0.0; }
      for (k=O_SAMP_PER_LINE; k<O_SAMP_PER_FFT; k++) { b[k].real = 0.0; b[k].imag = 0.0; }
      fftwf_execute(b_longf);

      fptmp = fopen("b.freq","w");
      for (k=0; k<O_SAMP_PER_LINE; k++) {
         tbuf[k] = sqrt(b[k].real*b[k].real+b[k].imag*b[k].imag)/(double)O_SAMP_PER_LINE;
         fprintf(fptmp,"%i %lf\n",k,tbuf[k]);
      }
      fclose(fptmp);

//      printf("Removing freq at location %i (value %lf)\n",O_SAMP_PER_FFT/4,
//                sqrt(b[O_SAMP_PER_FFT/4].real*b[O_SAMP_PER_FFT/4].real+
//		     b[O_SAMP_PER_FFT/4].imag*b[O_SAMP_PER_FFT/4].imag)/(double)O_SAMP_PER_FFT);
//      b[O_SAMP_PER_FFT/4].real = 0.0;
//      b[O_SAMP_PER_FFT/4].imag = 0.0;

//      printf("Removing freq at location %i (value %lf)\n",3*O_SAMP_PER_FFT/4,
//     		sqrt(b[3*O_SAMP_PER_FFT/4].real*b[3*O_SAMP_PER_FFT/4].real+
//		     b[3*O_SAMP_PER_FFT/4].imag*b[3*O_SAMP_PER_FFT/4].imag)/(double)O_SAMP_PER_FFT);
//      b[3*O_SAMP_PER_FFT/4].real = 0.0;
//      b[3*O_SAMP_PER_FFT/4].imag = 0.0;
  
      printf("Removing freq at location 4095 - 4097\n");
 
//      b[4093].real = 0.0;
//      b[4093].imag = 0.0;
//      b[4094].real = 0.0;
//      b[4094].imag = 0.0;
      b[4095].real = 0.0;
      b[4095].imag = 0.0;
      b[4096].real = 0.0;
      b[4096].imag = 0.0;
      b[4097].real = 0.0;
      b[4097].imag = 0.0;
//      b[4098].real = 0.0;
//      b[4098].imag = 0.0;
//      b[4099].real = 0.0;
//      b[4099].imag = 0.0;

      printf("Removing freq at location 1365\n");

      b[1365].real = 0.0;
      b[1365].imag = 0.0;

      fptmp = fopen("b.freq2","w");
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
      fptmp = fopen("b.time2","w");
      for (k=0; k<O_SAMP_PER_LINE; k++) {
         tbuf[k] -= sum;
         fprintf(fptmp,"%i %lf\n",k,tbuf[k]);
         ovpulse[ndwps-1][k] = tbuf[k];  /* save averaged pulses as we go */
      }
      fclose(fptmp);

      last_dwp = this_dwp;
      last_start = curr_line;
  }
  
  printf("\n");
  fclose(fpin_dat);
  fclose(fpin_hdr);
  fpin_dat = fopen(indat,"rb");
  if (fpin_dat == NULL) {printf("ERROR: Unable to open input data file %s\n",indat); exit(1);}  
  
  printf("Opening output file\n");
  fpout_dat = fopen(outdat,"wb");
  if (fpout_dat == NULL) {printf("ERROR: Unable to open output data file %s\n",outdat); exit(1);}  

  printf("Number of DWPs is %i\n",ndwps);

  /* Now that we have the pulses, apply them to the data file */
  sum = 0.0; curr_line =0;
  for (k=0; k<ndwps; k++) {

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

	
