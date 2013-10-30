/******************************************************************************
NAME:

SYNOPSIS:
	

DESCRIPTION:

EXTERNAL ASSOCIATES:
    NAME:               USAGE:
    ---------------------------------------------------------------

FILE REFERENCES:
    NAME:               USAGE:
    ---------------------------------------------------------------

PROGRAM HISTORY:
    VERS:   DATE:    AUTHOR:      PURPOSE:
    ---------------------------------------------------------------
    1.0	    9/18/13  T. Logan     Seasat project
    
HARDWARE/SOFTWARE LIMITATIONS:

ALGORITHM DESCRIPTION:

ALGORITHM REFERENCES:

BUGS:

******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>

#define SAMPLES_PER_LINE   13680	/* decoded samples per output line */
#define FFT_LEN		   2048
#define EXPANSION_FACTOR   8
#define  NPTS  1700
#define  pi     3.14159

typedef struct {
   float real;
   float imag;
} complexFloat;


void match_sinc(double *fbuf,double *fbuf2,double max, int best)
{
    int i, j;
    double diff[SAMPLES_PER_LINE];
    double min_diff;
    double x, y;
    double fs = 22765000.0;
    double pulsedur = 3.37999991E-05;
    int    npts=2.*fs*pulsedur;
    double ts=1./(2.*fs);
    int    k=0;
    int    k1start=0;
    int    k1end=(int) npts;
    int    k2start=0;
    int    k2end=(int) npts;
    double t, phase;
    double slope1 = 5.62130190E+11;
    
    static int first = 1;
    static double ref[NPTS];
    
    if (first == 1) {
      printf("npts = %i\n",npts);
      for (i=-npts/2; i<npts/2; i++)
       {
	   k=k+1;
	   t=(double)i*ts;
	   phase = pi*slope1*t*t+pi*fs*t;
	   // printf("phase %lf\n",phase);
	   ref[i+npts/2+1]=cos(phase);
       }

      FILE *fp;
      char tmp[256];
      sprintf(tmp,"synth_chirp%4i.time",best);
      fp = fopen(tmp,"w");
      for (i=0; i<SAMPLES_PER_LINE; i++) {
	if (i<best) fprintf(fp,"0.0\n");
	else if (i<best+npts) fprintf(fp,"%lf\n",ref[i-best]);
	else fprintf(fp,"0.0\n");
      }
      fclose(fp);

      complexFloat *b;
      fftwf_plan   b_longb, b_longf;
      b = (complexFloat *) fftwf_malloc(sizeof(fftwf_complex)*FFT_LEN);
      b_longb = fftwf_plan_dft_1d(FFT_LEN, (fftwf_complex *) b, (fftwf_complex *)b, FFTW_BACKWARD, FFTW_MEASURE);
      b_longf = fftwf_plan_dft_1d(FFT_LEN, (fftwf_complex *) b, (fftwf_complex *)b, FFTW_FORWARD, FFTW_MEASURE);
      
      for (k=0; k<npts; k++) { b[k].real = ref[k]; b[k].imag = 0.0; }
      for (k=npts; k<FFT_LEN; k++) { b[k].real = 0.0; b[k].imag = 0.0; }
      fftwf_execute(b_longf);

      sprintf(tmp,"synth_chirp%4i.freq",best);
      fp = fopen(tmp,"w");
      for (k=0; k<FFT_LEN; k++) {
         double tbuf = sqrt(b[k].real*b[k].real+b[k].imag*b[k].imag)/(double)FFT_LEN;
         fprintf(fp,"%i %lf\n",k,tbuf);
      }
      fclose(fp);
      first = 0;
    }

    min_diff = 10000000.0;
   
/*    
   for (j=6000; j<8200; j++) {
     diff[j] = 0.0;
     for (i=0; i<npts; i++) diff[j] += fbuf[i+j] - max*ref[i];
     if (fabs(diff[j]) < min_diff) { 
       min_diff = fabs(diff[j]); 
	best = j; 
 printf("reset min_diff to %lf at %i\n",min_diff,best);
    }
  }
    
     for (j=6000; j<8200; j++) printf("J = %i:  Diff is %lf\n",j,diff[j]);
    printf("Best fit is at an offset of %i (%lf)\n",best,min_diff);
*/
    
/*    best = 6384;  */
    for (j=0; j<SAMPLES_PER_LINE; j++) {
      if (j<best) fbuf2[j] = fbuf[j];
      else if (j<best+npts) {
          fbuf2[j] = max*((fbuf[j]/max)-ref[j-best]);
	  if (fabs(fbuf2[j]) > 16) {
	  	// printf("max = %lf, ref = %lf, fbuf = %lf, fbuf2 = %lf\n",max,ref[j-best],fbuf[j],fbuf2[j]);
		if (fbuf2[j]<0) fbuf2[j] = -16.0;
		else fbuf2[j] = 16.0;
	  }
      }
      else fbuf2[j] = fbuf[j];
    }
}
      


