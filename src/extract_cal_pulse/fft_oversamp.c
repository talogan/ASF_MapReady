/******************************************************************************
NAME: FFT Oversampling Subroutine

SYNOPSIS:   double* fft_oversamp(double *in, int len, int exp_factor)
	
	Given double array in of len, return oversampled
	double array of length (fft_len*exp_factor), oversampled 
	using FFTs and padding in the frequency domain.

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
    1.0	    9/16/13  T. Logan     Seasat project
    
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

typedef struct {
   float real;
   float imag;
} complexFloat;

double *fft_oversamp(double *in, int len, int exp_factor, double *out)
{
    static complexFloat *a, *b;
    static fftwf_plan   a_short, b_longb, b_longf;
    static int outlen;
    static int fft_len;
    static first = 1;

    int k;
    double kk;
    FILE *fpout;

    if (first == 1) { 
      printf("fft_oversamp: oversampling input of length %i to length %i\n",len,len*exp_factor);
      fft_len = 2;
      while (fft_len < len) fft_len *=2;
      printf("fft_oversamp: using ffts of length %i\n",fft_len);
      outlen = fft_len * exp_factor;
      printf("fft_oversamp: creating output of length %i\n",outlen);
    
      a = (complexFloat *) fftwf_malloc(sizeof(fftwf_complex)*fft_len);
      a_short = fftwf_plan_dft_1d(fft_len, (fftwf_complex *) a, (fftwf_complex *)a, FFTW_FORWARD, FFTW_MEASURE);

      b = (complexFloat *) fftwf_malloc(sizeof(fftwf_complex)*outlen);
      b_longb = fftwf_plan_dft_1d(outlen, (fftwf_complex *) b, (fftwf_complex *)b, FFTW_BACKWARD, FFTW_MEASURE);
      b_longf = fftwf_plan_dft_1d(outlen, (fftwf_complex *) b, (fftwf_complex *)b, FFTW_FORWARD, FFTW_MEASURE);
    }


if (first==1) {
   fpout = fopen("a.time","w");
   for (k=0; k<len; k++) {
      double mag = in[k];
      fprintf(fpout,"%i %lf\n",k,mag);
   }
   fclose(fpout);
}
    
    /* convert 'in' to complex, store in a, transform */
    for (k=0; k<len; k++) { kk = in[k]; a[k].real = kk; a[k].imag = 0.0; }
    for (k=len; k<fft_len; k++) { a[k].real = 0.0; a[k].imag = 0.0; }
    
    fftwf_execute(a_short);

if (first==1) {
   fpout = fopen("a.trans","w");
   for (k=0; k<fft_len; k++) {
      double mag = sqrt(a[k].real*a[k].real+a[k].imag*a[k].imag);
      fprintf(fpout,"%i %lf\n",k,mag);
   }
   fclose(fpout);
}

    for (k=0; k<outlen; k++) { b[k].real = 0.0; b[k].imag = 0.0; }

    /* Copy first half of a into b */
    for (k=0; k<fft_len/2; k++) { b[k].real = a[k].real; b[k].imag = a[k].imag;}

    /* Copy second half of a into b after the padding */
    for (k=fft_len/2; k<fft_len; k++) {
      b[k+fft_len*(exp_factor-1)].real = a[k].real;
      b[k+fft_len*(exp_factor-1)].imag = a[k].imag;
    }

if (first==1) {
    fpout = fopen("b.freq","w");
    for (k=0; k<len*exp_factor; k++) {
      out[k] = sqrt(b[k].real*b[k].real+b[k].imag*b[k].imag);
      fprintf(fpout,"%i %lf\n",k,out[k]);
    }
    fclose(fpout);
}

    /* perform the actual oversampling */
    fftwf_execute(b_longb);

if (first==1) fpout = fopen("b.time","w");
    for (k=0; k<len*exp_factor; k++) {
       out[k] = sqrt(b[k].real*b[k].real+b[k].imag*b[k].imag)/(double)len;
if (first==1)       fprintf(fpout,"%i %lf\n",k,out[k]);
    }
if (first==1)    fclose(fpout);


   first = 0;

}
