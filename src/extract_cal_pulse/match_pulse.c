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

#define SAMPLES_PER_LINE   13680	/* decoded samples per output line */
#define EXPANSION_FACTOR   4

void match_pulse(double *fbuf,double *ovpulse,double *fbuf2)
{
    double dtmp[SAMPLES_PER_LINE];
    int i, j, k, best;
    double diff[EXPANSION_FACTOR*8];
    double min_diff;
    FILE *fpout;
    static int first = 1;

    min_diff = 10000000.0;
    best = 0;
    for (j=0; j<EXPANSION_FACTOR*8; j++) {
      diff[j] = 0.0;
      for (i=0; i<SAMPLES_PER_LINE-10; i++) {
        diff[j] += fbuf[i] - ovpulse[i*EXPANSION_FACTOR+j];
      }
//      printf("MIN DIFF = %lf\n",min_diff);
      
      if (diff[j] < min_diff) { 
        min_diff = diff[j]; 
	best = j; 
//	printf("reset min_diff to %lf at %i\n",min_diff,best);
      }
    }
    
//    for (j=0; j<EXPANSION_FACTOR*8; j++) {
//	printf("J = %i:  Diff is %lf\n",j,diff[j]);
//    }
    
//    printf("Best fit is at an offset of %i (%lf)\n",best,min_diff);
    
    if (first==1) fpout = fopen("downsampled_chirp.txt","w");
    for (i=0; i<SAMPLES_PER_LINE; i++) {
      fbuf2[i] = fbuf[i] - ovpulse[i*EXPANSION_FACTOR+best];
      if (first==1) fprintf(fpout,"%lf\n",ovpulse[i*EXPANSION_FACTOR+best]);
    }
    if (first==1) {fclose(fpout); first = 0;}
}
      

