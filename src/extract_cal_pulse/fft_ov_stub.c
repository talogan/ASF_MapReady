
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define LEN 16
#define M_PI 3.14159265359
double *fft_oversamp(double *in, int len, int exp_factor);

main (int argc, char *argv[])
{

   double x[LEN];
   int i, k;
   FILE *fpout; 
   
   double *ret;
   double z;
   
   for (i=0; i< LEN; i++)
    {
       z = sin(((double)i/(double)LEN)*2.0*M_PI)*16.0+16.0;
       x[i] = z;
       printf("x[%i] = %f\n",i,x[i]);
    }
    
   fpout = fopen("a.org","w");
    for (k=0; k<LEN; k++) {
       fprintf(fpout,"%i %f\n",k,x[k]);
    }
   fclose(fpout);

   printf("CALLING OVERSAMP 1\n");
   ret = fft_oversamp(x,LEN,1);
   fpout = fopen("a1.new","w");
    for (k=0; k<LEN; k++) {
       fprintf(fpout,"%i %lf\n",k,ret[k]);
    }
   fclose(fpout);
   free(ret);
   
   printf("CALLING OVERSAMP 2\n");
   ret = fft_oversamp(x,LEN,2);
   fpout = fopen("a2.new","w");
   for (k=0; k<LEN*2; k++) {
       fprintf(fpout,"%i %lf\n",k,ret[k]);
    }
   fclose(fpout);
   free(ret);
      
   printf("CALLING OVERSAMP 4\n");
   ret = fft_oversamp(x,LEN,4);
   fpout = fopen("a4.new","w");
    for (k=0; k<LEN*4; k++) {
       fprintf(fpout,"%i %lf\n",k,ret[k]);
    }
   fclose(fpout);
   free(ret);
   
   printf("CALLING OVERSAMP 8\n");
   ret = fft_oversamp(x,LEN,8);
   fpout = fopen("a8.new","w");
    for (k=0; k<LEN*8; k++) {
       fprintf(fpout,"%i %lf\n",k,ret[k]);
    }
   fclose(fpout);
   free(ret);
}

   
   
