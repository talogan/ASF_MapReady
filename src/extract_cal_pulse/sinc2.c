

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define  NPTS  1700
#define  pi     3.14159
main (int argc, char *argv[])
{

double x, y;
int i,j;
double ref[NPTS];
double r001 = 848968.28241700004;
int isave = 1;
int nextend = 0;
int nlinesaz = 6840;
double wavl = 0.23499999999999999;
double azres = 4.0000;
double vel1 = 7617.9448;
double prf1 = 1647.0;
double fs = 22765000.0;
double re = 6375971.90;
double ht1 = 805223.87231400004;
double delr=299792456./fs/2.;
double r01 = r001 + (isave-nextend-1) * delr;
double r1  = r01 + delr*(nlinesaz-1);
double dxsamp1 = vel1/prf1;
double npfin1 = r1*wavl/(2.0*azres*dxsamp1)+2;
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


printf("npts = %i\n",npts);

for (i=-npts/2; i<npts/2; i++)
 {
   k=k+1;
   t=(double)i*ts;
   phase = pi*slope1*t*t+pi*fs*t;
   printf("phase %lf\n",phase);
   ref[i+npts/2+1]=cos(phase);
 }

FILE *fp = fopen("sinc.txt","w");
for (i=0; i<npts; i++) fprintf(fp,"%lf\n",ref[i]);
fclose(fp);

}


