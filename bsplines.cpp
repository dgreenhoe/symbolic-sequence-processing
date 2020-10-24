/*============================================================================
 * Daniel J. Greenhoe
 *============================================================================*/
/*=====================================
 * headers
 *=====================================*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "main.h"

/*---------------------------------------------------------------------------
 * plot auto-power spectrum data for plotting 
 * within LaTeX ps-tricks environment
 *---------------------------------------------------------------------------*/
int bspline_Sdat(void){
  time_t time1; time(&time1);   //starting time stamp (passed to plotting routine)
  struct tm *gmt;
  gmt = gmtime(&time1);
  char tbuffer[80];
  const double n=3;      // order of B-spline
  const long M=1024;     // approximate number of data points
  const long N=1000000;  // number of iterations per data point
  double w=0,s=0;
  long k=0;
  strftime(tbuffer,80,"%Y %B %d %A %I:%M:%S %p UTC",gmt);
  printf("%%============================================================================\n");
  printf("%% Daniel J. Greenhoe\n");
  printf("%% n=%1.0lf order B-spline auto-power spectrum data\n",n);
  printf("%% M = number of data points approx= %ld\n",M);
  printf("%% N = iterations per data point = %ld\n",N);
  printf("%% %s\n",tbuffer);
  printf("%% The command \fileplot{d4_phi.dat} may be used in a LaTeX environment for plotting data.\n");
  printf("%% \\fileplot is available in the LaTeX PSTricks package.\n");
  printf("%% Reference: http://www.ctan.org/pkg/pstricks\n");
  printf("%%============================================================================\n");
  printf("[\n");
  for(w=-8.0;w<=8.0;w+=16./(double)M){
    s=0;
    for(k=1; k<=N; k++){
      s += pow(1.0/(2*k-w/PI),2.*n+2.);
      s += pow(1.0/(2*k+w/PI),2.*n+2.);
      }
    s*=pow(sin(w/2.)/(PI/2.),2.*n+2.);
    if(w==0) s+=1; // use l'Hopital's rule
    else     s+=pow(sin(w/2.)/(w/2.),2.*n+2.);
    printf("(%lf, %lf)\n",w,s);
    }
  printf("]\n");
  return 0;
  }
