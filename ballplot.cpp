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
#include "r1.h"
#include "r2.h"
#include "r3.h"
#include "r4.h"
#include "r6.h"
#include "elliptic.h"
#include "larc.h"
#include "mca.h"
#include "euclid.h"

/*=====================================
 * function prototypes
 *=====================================*/
int plot_print_header(FILE *fptr, const time_t time1, const char *comment, FILE *lptr);
int plot_log_header(FILE *fptr);
int plot_dna_seq(const unsigned seed, int N, const char *filename);

/*-------------------------------------------------------------------------
 * compute metric ball with center <p> and radius <r> 
 * using Balloon metric (deprecated metric)
 *-------------------------------------------------------------------------*/
int plot_balloon_metric_ball(vectR2 p, double r, FILE *fptr){
  ellipsec ellipse;
  time_t time1,time2; 
  struct tm *gmt;
  double seconds;
  double x,y,d;
  vectR2 q;

  time(&time1);
  gmt = gmtime(&time1);

  fprintf(fptr,"%%========================================================================\n");
  fprintf(fptr,"%% Daniel J. Greenhoe \n");
  fprintf(fptr,"%% The command \\fileplot{<filename>} may be used in a LaTeX environment\n");
  fprintf(fptr,"%% for plotting data.\n");
  fprintf(fptr,"%% \\fileplot is available in the LaTeX PSTricks package.\n");
  fprintf(fptr,"%% References: http://www.ctan.org/pkg/pstricks-base\n");
  fprintf(fptr,"%%             http://www.ctan.org/pkg/pstricks-add\n");
//fprintf(fptr,"%% %4d %2d %2d %2d:%2d\n",gmt->tm_year,gmt->tm_mon,gmt->tm_hour,gmt->tm_hour,gmt->tm_min);
  fprintf(fptr,"%% GMT %s",asctime(gmt));
  fprintf(fptr,"%%========================================================================\n");
  fprintf(fptr,"[\n");


  for (x=-5;x<=5;x+=0.03){
    for (y=-5;y<=5;y+=0.03){
      q.put(x,y);
      d=metric_balloon(p,q);
      if(d<=r)fprintf(fptr,"(%9.6lf,%9.6lf) %%=(x,y); p=(%9.6lf,%9.6lf) d(p,q)=%9.6lf\n",q.getx(),q.gety(),p.getx(),p.gety(),d);
      }
    printf(":");
    }
  fprintf(fptr,"]\n");
  time(&time2);
  gmt = gmtime(&time2);
  seconds = difftime(time2,time1);
  fprintf(fptr,"%%========================================================================\n");
  fprintf(fptr,"%% end data \n");
  fprintf(fptr,"%% GMT %s",asctime(gmt));
  fprintf(fptr,"%% (%.0lf seconds ellapsed)\n",seconds);
  fprintf(fptr,"%%========================================================================\n");
  return 0;
}

/*-------------------------------------------------------------------------
 * compute Lebesgue arc metric ball in R^2 with center <p> and radius <r> 
 * using Lebesgue arc metric 
 *-------------------------------------------------------------------------*/
int plot_larc_ball(vectR2 p, double r, FILE *fptr){
  ellipsec ellipse;
  time_t time1,time2; 
  struct tm *gmt;
  double seconds;
  double x,y,d;
  vectR2 q;

  time(&time1);
  gmt = gmtime(&time1);

  fprintf(fptr,"%%========================================================================\n");
  fprintf(fptr,"%% Daniel J. Greenhoe \n");
  fprintf(fptr,"%% data file for plotting Lebesgue Arc metric ball\n");
  fprintf(fptr,"%% The command \\fileplot{<filename>} may be used in a LaTeX environment\n");
  fprintf(fptr,"%% for plotting data.\n");
  fprintf(fptr,"%% \\fileplot is available in the LaTeX PSTricks package.\n");
  fprintf(fptr,"%% Reference: http://www.ctan.org/pkg/pstricks-base\n");
//fprintf(fptr,"%% %4d %2d %2d %2d:%2d\n",gmt->tm_year,gmt->tm_mon,gmt->tm_hour,gmt->tm_hour,gmt->tm_min);
  fprintf(fptr,"%% GMT %s",asctime(gmt));
  fprintf(fptr,"%%========================================================================\n");
  fprintf(fptr,"[\n");


  for (x=-3;x<=5;x+=0.025){
    for (y=-3;y<=5;y+=0.025){
      q.put(x,y);
      d=larc_metric(p,q);
      if(d<=r)fprintf(fptr,"(%9.6lf,%9.6lf) %%=(x,y); p=(%9.6lf,%9.6lf) d(p,q)=%9.6lf\n",q.getx(),q.gety(),p.getx(),p.gety(),d);
      }
    printf(":");
    }
  fprintf(fptr,"]\n");
  time(&time2);
  gmt = gmtime(&time2);
  seconds = difftime(time2,time1);
  fprintf(fptr,"%%========================================================================\n");
  fprintf(fptr,"%% end data \n");
  fprintf(fptr,"%% GMT %s",asctime(gmt));
  fprintf(fptr,"%% (%.0lf seconds ellapsed)\n",seconds);
  fprintf(fptr,"%%========================================================================\n");
  return 0;
}

/*-------------------------------------------------------------------------
 * compute Lebesgue arc metric ball in R^3 with center <p> and radius <r> 
 * using Lebesgue arc metric 
 *-------------------------------------------------------------------------*/
int plot_larc_ball(vectR3 p, double r, FILE *fptr){
  time_t time1,time2; 
  struct tm *gmt;
  double seconds;
  double theta,phi;
  vectR3 q;
  const double minrq    = 0.0;
  const double maxrq    = 1.0e6;
  const double maxerror = 1.0e6;
  const long int N      = 1000;

  time(&time1);
  gmt = gmtime(&time1);

  fprintf(fptr,"%%========================================================================\n");
  fprintf(fptr,"%% Daniel J. Greenhoe \n");
  fprintf(fptr,"%% data file for plotting Lebesgue Arc metric ball in R^3\n");
  fprintf(fptr,"%% The tikzpicture environment may be used in LaTeX\n");
  fprintf(fptr,"%% for plotting data.\n");
  fprintf(fptr,"%% tikzpicture is available in the LaTeX pgfplots package.\n");
  fprintf(fptr,"%% Reference:\n");
  fprintf(fptr,"%% ftp://ftp.ccu.edu.tw/pub/tex/graphics/pgf/contrib/pgfplots/doc/pgfplots.pdf\n");
  fprintf(fptr,"%% GMT %s",asctime(gmt));
  fprintf(fptr,"%%========================================================================\n");


  for(phi=PI/2; phi>=-PI/2; phi-=PI/18){
    for(theta=0; theta<2*PI; theta+=2*PI/36){
    //q = larc_findq( p, theta, phi, r, 1000);
      q = larc_findq( p, theta, phi, r, minrq, maxrq, maxerror, N );
      fprintf(fptr,"%9.6lf %9.6lf %9.6lf %%=q; p=(%9.6lf,%9.6lf,%9.6lf)\n",q.getx(),q.gety(),q.getz(),p.getx(),p.gety(),p.getz());
      }
    fprintf(fptr,"\n");
    printf(":");
    }
  time(&time2);
  gmt = gmtime(&time2);
  seconds = difftime(time2,time1);
  fprintf(fptr,"%%========================================================================\n");
  fprintf(fptr,"%% end data \n");
  fprintf(fptr,"%% GMT %s",asctime(gmt));
  fprintf(fptr,"%% (%.0lf seconds ellapsed)\n",seconds);
  fprintf(fptr,"%%========================================================================\n");
  return 0;
}

/*-------------------------------------------------------------------------
 * compute alpha-scaled Euclidean metric ball in R^3 
 * with center <p> and radius <r> using Euclidean metric 
 *-------------------------------------------------------------------------*/
int plot_euclidean_metric_ball(double alpha, vectR3 p, double r, FILE *fptr){
  time_t time1,time2; 
  struct tm *gmt;
  double seconds;
  double theta,phi;
  vectR3 q;

  time(&time1);
  gmt = gmtime(&time1);

  fprintf(fptr,"%%========================================================================\n");
  fprintf(fptr,"%% Daniel J. Greenhoe \n");
  fprintf(fptr,"%% data file for plotting an %lf-scaled Euclidean metric ball in R^3\n",alpha);
  fprintf(fptr,"%% The tikzpicture environment may be used in LaTeX\n");
  fprintf(fptr,"%% for plotting data.\n");
  fprintf(fptr,"%% tikzpicture is available in the LaTeX pgfplots package.\n");
  fprintf(fptr,"%% Reference:\n");
  fprintf(fptr,"%% ftp://ftp.ccu.edu.tw/pub/tex/graphics/pgf/contrib/pgfplots/doc/pgfplots.pdf\n");
  fprintf(fptr,"%% GMT %s",asctime(gmt));
  fprintf(fptr,"%%========================================================================\n");
  //fprintf(fptr,"\n");


  for(phi=PI/2; phi>=-PI/2; phi-=PI/18){
    for(theta=0; theta<2*PI; theta+=2*PI/36){
      q=ae_findq(alpha,p,theta,phi,r,1000);
      //fprintf(fptr,"%9.6lf %9.6lf %9.6lf\n",q.getx(),q.gety(),q.getz());
      fprintf(fptr,"%9.6lf %9.6lf %9.6lf %%=q; p=(%9.6lf,%9.6lf,%9.6lf)\n",q.getx(),q.gety(),q.getz(),p.getx(),p.gety(),p.getz());
      }
    fprintf(fptr,"\n");
    printf(":");
    }
  time(&time2);
  gmt = gmtime(&time2);
  seconds = difftime(time2,time1);
  fprintf(fptr,"%%========================================================================\n");
  fprintf(fptr,"%% end data \n");
  fprintf(fptr,"%% GMT %s",asctime(gmt));
  fprintf(fptr,"%% (%.0lf seconds ellapsed)\n",seconds);
  fprintf(fptr,"%%========================================================================\n");
  return 0;
}

/*-------------------------------------------------------------------------
 * compute mean circular arc metric ball in R^2 with center <p> and radius <r> 
 *-------------------------------------------------------------------------*/
int plot_mca_metric_ball(vectR2 p, double r, FILE *fptr){
  ellipsec ellipse;
  time_t time1,time2; 
  struct tm *gmt;
  double seconds;
  double x,y,d;
  vectR2 q;

  time(&time1);
  gmt = gmtime(&time1);

  fprintf(fptr,"%%========================================================================\n");
  fprintf(fptr,"%% Daniel J. Greenhoe \n");
  fprintf(fptr,"%% data file for plotting Mean Circular Arc metric ball\n");
  fprintf(fptr,"%% The command \\fileplot{<filename>} may be used in a LaTeX environment\n");
  fprintf(fptr,"%% for plotting data.\n");
  fprintf(fptr,"%% \\fileplot is available in the LaTeX PSTricks package.\n");
  fprintf(fptr,"%% Reference: http://www.ctan.org/pkg/pstricks-base\n");
  fprintf(fptr,"%% GMT %s",asctime(gmt));
  fprintf(fptr,"%%========================================================================\n");
  fprintf(fptr,"[\n");


  for (x=-3;x<=5;x+=0.025){
    for (y=-3;y<=5;y+=0.025){
      q.put(x,y);
      d=mca_metric(p,q);
      if(d<=r)fprintf(fptr,"(%9.6lf,%9.6lf) %%=(x,y); p=(%9.6lf,%9.6lf) d(p,q)=%9.6lf\n",q.getx(),q.gety(),p.getx(),p.gety(),d);
      }
    printf(":");
    }
  fprintf(fptr,"]\n");
  time(&time2);
  gmt = gmtime(&time2);
  seconds = difftime(time2,time1);
  fprintf(fptr,"%%========================================================================\n");
  fprintf(fptr,"%% end data \n");
  fprintf(fptr,"%% GMT %s",asctime(gmt));
  fprintf(fptr,"%% (%.0lf seconds ellapsed)\n",seconds);
  fprintf(fptr,"%%========================================================================\n");
  return 0;
}




