/*============================================================================
 * Daniel J. Greenhoe
 * elliptic routines
 * "ellipse" here is defined as all the points (x,y) in R^2 that satisfy
 *   x^2    y^2   
 *   ---  + ---  = 1  
 *   a^2    b^2
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
#include "symseq.h"
#include "r1.h"
#include "r2.h"
#include "r3.h"
#include "r4.h"
#include "r6.h"
#include "c1.h"
#include "c4.h"
#include "c6.h"
#include "r1op.h"
#include "r2op.h"
#include "r6op.h"
#include "elliptic.h"
#include "larc.h"
#include "mca.h"
#include "euclid.h"
#include "die.h"
#include "realdie.h"
#include "fairdie.h"
#include "spinner.h"
#include "dnan.h"
#include "dft.h"
#include "test.h"
#include "gtest/gtest.h"  // https://github.com/google/googletest/blob/master/googletest/docs/primer.md

//-----------------------------------------------------------------------------
//! \brief Test complex operations
//-----------------------------------------------------------------------------
TEST( TestSuiteGeneral, C1 )
{
  complex p(3,-4);
  complex q,s;

  ASSERT_EQ( p.getx(),  3 );
  ASSERT_EQ( p.gety(), -4 );
  ASSERT_EQ( p.norm(),  5 );
  p.put(sqrt(3), 1);
  ASSERT_DOUBLE_EQ( p.theta()/M_PI*180, 30.0 );
  q=-p;
  ASSERT_EQ( q.getx(), -sqrt(3) );
  ASSERT_EQ( q.gety(), -1       );
  p.put(  3, -5 ); 
  q.put( -2,  7 ); 
  s = p + q;
  ASSERT_EQ( s.getx(),  3 - 2 );
  ASSERT_EQ( s.gety(), -5 + 7 );
  s = p - q;
  ASSERT_EQ( s.getx(),  3 + 2 );
  ASSERT_EQ( s.gety(), -5 - 7 );
  p.put(  2, -3 ); 
  q.put( -5,  7 ); 
  s = p * q;
  ASSERT_EQ( s.getx(),  (2 * (-5)) - ((-3) * ( 7)) );
  ASSERT_EQ( s.gety(),  (2 * ( 7)) + ((-5) * (-3)) );
  q = p; 
  q &= (PI/2);
  ASSERT_DOUBLE_EQ( q.getx(),  3 );
  ASSERT_DOUBLE_EQ( q.gety(),  2 );
  q = p; 
  q &= M_PI;
  ASSERT_DOUBLE_EQ( q.getx(), -2 );
  ASSERT_DOUBLE_EQ( q.gety(),  3 );
  p.clear();
  ASSERT_EQ( p.getx(),  0 );
  ASSERT_EQ( p.gety(),  0 );
}

/*-------------------------------------------------------------------------
 * test complex operations
 *-------------------------------------------------------------------------*/
int test_seqR2(void){
  printf("\ntest seqR2 operations:\n");
  printf("-----------------------\n");

  seqR2 x(6);       printf("construct: "); x.list();
  x.fill(2);        printf("\nfill(2): "); x.list();
  x.inc(2,-5,1,-5); printf("\ninc: ");     x.list();
  x.max('p');
  x.clear();        printf("\nclear: ");   x.list();
  return 0;
  }


/*-------------------------------------------------------------------------
 * test Lagrange arc metric in R^2
 *-------------------------------------------------------------------------*/
int test_pqtheta(void){
  vectR2 p,q;
  double tr,td,dc;
  int fails=0;
  printf("\n(p,q) theta tests in R2\n");
  printf("-----------------------\n");
  p.put(0,0);                   q.put(2,0);                    tr=pqtheta(p,q); td=tr*180/PI; dc= -1*180/PI; printf("p=(%6.3lf,%6.3lf) q=(%6.3lf,%6.3lf) theta=%6.2lf radians %6.2lf degrees ",p.getx(),p.gety(),q.getx(),q.gety(),tr,td);if(td>dc*0.999||td<dc*1.001){fails++; printf("FAIL %lf!=%lf\n",td,dc);}else printf("ok\n");
  p.put(0,1);                   q.put(0,0);                    tr=pqtheta(p,q); td=tr*180/PI; dc= -2*180/PI; printf("p=(%6.3lf,%6.3lf) q=(%6.3lf,%6.3lf) theta=%6.2lf radians %6.2lf degrees ",p.getx(),p.gety(),q.getx(),q.gety(),tr,td);if(td>dc*0.999||td<dc*1.001){fails++; printf("FAIL %lf!=%lf\n",td,dc);}else printf("ok\n");
  p.put(0,2);                   q.put(2,0);                    tr=pqtheta(p,q); td=tr*180/PI; dc= 90; printf("p=(%6.3lf,%6.3lf) q=(%6.3lf,%6.3lf) theta=%6.2lf radians %6.2lf degrees ",p.getx(),p.gety(),q.getx(),q.gety(),tr,td);if(td<dc*0.999||td>dc*1.001){fails++; printf("FAIL %lf!=%lf\n",td,dc);}else printf("ok\n");
  p.put(0,1);                   q.put(1,0);                    tr=pqtheta(p,q); td=tr*180/PI; dc= 90; printf("p=(%6.3lf,%6.3lf) q=(%6.3lf,%6.3lf) theta=%6.2lf radians %6.2lf degrees ",p.getx(),p.gety(),q.getx(),q.gety(),tr,td);if(td<dc*0.999||td>dc*1.001){fails++; printf("FAIL %lf!=%lf\n",td,dc);}else printf("ok\n");
  p.put(1,3);                   q.put(-1,-3);                  tr=pqtheta(p,q); td=tr*180/PI; dc=180; printf("p=(%6.3lf,%6.3lf) q=(%6.3lf,%6.3lf) theta=%6.2lf radians %6.2lf degrees ",p.getx(),p.gety(),q.getx(),q.gety(),tr,td);if(td<dc*0.999||td>dc*1.001){fails++; printf("FAIL %lf!=%lf\n",td,dc);}else printf("ok\n");
  p.put(-2,5);                  q.put(2,-5);                   tr=pqtheta(p,q); td=tr*180/PI; dc=180; printf("p=(%6.3lf,%6.3lf) q=(%6.3lf,%6.3lf) theta=%6.2lf radians %6.2lf degrees ",p.getx(),p.gety(),q.getx(),q.gety(),tr,td);if(td<dc*0.999||td>dc*1.001){fails++; printf("FAIL %lf!=%lf\n",td,dc);}else printf("ok\n");
  p.put(0,1);                   q.put( cos(PI/4), sin(PI/4));  tr=pqtheta(p,q); td=tr*180/PI; dc= 45; printf("p=(%6.3lf,%6.3lf) q=(%6.3lf,%6.3lf) theta=%6.2lf radians %6.2lf degrees ",p.getx(),p.gety(),q.getx(),q.gety(),tr,td);if(td<dc*0.999||td>dc*1.001){fails++; printf("FAIL %lf!=%lf\n",td,dc);}else printf("ok\n");
  p.put(0,1);                   q.put( cos(PI/4),-sin(PI/4));  tr=pqtheta(p,q); td=tr*180/PI; dc=135; printf("p=(%6.3lf,%6.3lf) q=(%6.3lf,%6.3lf) theta=%6.2lf radians %6.2lf degrees ",p.getx(),p.gety(),q.getx(),q.gety(),tr,td);if(td<dc*0.999||td>dc*1.001){fails++; printf("FAIL %lf!=%lf\n",td,dc);}else printf("ok\n");
  p.put(0,1);                   q.put(-cos(PI/4), sin(PI/4));  tr=pqtheta(p,q); td=tr*180/PI; dc= 45; printf("p=(%6.3lf,%6.3lf) q=(%6.3lf,%6.3lf) theta=%6.2lf radians %6.2lf degrees ",p.getx(),p.gety(),q.getx(),q.gety(),tr,td);if(td<dc*0.999||td>dc*1.001){fails++; printf("FAIL %lf!=%lf\n",td,dc);}else printf("ok\n");
  p.put(0,1);                   q.put(-cos(PI/4),-sin(PI/4));  tr=pqtheta(p,q); td=tr*180/PI; dc=135; printf("p=(%6.3lf,%6.3lf) q=(%6.3lf,%6.3lf) theta=%6.2lf radians %6.2lf degrees ",p.getx(),p.gety(),q.getx(),q.gety(),tr,td);if(td<dc*0.999||td>dc*1.001){fails++; printf("FAIL %lf!=%lf\n",td,dc);}else printf("ok\n");
  p.put( 7*cos(PI/4), 7*sin(PI/4)); q.put(0,1);                    tr=pqtheta(p,q); td=tr*180/PI; dc= 45; printf("p=(%6.3lf,%6.3lf) q=(%6.3lf,%6.3lf) theta=%6.2lf radians %6.2lf degrees ",p.getx(),p.gety(),q.getx(),q.gety(),tr,td);if(td<dc*0.999||td>dc*1.001){fails++; printf("FAIL %lf!=%lf\n",td,dc);}else printf("ok\n");
  p.put( 7*cos(PI/4),-7*sin(PI/4)); q.put(0,1);                    tr=pqtheta(p,q); td=tr*180/PI; dc=135; printf("p=(%6.3lf,%6.3lf) q=(%6.3lf,%6.3lf) theta=%6.2lf radians %6.2lf degrees ",p.getx(),p.gety(),q.getx(),q.gety(),tr,td);if(td<dc*0.999||td>dc*1.001){fails++; printf("FAIL %lf!=%lf\n",td,dc);}else printf("ok\n");
  p.put(-7*cos(PI/4), 7*sin(PI/4)); q.put(0,1);                    tr=pqtheta(p,q); td=tr*180/PI; dc= 45; printf("p=(%6.3lf,%6.3lf) q=(%6.3lf,%6.3lf) theta=%6.2lf radians %6.2lf degrees ",p.getx(),p.gety(),q.getx(),q.gety(),tr,td);if(td<dc*0.999||td>dc*1.001){fails++; printf("FAIL %lf!=%lf\n",td,dc);}else printf("ok\n");
  p.put(-7*cos(PI/4),-7*sin(PI/4)); q.put(0,1);                    tr=pqtheta(p,q); td=tr*180/PI; dc=135; printf("p=(%6.3lf,%6.3lf) q=(%6.3lf,%6.3lf) theta=%6.2lf radians %6.2lf degrees ",p.getx(),p.gety(),q.getx(),q.gety(),tr,td);if(td<dc*0.999||td>dc*1.001){fails++; printf("FAIL %lf!=%lf\n",td,dc);}else printf("ok\n");
  printf("number of fails = %d\n",fails);
  return fails;
  }

/*-------------------------------------------------------------------------
 * half circle test
 *-------------------------------------------------------------------------*/
int test_halfcircle(void){
  ellipsec circle(1,1);
  vectR2 q;
  vectR2 (ellipsec::*memfptr)(double t);
  double length, errorl,perror;
  
  memfptr = &ellipsec::xy;
  q=(circle.*memfptr)(1); printf("q=(%lf,%lf)\n",q.getx(),q.gety());

  length = circle.pathlength(PI/2,-PI/2, 100000);
  errorl  = length-PI;
  perror  = fabs(100*errorl/PI);
  printf("half circle = %.16lf\n",length);
  printf("pi          = %.16lf error=%lf (%lf%%)\n",PI, errorl, perror);
  if (perror<0.001)  return 1;
  else             return 0;
  }

/*-------------------------------------------------------------------------
 * perimeter of circle test
 *-------------------------------------------------------------------------*/
int test_circle(void){
  ellipsec circle(1,1);
  double x,y;
  double length, errorl,perror;
  int fails=0;
  long N=1000;

  length = circle.pathlength(0,2*PI,N);
  errorl  = length-2*PI;
  perror  = fabs(100*errorl/(2*PI));
  printf("\ntest circle.pathlength(0,2pi,%ld):\n",N);
  printf("---------------------------------------\n");

  x=0; y=0;
  circle.set(1,1,0,x,y);
  printf("L circle(%6.3lf,%6.3lf)=%.16lf  error=%lf (%lf%%)",x,y,length, errorl, perror); 
  if (perror<0.001)printf("  ok\n"); 
  else{fails++;  printf("  FAIL\n");}

  x=1; y=-2;
  circle.set(1,1,0,x,y);
  printf("L circle(%6.3lf,%6.3lf)=%.16lf  error=%lf (%lf%%)",x,y,length, errorl, perror); 
  if (perror<0.001)printf("  ok\n"); 
  else{fails++;  printf("  FAIL\n");}

  x=-3; y=-5;
  circle.set(1,1,0,x,y);
  printf("L circle(%6.3lf,%6.3lf)=%.16lf  error=%lf (%lf%%)",x,y,length, errorl, perror); 
  if (perror<0.001)printf("  ok\n"); 
  else{fails++;  printf("  FAIL\n");}
  return fails;
  }

/*-------------------------------------------------------------------------
 * find points on unit circle that are distance 1 
 * from the point at pi/4 radians on unit circle
 *-------------------------------------------------------------------------*/
int test_circle_d1(void){
  ellipsec circle(1,1);
  double phi1, phi2;
  double err1, err2;
  double error1, error2;
  double perror1, perror2;

  printf("\nfind points on unit circle that are distance 1 from the point at pi/4 radians on unit circle\n");
  printf("--------------------------------------------------------------\n");
  circle.findt_dfroms(PI/4, 1, +1, 1000, &phi1, &err1);
  circle.findt_dfroms(PI/4, 1, -1, 1000, &phi2, &err2);
  error1 = phi1 - (PI/4+1);
  error2 = phi2 - (PI/4-1);
  perror1 = fabs(100*error1/(PI/4+1));
  perror2 = fabs(100*error2/(PI/4-1));
  printf("(x,y) = (%9.6lf, %9.6lf) phi = %9.6lf error=%9.6lf(%lf%%)\n",circle.x(phi1),circle.y(phi1),phi1,error1,perror1);
  printf("(x,y) = (%9.6lf, %9.6lf) phi = %9.6lf error=%9.6lf(%lf%%)\n",circle.x(phi2),circle.y(phi2),phi2,error2,perror2);
  if(perror1>0.1) return -1;
  if(perror2>0.5) return -5;
  return 0;
  }

/*-------------------------------------------------------------------------
 * test elliptic normalize
 *-------------------------------------------------------------------------*/
int test_normalize(void){
  vectR2 p,q;
  double phi;
  p.put(0,2); phi=PI/2; q= p & phi;  printf("p=(%lf,%lf) q=(%lf,%lf) phi=%lf\n",p.getx(),p.gety(),q.getx(),q.gety(), phi);
  p.put(0,2); phi=PI/2; q= p & (-phi);  printf("p=(%lf,%lf) q=(%lf,%lf) phi=%lf\n",p.getx(),p.gety(),q.getx(),q.gety(), -phi);
  p.put(2,0); phi=PI/2; q= p & phi;  printf("p=(%lf,%lf) q=(%lf,%lf) phi=%lf\n",p.getx(),p.gety(),q.getx(),q.gety(), phi);
  p.put(2,0); phi=PI/2; q= p & (-phi);  printf("p=(%lf,%lf) q=(%lf,%lf) phi=%lf\n",p.getx(),p.gety(),q.getx(),q.gety(), -phi);
  p.put(2,0); phi=PI/2; p &= phi;  printf("p=(%lf,%lf) phi=%lf\n",p.getx(),p.gety(), phi);
  return 0;
  }

/*-------------------------------------------------------------------------
 * find points on ellipse(0.8,1) that are distance 1 
 * from the point at pi/2 radians on the ellipse
 *-------------------------------------------------------------------------*/
int test_ellipse_d1(void){
  ellipsec ellipse(0.8,1);
  double phi1, phi2;
  double err1, err2;
  double perror1, perror2;

  printf("find points on ellipse(0.8,1) that are distance 1 from the point at pi/2 radians on the ellipse:\n");
  printf("-------------------------------------\n");
  ellipse.findt_dfroms(PI/2, 1, +1, 1000, &phi1, &err1);
  ellipse.findt_dfroms(PI/2, 1, -1, 1000, &phi2, &err2);
  perror1 = fabs(100*err1);
  perror2 = fabs(100*err2);
  printf("(x,y) = (%9.6lf, %9.6lf) phi = %9.6lf error = %9.6lf (%lf%%)\n", ellipse.x(phi1), ellipse.y(phi1), phi1, err1, perror1 );
  printf("(x,y) = (%9.6lf, %9.6lf) phi = %9.6lf error = %9.6lf (%lf%%)\n", ellipse.x(phi2), ellipse.y(phi2), phi2, err2, perror2 );
  if(perror1>0.2) return -1;
  if(perror2>0.2) return -2;
  return 0;
  }

/*-------------------------------------------------------------------------
 * find points on ellipse(0.8,1) that are distance 1 
 * from the point at pi/2 radians on the ellipse
 *-------------------------------------------------------------------------*/
int test_findt(void){
  vectR2 p(cos(PI/4),sin(PI/4));
  ellipsec ellipse;
  double t;
  if(ellipse.setab_givenxyb(p,2)){
    if(ellipse.tgivenxy(p,&t)==0)return 0;
    printf("t=%lf (x,y)=(%lf,%lf) ellipse(a,b)=(%lf,%lf)\n",t, p.getx(),p.gety(),ellipse.geta(),ellipse.getb());
    }
  return 1;
  }

/*-------------------------------------------------------------------------
 * find the perimeter of an ellipse
 *-------------------------------------------------------------------------*/
int test_perimeter(void){
  ellipsec ellipse(0.8,2);
  double perim=ellipse.perimeter(1000);
  double est  =ellipse.estimate();
  printf("ellipse(%lf,%lf) perimeter=%9.6lf estimate=%9.6lf\n",ellipse.geta(),ellipse.getb(),perim,est);
  return 1;
  }

/*-------------------------------------------------------------------------
 * test DNA metric
 *-------------------------------------------------------------------------*/
int test_dna_metric(void){
  const char domain[6]="0ATCG";
  const int N=strlen(domain);
  int n,m;
  double d;

  printf("dna_metric(a,b) test over domain of N=%d symbols\n",N);
  printf("------------------------------------------------\n");
  printf("  ");
  for(m=0;m<N;m++)printf("     %c",domain[m]);//print column headers: 0 A ...
  putchar('\n');
  for(n=0;n<N;n++){
    printf("%c  ",domain[n]);// print row headers
    for(m=0;m<N;m++){
      d=dnan_metric(domain[n],domain[m]);
      if(d==(double)abs(d)) printf("%5.0lf ",d);//if d is integer
      else                  printf("%5.3lf ",d);//if d is non-integer
      }
    putchar('\n');
    }
  return 1;
  }

/*-------------------------------------------------------------------------
 * test DNA-N metric
 *-------------------------------------------------------------------------*/
int test_dnan_metric(void){
  const char domain[7]="0ATCGN";
  const int N=strlen(domain);
  int n,m;
  double d;

  printf("dnan_metric(a,b) test over domain of N=%d symbols\n",N);
  printf("-------------------------------------------------\n");
  printf("  ");
  for(m=0;m<N;m++)printf("     %c",domain[m]);//print column headers: 0 A ...
  putchar('\n');
  for(n=0;n<N;n++){
    printf("%c  ",domain[n]);// print row headers
    for(m=0;m<N;m++){
      d=dnan_metric(domain[n],domain[m]);
      if(d==(double)abs(d)) printf("%5.0lf ",d);//if d is integer
      else                  printf("%5.3lf ",d);//if d is non-integer
      }
    putchar('\n');
    }
  return 1;
  }

/*-------------------------------------------------------------------------
 * test balloon metric
 *-------------------------------------------------------------------------*/
int test_balloon_metric(void){
  vectR2 p,q;
  double d;
  printf("Balloon metric tests\n");
  printf("--------------------\n");
  p.put(0,2);                   q.put(2,0);                    d=metric_balloon(p,q);  printf("d((%9.6lf,%9.6lf),(%9.6lf,%9.6lf))=%9.6lf\n",p.getx(),p.gety(),q.getx(),q.gety(),d);
  p.put(0,1);                   q.put(1,0);                    d=metric_balloon(p,q);  printf("d((%9.6lf,%9.6lf),(%9.6lf,%9.6lf))=%9.6lf\n",p.getx(),p.gety(),q.getx(),q.gety(),d);
  p.put(0,1);                   q.put( cos(PI/4), sin(PI/4));  d=metric_balloon(p,q);  printf("d((%9.6lf,%9.6lf),(%9.6lf,%9.6lf))=%9.6lf\n",p.getx(),p.gety(),q.getx(),q.gety(),d);
  p.put(0,1);                   q.put(-cos(PI/4), sin(PI/4));  d=metric_balloon(p,q);  printf("d((%9.6lf,%9.6lf),(%9.6lf,%9.6lf))=%9.6lf\n",p.getx(),p.gety(),q.getx(),q.gety(),d);
  p.put(0,1);                   q.put(-cos(PI/4),-sin(PI/4));  d=metric_balloon(p,q);  printf("d((%9.6lf,%9.6lf),(%9.6lf,%9.6lf))=%9.6lf\n",p.getx(),p.gety(),q.getx(),q.gety(),d);
  p.put(0,1);                   q.put( cos(PI/4),-sin(PI/4));  d=metric_balloon(p,q);  printf("d((%9.6lf,%9.6lf),(%9.6lf,%9.6lf))=%9.6lf\n",p.getx(),p.gety(),q.getx(),q.gety(),d);
  p.put(1,0);                   q.put(-0.5,-0.5);              d=metric_balloon(p,q);  printf("d((%9.6lf,%9.6lf),(%9.6lf,%9.6lf))=%9.6lf\n",p.getx(),p.gety(),q.getx(),q.gety(),d);
  p.put(1,0);                   q.put(-2,-2);                  d=metric_balloon(p,q);  printf("d((%9.6lf,%9.6lf),(%9.6lf,%9.6lf))=%9.6lf\n",p.getx(),p.gety(),q.getx(),q.gety(),d);
  p.put(1,0);                   q.put(0,2);                    d=metric_balloon(p,q);  printf("d((%9.6lf,%9.6lf),(%9.6lf,%9.6lf))=%9.6lf\n",p.getx(),p.gety(),q.getx(),q.gety(),d);
  p.put( cos(PI/4),-sin(PI/4)); q.put(-2,1);                   d=metric_balloon(p,q);  printf("d((%9.6lf,%9.6lf),(%9.6lf,%9.6lf))=%9.6lf\n",p.getx(),p.gety(),q.getx(),q.gety(),d);
  p.put( cos(PI/4),-sin(PI/4)); q.put(-1.63,1.33);             d=metric_balloon(p,q);  printf("d((%9.6lf,%9.6lf),(%9.6lf,%9.6lf))=%9.6lf\n",p.getx(),p.gety(),q.getx(),q.gety(),d);

  return 1;
  }

/*-------------------------------------------------------------------------
 * test Mean Cicular Arc metric in R^2
 *-------------------------------------------------------------------------*/
int test_mca_metric(void){
  vectR2 p,q;
  double d,d2;
  printf("Lagrange arc metric tests in R2\n");
  printf("--------------------\n");
  p.put(0,2);                   q.put(2,0);                    d=mca_metric(p,q);  d2=larc_metric(p,q);  printf("d((%9.6lf,%9.6lf),(%9.6lf,%9.6lf))=%9.6lf %9.6lf t=%6.2lf\n",p.getx(),p.gety(),q.getx(),q.gety(),d,d2,pqtheta(p,q)*180/PI);
  p.put(0,1);                   q.put(1,0);                    d=mca_metric(p,q);  d2=larc_metric(p,q);  printf("d((%9.6lf,%9.6lf),(%9.6lf,%9.6lf))=%9.6lf %9.6lf t=%6.2lf\n",p.getx(),p.gety(),q.getx(),q.gety(),d,d2,pqtheta(p,q)*180/PI);
  p.put(0,1);                   q.put( cos(PI/4), sin(PI/4));  d=mca_metric(p,q);  d2=larc_metric(p,q);  printf("d((%9.6lf,%9.6lf),(%9.6lf,%9.6lf))=%9.6lf %9.6lf t=%6.2lf\n",p.getx(),p.gety(),q.getx(),q.gety(),d,d2,pqtheta(p,q)*180/PI);
  p.put(0,1);                   q.put(-cos(PI/4), sin(PI/4));  d=mca_metric(p,q);  d2=larc_metric(p,q);  printf("d((%9.6lf,%9.6lf),(%9.6lf,%9.6lf))=%9.6lf %9.6lf t=%6.2lf\n",p.getx(),p.gety(),q.getx(),q.gety(),d,d2,pqtheta(p,q)*180/PI);
  p.put(0,1);                   q.put(-cos(PI/4),-sin(PI/4));  d=mca_metric(p,q);  d2=larc_metric(p,q);  printf("d((%9.6lf,%9.6lf),(%9.6lf,%9.6lf))=%9.6lf %9.6lf t=%6.2lf\n",p.getx(),p.gety(),q.getx(),q.gety(),d,d2,pqtheta(p,q)*180/PI);
  p.put(0,1);                   q.put( cos(PI/4),-sin(PI/4));  d=mca_metric(p,q);  d2=larc_metric(p,q);  printf("d((%9.6lf,%9.6lf),(%9.6lf,%9.6lf))=%9.6lf %9.6lf t=%6.2lf\n",p.getx(),p.gety(),q.getx(),q.gety(),d,d2,pqtheta(p,q)*180/PI);
  p.put(1,0);                   q.put(-0.5,-0.5);              d=mca_metric(p,q);  d2=larc_metric(p,q);  printf("d((%9.6lf,%9.6lf),(%9.6lf,%9.6lf))=%9.6lf %9.6lf t=%6.2lf\n",p.getx(),p.gety(),q.getx(),q.gety(),d,d2,pqtheta(p,q)*180/PI);
  p.put(1,0);                   q.put(-2,-2);                  d=mca_metric(p,q);  d2=larc_metric(p,q);  printf("d((%9.6lf,%9.6lf),(%9.6lf,%9.6lf))=%9.6lf %9.6lf t=%6.2lf\n",p.getx(),p.gety(),q.getx(),q.gety(),d,d2,pqtheta(p,q)*180/PI);
  p.put(1,0);                   q.put(0,2);                    d=mca_metric(p,q);  d2=larc_metric(p,q);  printf("d((%9.6lf,%9.6lf),(%9.6lf,%9.6lf))=%9.6lf %9.6lf t=%6.2lf\n",p.getx(),p.gety(),q.getx(),q.gety(),d,d2,pqtheta(p,q)*180/PI);
  p.put( cos(PI/4),-sin(PI/4)); q.put(-2,1);                   d=mca_metric(p,q);  d2=larc_metric(p,q);  printf("d((%9.6lf,%9.6lf),(%9.6lf,%9.6lf))=%9.6lf %9.6lf t=%6.2lf\n",p.getx(),p.gety(),q.getx(),q.gety(),d,d2,pqtheta(p,q)*180/PI);
  p.put( cos(PI/4),-sin(PI/4)); q.put(-1.63,1.33);             d=mca_metric(p,q);  d2=larc_metric(p,q);  printf("d((%9.6lf,%9.6lf),(%9.6lf,%9.6lf))=%9.6lf %9.6lf t=%6.2lf\n",p.getx(),p.gety(),q.getx(),q.gety(),d,d2,pqtheta(p,q)*180/PI);
  p.put( cos(PI/4), sin(PI/4)); q.put( 1, 1);                  d=mca_metric(p,q);  d2=larc_metric(p,q);  printf("d((%9.6lf,%9.6lf),(%9.6lf,%9.6lf))=%9.6lf %9.6lf t=%6.2lf\n",p.getx(),p.gety(),q.getx(),q.gety(),d,d2,pqtheta(p,q)*180/PI);
  p.put(1,0);                   q.put(2,0);                    d=mca_metric(p,q);  d2=larc_metric(p,q);  printf("d((%9.6lf,%9.6lf),(%9.6lf,%9.6lf))=%9.6lf %9.6lf t=%6.2lf\n",p.getx(),p.gety(),q.getx(),q.gety(),d,d2,pqtheta(p,q)*180/PI);
  p.put(1,0);                   q.put(-1,0);                   d=mca_metric(p,q);  d2=larc_metric(p,q);  printf("d((%9.6lf,%9.6lf),(%9.6lf,%9.6lf))=%9.6lf %9.6lf t=%6.2lf\n",p.getx(),p.gety(),q.getx(),q.gety(),d,d2,pqtheta(p,q)*180/PI);
  p.put(1,1);                   q.put(0,0);                    d=mca_metric(p,q);  d2=larc_metric(p,q);  printf("d((%9.6lf,%9.6lf),(%9.6lf,%9.6lf))=%9.6lf %9.6lf t=%6.2lf\n",p.getx(),p.gety(),q.getx(),q.gety(),d,d2,pqtheta(p,q)*180/PI);
  p.put(1,0);                   q.put(0,-5);                   d=mca_metric(p,q);  d2=larc_metric(p,q);  printf("d((%9.6lf,%9.6lf),(%9.6lf,%9.6lf))=%9.6lf %9.6lf t=%6.2lf\n",p.getx(),p.gety(),q.getx(),q.gety(),d,d2,pqtheta(p,q)*180/PI);
  p.put(0,1);                   q.put(1,0);                    d=mca_metric(p,q);  d2=larc_metric(p,q);  printf("d((%9.6lf,%9.6lf),(%9.6lf,%9.6lf))=%9.6lf %9.6lf t=%6.2lf\n",p.getx(),p.gety(),q.getx(),q.gety(),d,d2,pqtheta(p,q)*180/PI);
  p.put(0,1);                   q.put(2,0);                    d=mca_metric(p,q);  d2=larc_metric(p,q);  printf("d((%9.6lf,%9.6lf),(%9.6lf,%9.6lf))=%9.6lf %9.6lf t=%6.2lf\n",p.getx(),p.gety(),q.getx(),q.gety(),d,d2,pqtheta(p,q)*180/PI);
  p.put(0,1);                   q.put(0,-2);                   d=mca_metric(p,q);  d2=larc_metric(p,q);  printf("d((%9.6lf,%9.6lf),(%9.6lf,%9.6lf))=%9.6lf %9.6lf t=%6.2lf\n",p.getx(),p.gety(),q.getx(),q.gety(),d,d2,pqtheta(p,q)*180/PI);
  p.put(0,1);                   q.put(-2,1);                   d=mca_metric(p,q);  d2=larc_metric(p,q);  printf("d((%9.6lf,%9.6lf),(%9.6lf,%9.6lf))=%9.6lf %9.6lf t=%6.2lf\n",p.getx(),p.gety(),q.getx(),q.gety(),d,d2,pqtheta(p,q)*180/PI);
  p.put(1,0);                   q.put(2,0);                    d=mca_metric(p,q);  d2=larc_metric(p,q);  printf("d((%9.6lf,%9.6lf),(%9.6lf,%9.6lf))=%9.6lf %9.6lf t=%6.2lf\n",p.getx(),p.gety(),q.getx(),q.gety(),d,d2,pqtheta(p,q)*180/PI);
  p.put(1.5,1.5);               q.put(1.75,1.75);              d=mca_metric(p,q);  d2=larc_metric(p,q);  printf("d((%9.6lf,%9.6lf),(%9.6lf,%9.6lf))=%9.6lf %9.6lf t=%6.2lf\n",p.getx(),p.gety(),q.getx(),q.gety(),d,d2,pqtheta(p,q)*180/PI);
  p.put(0,0.5);                 q.put(2.75,5);                 d=mca_metric(p,q);  d2=larc_metric(p,q);  printf("d((%9.6lf,%9.6lf),(%9.6lf,%9.6lf))=%9.6lf %9.6lf t=%6.2lf\n",p.getx(),p.gety(),q.getx(),q.gety(),d,d2,pqtheta(p,q)*180/PI);
  p.put(0,0.5);                 q.put(1,100);                  d=mca_metric(p,q);  d2=larc_metric(p,q);  printf("d((%9.6lf,%9.6lf),(%9.6lf,%9.6lf))=%9.6lf %9.6lf t=%6.2lf\n",p.getx(),p.gety(),q.getx(),q.gety(),d,d2,pqtheta(p,q)*180/PI);
  return 1;
  }

/*-------------------------------------------------------------------------
 * test spinner
 *-------------------------------------------------------------------------*/
int test_conj(void){
  complex z,zc;

  z.put(1.0,1.0);   zc=z.conj();  z.list();  zc.list(); putchar('\n');
  z.put(-3.2,4.7);  zc=z.conj();  z.list();  zc.list(); putchar('\n');
  return 0;
  }

/*-------------------------------------------------------------------------
 * test spinner
 *-------------------------------------------------------------------------*/
int test_spinner(void){
  int n,m;
  const long N=160;
  spinseq x(N);
  seqR1 cseq(2*N+1);

  printf("Test spinner routines\n");
  printf("---------------------\n");
  printf("%3.1lf ",spin_metric('A','A'));
  for(m=0;m<6;m++)printf("%3.1lf ",spin_metric('A',(char)('A'+m)));
  putchar('\n');
  for(n=0;n<6;n++){
    printf("%3.1lf ",spin_metric('A','A'+n));
    for(m=0;m<6;m++)printf("%3.1lf ",spin_metric((char)('A'+n),(char)('A'+m)));
    putchar('\n');
    }

  x.randomize(0x5EED);
  //cseq=spin_correlation(x,x,':'); //auto-correlation for spinner sequence

  x.list(          ); putchar('\n');
  cseq.list(0,  1*N); putchar('\n');
  cseq.list(N+1,2*N); putchar('\n');

  return 1;
  }

/*-------------------------------------------------------------------------
 * test die seqR1 functions
 *-------------------------------------------------------------------------*/
int test_die(void){
  const long N=100;
  dieseq x(N);

  printf("Test die seqR1 functions\n");
  printf("---------------------------\n");
  x.randomize(0x5EED);
  x.list();putchar('\n');
  return 0;
  }

/*-------------------------------------------------------------------------
 * test real die seqR1 functions
 *-------------------------------------------------------------------------*/
int test_rdie(void){
  long m;
  const long N=100;
  rdieseq x(N);
  seqR1 Rxx(2*N+1);
  seqR1 histo(8);

  printf("Test real die seqR1 operations\n");
  printf("---------------------------------\n");
  x.randomize(0x5EED);
  x.list();putchar('\n');

  printf("Test real die histogram operation\n");
  printf("---------------------------------\n");
  histo=x.histogram();
  histo.list();putchar('\n');

  printf("Test rdie_metric d(a,b)\n");
  printf("-----------------------\n");
  x.metrictbl();

  printf("Test auto-correlation Rxx\n");
  printf("-------------------------\n");
  m=-N; printf("Rxx(%4ld)=%8.3lf\n",m,x.Rxx(m));
  m=-1; printf("Rxx(%4ld)=%8.3lf\n",m,x.Rxx(m));
  m= 0; printf("Rxx(%4ld)=%8.3lf\n",m,x.Rxx(m));
  m= 1; printf("Rxx(%4ld)=%8.3lf\n",m,x.Rxx(m));
  m= N; printf("Rxx(%4ld)=%8.3lf\n",m,x.Rxx(m));

  if(x.Rxx(&Rxx,1)){fprintf(stderr,"ERROR computing Rxx for seqR1 x.\n"); return -1;} //auto-correlation seqR1 of truncated xR3hl
  // |      |   |____________switch to turn on counting display
  // |      |________________pointer to output correlation sequence
  // |_______________________input real die sequence

  Rxx.list();

  return 0;
  }

/*-------------------------------------------------------------------------
 * test die to C^1 mapping
 *-------------------------------------------------------------------------*/
TEST( TestSuiteDie, C1 )
{
  ASSERT_DOUBLE_EQ( (die_dietoC1c('0')).real(), 0.0 );
  ASSERT_DOUBLE_EQ( (die_dietoC1c('0')).imag(), 0.0 );
  double theta = 30.0;
  ASSERT_DOUBLE_EQ( (die_dietoC1c('A')).real(), cos( theta / 180. * M_PI ) );
  ASSERT_DOUBLE_EQ( (die_dietoC1c('A')).imag(), sin( theta / 180. * M_PI ) );
  theta += 60;
  ASSERT_DOUBLE_EQ( (die_dietoC1c('B')).real(), cos( theta / 180. * M_PI ) );
  ASSERT_DOUBLE_EQ( (die_dietoC1c('B')).imag(), sin( theta / 180. * M_PI ) );
  theta += 60;
  ASSERT_DOUBLE_EQ( (die_dietoC1c('C')).real(), cos( theta / 180. * M_PI ) );
  ASSERT_DOUBLE_EQ( (die_dietoC1c('C')).imag(), sin( theta / 180. * M_PI ) );
  theta += 60;
  ASSERT_FLOAT_EQ(  (die_dietoC1c('D')).real(), cos( theta / 180. * M_PI ) );
  ASSERT_FLOAT_EQ(  (die_dietoC1c('D')).imag(), sin( theta / 180. * M_PI ) );
  theta += 60;
  ASSERT_DOUBLE_EQ( (die_dietoC1c('E')).real(), cos( theta / 180. * M_PI ) );
  ASSERT_DOUBLE_EQ( (die_dietoC1c('E')).imag(), sin( theta / 180. * M_PI ) );
}

/*-------------------------------------------------------------------------
 * test expi function
 *-------------------------------------------------------------------------*/
int test_expi(void){
  const long N=12;
  long n;
  complex y;
  double theta;
  for(n=0; n<=N; n++){
    theta = 2.0*PI*(double)n/(double)N;
    y = expi(theta);
    printf("n=%2ld, theta=%lf (%3.0lf), expi(theta)=(%+9.6lf,%+9.6lf) mag=%lf\n",n,theta,theta/PI*180,y.real(),y.imag(), y.mag());
    }
  return 0;
  }

/*-------------------------------------------------------------------------
 * test DFT operator
 *-------------------------------------------------------------------------*/
int test_dft_R1(void){
  long N=100;
  seqR1 y(N);
  seqC1 Dy(N);
  seqR1 mDy(N);
  complex yn;
  long n;
  printf("\ntest DFT operations with cos(2pi n/10):\n");
  printf("---------------------------------------\n");
  cos((2.0*PI/10.0),&y);
  dft(&y,&Dy);
  mag(&Dy,&mDy);
  printf("\n(y_n)=cos(2 pi n/10):\n");  y.list();
  printf("\nDFT(y_n):\n");              Dy.list();
  printf("\n|DFT(y_n)|:\n");            mDy.list();
  printf("\n");                         mDy.max('p');

  printf("\ntest DFT operations with sin(2pi n/10):\n");
  printf("---------------------------------------\n");
  sin((2.0*PI/10.0),&y);
  dft(&y,&Dy);
  mag(&Dy,&mDy);
  printf("\n(y_n)=sin(2 pi n/10):\n");  y.list();
  printf("\nDFT(y_n):\n");              Dy.list();
  printf("\n|DFT(y_n)|:\n");            mDy.list();
  printf("\n");                         mDy.max('p');

  n=10; yn=dftn(&y,n); printf("y[%02ld]=(%+lf,%+lf) |y[%02ld]|=%lf\n",n,yn.real(),yn.imag(),n,yn.mag());
  n=90; yn=dftn(&y,n); printf("y[%02ld]=(%+lf,%+lf) |y[%02ld]|=%lf\n",n,yn.real(),yn.imag(),n,yn.mag());

  return 0;
  }