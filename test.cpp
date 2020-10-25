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

/*-------------------------------------------------------------------------
 * test opair operations
 *-------------------------------------------------------------------------*/
int test_opair(void){
  opair ab(3,-4);
  printf("\ntest opair operations:\n");
  printf("----------------------\n");
  printf("ab = (%lf,%lf)\n",ab.getx(),ab.gety());
  ab.clear();
  printf("clear-->(%lf,%lf)\n",ab.getx(),ab.gety());
  return 0;
  }

/*-------------------------------------------------------------------------
 * test vectR2 operations
 *-------------------------------------------------------------------------*/
int test_vectR2(void){
  vectR2 p(3,-4);
  vectR2 q,s;

  printf("\ntest vectR2 operations:\n");
  printf("----------------------\n");
                                   printf("p = (%lf,%lf)\n",p.getx(),p.gety());
                                   printf("norm(%lf,%lf)=%lf\n",p.getx(),p.gety(),p.norm());
  p.put(sqrt(3),1);                printf("theta(%lf,%lf)=%lf degrees\n",p.getx(),p.gety(),p.theta()/PI*180);
  q=-p;                            printf("-(%lf,%lf)=(%lf,%lf)\n",p.getx(),p.gety(),q.getx(),q.gety());
  p.put(3,-5); q.put(-2,7); s=p+q; printf("(%lf,%lf)+(%lf,%lf)=(%lf,%lf)\n",p.getx(),p.gety(),q.getx(),q.gety(),s.getx(),s.gety());
  p.put(3,-5); q.put(-2,7); s=p-q; printf("(%lf,%lf)-(%lf,%lf)=(%lf,%lf)\n",p.getx(),p.gety(),q.getx(),q.gety(),s.getx(),s.gety());
  p.put(2,-3); q.put(-5,7);        printf("(%lf,%lf)*(%lf,%lf)=%lf\n",p.getx(),p.gety(),q.getx(),q.gety(),p^q);
  q=p; q&=(PI/2);                  printf("rotate(%lf,%lf)by 90degrees = (%lf,%lf)\n",p.getx(),p.gety(),q.getx(),q.gety());
  q=p; q&=(PI);                    printf("rotate(%lf,%lf)by 180degrees = (%lf,%lf)\n",p.getx(),p.gety(),q.getx(),q.gety());
  p.clear();                       printf("clear operation--> (%lf,%lf)\n",p.getx(),p.gety());
  return 0;
  }

/*-------------------------------------------------------------------------
 * test complex operations
 *-------------------------------------------------------------------------*/
int test_complex(void){
  complex p(3,-4);
  complex q,s;

  printf("\ntest complex operations:\n");
  printf("-----------------------\n");
  printf("p = (%lf,%lf)\n",p.getx(),p.gety());
  printf("norm(%lf,%lf)=%lf\n",p.getx(),p.gety(),p.norm());
  p.put(10*sqrt(3),2);             printf("theta(%lf,%lf)=%lf degrees\n",           p.getx(), p.gety(), p.theta()/PI*180);
  q=-p;                            printf("-(%lf,%lf)=(%lf,%lf)\n",p.getx(),           p.gety(), q.getx(), q.gety());
  p.put(3,-5); q.put(-2,7); s=p+q; printf("(%lf,%lf)+(%lf,%lf)=(%lf,%lf)\n",           p.getx(), p.gety(), q.getx(), q.gety(), s.getx(), s.gety());
  p.put(3,-5); q.put(-2,7); s=p-q; printf("(%lf,%lf)-(%lf,%lf)=(%lf,%lf)\n",           p.getx(), p.gety(), q.getx(), q.gety(), s.getx(), s.gety());
  p.put(2,-3); q.put(-5,7); s=p*q; printf("(%lf,%lf)*(%lf,%lf)=(%lf,%lf)\n",           p.getx(), p.gety(), q.getx(), q.gety(), s.getx(), s.gety());
  q=p; q&=(PI/2);                  printf("rotate(%lf,%lf)by 90degrees = (%lf,%lf)\n", p.getx(), p.gety(), q.getx(), q.gety() );
  q=p; q&=(PI);                    printf("rotate(%lf,%lf)by 180degrees = (%lf,%lf)\n",p.getx(), p.gety(), q.getx(), q.gety() );
  p.clear();                       printf("clear operation--> (%lf,%lf)\n",            p.getx(), p.gety());
  return 0;
  }

/*-------------------------------------------------------------------------
 * test otriple operations
 *-------------------------------------------------------------------------*/
int test_otriple(void){
  otriple abc(3,-4,5);
  printf("\ntest otriple operations:\n");
  printf("---------------------------\n");
  printf("abc = (%lf,%lf,%lf)\n",abc.getx(),abc.gety(),abc.getz());
  abc.clear();
  printf("clear-->(%lf,%lf,%lf)\n",abc.getx(),abc.gety(),abc.getz());
  return 0;
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
 * test osix operations
 *-------------------------------------------------------------------------*/
int test_osix(void){
  osix x(1,-2,3,-5,7,-11);
  osix y,z;
  printf("\ntest osix operations:\n");
  printf("---------------------\n");
  y=x;
  z.put(y);
  //v=z.get();
  printf("x construct:(%lf,%lf,%lf,%lf,%lf,%lf)\n",x.get1(),x.get2(),x.get3(),x.get4(),x.get5(),x.get6());
  printf("y=x:        (%lf,%lf,%lf,%lf,%lf,%lf)\n",y.get1(),y.get2(),y.get3(),y.get4(),y.get5(),y.get6());
  printf("z.put(y):   (%lf,%lf,%lf,%lf,%lf,%lf)\n",z.get1(),z.get2(),z.get3(),z.get4(),z.get5(),z.get6());
  printf("z.put(y):   (%lf,%lf,%lf,%lf,%lf,%lf)\n",z.get(0),z.get(1),z.get(2),z.get(3),z.get(4),z.get(5));
  x.clear();
  //x=(1,2,3,4,5,6);
  printf("clear x:    (%lf,%lf,%lf,%lf,%lf,%lf)\n",x.get1(),x.get2(),x.get3(),x.get4(),x.get5(),x.get6());
  return 0;
  }

/*-------------------------------------------------------------------------
 * test osix operations
 *-------------------------------------------------------------------------*/
int test_vectR6(void){
  vectR6 x(1,-2,3,-5,7,-11);
  vectR6 y,z,v;
  printf("\ntest vectR6 operations:\n");
  printf("-----------------------\n");
             printf("x construct: (%9.6lf,%9.6lf,%9.6lf,%9.6lf,%9.6lf,%9.6lf)\n",x.get1(),x.get2(),x.get3(),x.get4(),x.get5(),x.get6());
  y=x;       printf("y=x:         (%9.6lf,%9.6lf,%9.6lf,%9.6lf,%9.6lf,%9.6lf)\n",y.get1(),y.get2(),y.get3(),y.get4(),y.get5(),y.get6());
  z.put(y);  printf("z.put(y):    (%9.6lf,%9.6lf,%9.6lf,%9.6lf,%9.6lf,%9.6lf)\n",z.get1(),z.get2(),z.get3(),z.get4(),z.get5(),z.get6());
             printf("z.put(y):    (%9.6lf,%9.6lf,%9.6lf,%9.6lf,%9.6lf,%9.6lf)\n",z.get(0),z.get(1),z.get(2),z.get(3),z.get(4),z.get(5));
             printf("z.list():    ");z.list();printf("\n");
  x.clear(); printf("clear x:     ");x.list();printf("\n");
  x.put(1,2,3,4,5,6); 
  y.put(1,1./2.,1./3.,1./4.,1./5.,1./6.);
             x.list("x.put(1..6): ","\n");
             y.list("y.put(1/x..):");
  z=x+y;     z.list("z=x+y:       ");
  v=z-y;     v.list("v=z-y:       ");
  v+=y;      v.list("v+=y:        ");
  v-=y;      v.list("v-=y:        ");
  y=2*v;     y.list("y=2*v:       ");
  v*=-3;     v.list("v*=-3:       ");
  printf("x.mag()=%lf  x.r()=%lf   x.norm()=%lf\n",x.mag(),x.r(),x.norm());
  printf("x dot y = %lf\n",x^y);
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
 * test Lagrange arc metric in R^2
 *-------------------------------------------------------------------------*/
int test_larc_metric_R2(void){
  vectR2 p,q;
  double d,dN;
  int fails=0;
  printf("Lagrange arc metric tests in R2\n");
  printf("--------------------\n");
  p.put(1,0);                   q.put(-1,0);                   d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  printf("d((%6.3lf,%6.3lf),(%6.3lf,%6.3lf))=%13.10lf ~= %13.10lf  ",p.getx(),p.gety(),q.getx(),q.gety(),d,dN);if(dN<d*0.999||dN>d*1.001){printf("FAIL\n");fails++;}else printf("ok\n");
  q.put(1,0);                   p.put(-1,0);                   d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  printf("d((%6.3lf,%6.3lf),(%6.3lf,%6.3lf))=%13.10lf ~= %13.10lf  ",p.getx(),p.gety(),q.getx(),q.gety(),d,dN);if(dN<d*0.999||dN>d*1.001){printf("FAIL\n");fails++;}else printf("ok\n");
  p.put(0,2);                   q.put(2,0);                    d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  printf("d((%6.3lf,%6.3lf),(%6.3lf,%6.3lf))=%13.10lf ~= %13.10lf  ",p.getx(),p.gety(),q.getx(),q.gety(),d,dN);if(dN<d*0.999||dN>d*1.001){printf("FAIL\n");fails++;}else printf("ok\n");
  p.put(0,1);                   q.put(1,0);                    d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  printf("d((%6.3lf,%6.3lf),(%6.3lf,%6.3lf))=%13.10lf ~= %13.10lf  ",p.getx(),p.gety(),q.getx(),q.gety(),d,dN);if(dN<d*0.999||dN>d*1.001){printf("FAIL\n");fails++;}else printf("ok\n");
  p.put(0,1);                   q.put( cos(PI/4), sin(PI/4));  d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  printf("d((%6.3lf,%6.3lf),(%6.3lf,%6.3lf))=%13.10lf ~= %13.10lf  ",p.getx(),p.gety(),q.getx(),q.gety(),d,dN);if(dN<d*0.999||dN>d*1.001){printf("FAIL\n");fails++;}else printf("ok\n");
  p.put(0,1);                   q.put(-cos(PI/4), sin(PI/4));  d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  printf("d((%6.3lf,%6.3lf),(%6.3lf,%6.3lf))=%13.10lf ~= %13.10lf  ",p.getx(),p.gety(),q.getx(),q.gety(),d,dN);if(dN<d*0.999||dN>d*1.001){printf("FAIL\n");fails++;}else printf("ok\n");
  p.put(0,1);                   q.put(-cos(PI/4),-sin(PI/4));  d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  printf("d((%6.3lf,%6.3lf),(%6.3lf,%6.3lf))=%13.10lf ~= %13.10lf  ",p.getx(),p.gety(),q.getx(),q.gety(),d,dN);if(dN<d*0.999||dN>d*1.001){printf("FAIL\n");fails++;}else printf("ok\n");
  p.put(0,1);                   q.put( cos(PI/4),-sin(PI/4));  d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  printf("d((%6.3lf,%6.3lf),(%6.3lf,%6.3lf))=%13.10lf ~= %13.10lf  ",p.getx(),p.gety(),q.getx(),q.gety(),d,dN);if(dN<d*0.999||dN>d*1.001){printf("FAIL\n");fails++;}else printf("ok\n");
  p.put(1,0);                   q.put(-0.5,-0.5);              d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  printf("d((%6.3lf,%6.3lf),(%6.3lf,%6.3lf))=%13.10lf ~= %13.10lf  ",p.getx(),p.gety(),q.getx(),q.gety(),d,dN);if(dN<d*0.999||dN>d*1.001){printf("FAIL\n");fails++;}else printf("ok\n");
  p.put(1,0);                   q.put(-2,-2);                  d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  printf("d((%6.3lf,%6.3lf),(%6.3lf,%6.3lf))=%13.10lf ~= %13.10lf  ",p.getx(),p.gety(),q.getx(),q.gety(),d,dN);if(dN<d*0.999||dN>d*1.001){printf("FAIL\n");fails++;}else printf("ok\n");
  p.put(1,0);                   q.put(0,2);                    d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  printf("d((%6.3lf,%6.3lf),(%6.3lf,%6.3lf))=%13.10lf ~= %13.10lf  ",p.getx(),p.gety(),q.getx(),q.gety(),d,dN);if(dN<d*0.999||dN>d*1.001){printf("FAIL\n");fails++;}else printf("ok\n");
  p.put( cos(PI/4),-sin(PI/4)); q.put(-2,1);                   d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  printf("d((%6.3lf,%6.3lf),(%6.3lf,%6.3lf))=%13.10lf ~= %13.10lf  ",p.getx(),p.gety(),q.getx(),q.gety(),d,dN);if(dN<d*0.999||dN>d*1.001){printf("FAIL\n");fails++;}else printf("ok\n");
  p.put( cos(PI/4),-sin(PI/4)); q.put(-1.63,1.33);             d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  printf("d((%6.3lf,%6.3lf),(%6.3lf,%6.3lf))=%13.10lf ~= %13.10lf  ",p.getx(),p.gety(),q.getx(),q.gety(),d,dN);if(dN<d*0.999||dN>d*1.001){printf("FAIL\n");fails++;}else printf("ok\n");
  p.put( cos(PI/4), sin(PI/4)); q.put( 1, 1);                  d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  printf("d((%6.3lf,%6.3lf),(%6.3lf,%6.3lf))=%13.10lf ~= %13.10lf  ",p.getx(),p.gety(),q.getx(),q.gety(),d,dN);if(dN<d*0.999||dN>d*1.001){printf("FAIL\n");fails++;}else printf("ok\n");
  p.put(1,0);                   q.put(2,0);                    d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  printf("d((%6.3lf,%6.3lf),(%6.3lf,%6.3lf))=%13.10lf ~= %13.10lf  ",p.getx(),p.gety(),q.getx(),q.gety(),d,dN);if(dN<d*0.999||dN>d*1.001){printf("FAIL\n");fails++;}else printf("ok\n");
  p.put(1,0);                   q.put(-1,0);                   d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  printf("d((%6.3lf,%6.3lf),(%6.3lf,%6.3lf))=%13.10lf ~= %13.10lf  ",p.getx(),p.gety(),q.getx(),q.gety(),d,dN);if(dN<d*0.999||dN>d*1.001){printf("FAIL\n");fails++;}else printf("ok\n");
  p.put(1,1);                   q.put(0,0);                    d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  printf("d((%6.3lf,%6.3lf),(%6.3lf,%6.3lf))=%13.10lf ~= %13.10lf  ",p.getx(),p.gety(),q.getx(),q.gety(),d,dN);if(dN<d*0.999||dN>d*1.001){printf("FAIL\n");fails++;}else printf("ok\n");
  p.put(1,0);                   q.put(0,-1);                   d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  printf("d((%6.3lf,%6.3lf),(%6.3lf,%6.3lf))=%13.10lf ~= %13.10lf  ",p.getx(),p.gety(),q.getx(),q.gety(),d,dN);if(dN<d*0.999||dN>d*1.001){printf("FAIL\n");fails++;}else printf("ok\n");
  p.put(0,1);                   q.put(-cos(PI/4), sin(PI/4));  d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  printf("d((%6.3lf,%6.3lf),(%6.3lf,%6.3lf))=%13.10lf ~= %13.10lf  ",p.getx(),p.gety(),q.getx(),q.gety(),d,dN);if(dN<d*0.999||dN>d*1.001){printf("FAIL\n");fails++;}else printf("ok\n");
  p.put(0,1);                   q.put( cos(PI/4), sin(PI/4));  d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  printf("d((%6.3lf,%6.3lf),(%6.3lf,%6.3lf))=%13.10lf ~= %13.10lf  ",p.getx(),p.gety(),q.getx(),q.gety(),d,dN);if(dN<d*0.999||dN>d*1.001){printf("FAIL\n");fails++;}else printf("ok\n");
  p.put(1,0);                   q.put(-2,0);                   d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  printf("d((%6.3lf,%6.3lf),(%6.3lf,%6.3lf))=%13.10lf ~= %13.10lf  ",p.getx(),p.gety(),q.getx(),q.gety(),d,dN);if(dN<d*0.999||dN>d*1.001){printf("FAIL\n");fails++;}else printf("ok\n");
  q.put(1,0);                   p.put(-2,0);                   d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  printf("d((%6.3lf,%6.3lf),(%6.3lf,%6.3lf))=%13.10lf ~= %13.10lf  ",p.getx(),p.gety(),q.getx(),q.gety(),d,dN);if(dN<d*0.999||dN>d*1.001){printf("FAIL\n");fails++;}else printf("ok\n");
  p.put(1,0);                   q.put(-0.5,0);                 d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  printf("d((%6.3lf,%6.3lf),(%6.3lf,%6.3lf))=%13.10lf ~= %13.10lf  ",p.getx(),p.gety(),q.getx(),q.gety(),d,dN);if(dN<d*0.999||dN>d*1.001){printf("FAIL\n");fails++;}else printf("ok\n");
  printf("Theorem 3.9 (1): triangle inequality fails\n");
  p.put(1,0);                   q.put(-0.5,0);                 d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  printf("d((%6.3lf,%6.3lf),(%6.3lf,%6.3lf))=%13.10lf ~= %13.10lf  ",p.getx(),p.gety(),q.getx(),q.gety(),d,dN);if(dN<d*0.999||dN>d*1.001){printf("FAIL\n");fails++;}else printf("ok\n");
  p.put(1,0);                   q.put(-0.5,0.2);               d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  printf("d((%6.3lf,%6.3lf),(%6.3lf,%6.3lf))=%13.10lf ~= %13.10lf  ",p.getx(),p.gety(),q.getx(),q.gety(),d,dN);if(dN<d*0.999||dN>d*1.001){printf("FAIL\n");fails++;}else printf("ok\n");
  p.put(-0.5,0.2);              q.put(-0.5,0);                 d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  printf("d((%6.3lf,%6.3lf),(%6.3lf,%6.3lf))=%13.10lf ~= %13.10lf  ",p.getx(),p.gety(),q.getx(),q.gety(),d,dN);if(dN<d*0.999||dN>d*1.001){printf("FAIL\n");fails++;}else printf("ok\n");
  printf("Theorem 3.9 (2): translation invariance fails\n");
  p.put(1,0.5);                 q.put(0.5,1);                  d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  printf("d((%6.3lf,%6.3lf),(%6.3lf,%6.3lf))=%13.10lf ~= %13.10lf  ",p.getx(),p.gety(),q.getx(),q.gety(),d,dN);if(dN<d*0.999||dN>d*1.001){printf("FAIL\n");fails++;}else printf("ok\n");
  p.put(0.5,0);                 q.put(0,0.5);                  d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  printf("d((%6.3lf,%6.3lf),(%6.3lf,%6.3lf))=%13.10lf ~= %13.10lf  ",p.getx(),p.gety(),q.getx(),q.gety(),d,dN);if(dN<d*0.999||dN>d*1.001){printf("FAIL\n");fails++;}else printf("ok\n");
  printf("Theorem 3.9 (5): balls are not convex\n");
  p.put(0,1);                   q.put(-0.70,-1.12);            d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  printf("d((%6.3lf,%6.3lf),(%6.3lf,%6.3lf))=%13.10lf ~= %13.10lf  ",p.getx(),p.gety(),q.getx(),q.gety(),d,dN);if(dN<d*0.999||dN>d*1.001){printf("FAIL\n");fails++;}else printf("ok\n");
  q.put(0,1);                   p.put(-0.70,-1.12);            d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  printf("d((%6.3lf,%6.3lf),(%6.3lf,%6.3lf))=%13.10lf ~= %13.10lf  ",p.getx(),p.gety(),q.getx(),q.gety(),d,dN);if(dN<d*0.999||dN>d*1.001){printf("FAIL\n");fails++;}else printf("ok\n");
  p.put(0,1);                   q.put( 0.70,-1.12);            d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  printf("d((%6.3lf,%6.3lf),(%6.3lf,%6.3lf))=%13.10lf ~= %13.10lf  ",p.getx(),p.gety(),q.getx(),q.gety(),d,dN);if(dN<d*0.999||dN>d*1.001){printf("FAIL\n");fails++;}else printf("ok\n");
  q.put(0,1);                   p.put( 0.70,-1.12);            d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  printf("d((%6.3lf,%6.3lf),(%6.3lf,%6.3lf))=%13.10lf ~= %13.10lf  ",p.getx(),p.gety(),q.getx(),q.gety(),d,dN);if(dN<d*0.999||dN>d*1.001){printf("FAIL\n");fails++;}else printf("ok\n");
  p.put(0,1);                   q.put( 0,   -1.12);            d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  printf("d((%6.3lf,%6.3lf),(%6.3lf,%6.3lf))=%13.10lf ~= %13.10lf  ",p.getx(),p.gety(),q.getx(),q.gety(),d,dN);if(dN<d*0.999||dN>d*1.001){printf("FAIL\n");fails++;}else printf("ok\n");
  q.put(0,1);                   p.put( 0,   -1.12);            d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  printf("d((%6.3lf,%6.3lf),(%6.3lf,%6.3lf))=%13.10lf ~= %13.10lf  ",p.getx(),p.gety(),q.getx(),q.gety(),d,dN);if(dN<d*0.999||dN>d*1.001){printf("FAIL\n");fails++;}else printf("ok\n");
  printf("Theorem 3.10: Lagrange arc distance versus Euclidean metric\n");
  p.put(1,0);                   q.put(-0.50,0);                d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  printf("d((%6.3lf,%6.3lf),(%6.3lf,%6.3lf))=%13.10lf ~= %13.10lf  ",p.getx(),p.gety(),q.getx(),q.gety(),d,dN);if(dN<d*0.999||dN>d*1.001){printf("FAIL\n");fails++;}else printf("ok\n");
  q.put(1,0);                   p.put(-0.50,0);                d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  printf("d((%6.3lf,%6.3lf),(%6.3lf,%6.3lf))=%13.10lf ~= %13.10lf  ",p.getx(),p.gety(),q.getx(),q.gety(),d,dN);if(dN<d*0.999||dN>d*1.001){printf("FAIL\n");fails++;}else printf("ok\n");
  p.put(1,0);                   q.put(-0.50,0.75);             d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  printf("d((%6.3lf,%6.3lf),(%6.3lf,%6.3lf))=%13.10lf ~= %13.10lf  ",p.getx(),p.gety(),q.getx(),q.gety(),d,dN);if(dN<d*0.999||dN>d*1.001){printf("FAIL\n");fails++;}else printf("ok\n");
  q.put(1,0);                   p.put(-0.50,0.75);             d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  printf("d((%6.3lf,%6.3lf),(%6.3lf,%6.3lf))=%13.10lf ~= %13.10lf  ",p.getx(),p.gety(),q.getx(),q.gety(),d,dN);if(dN<d*0.999||dN>d*1.001){printf("FAIL\n");fails++;}else printf("ok\n");
  printf("Example 3.12: Lagrange arc distance in R^2\n");
  p.put( 0,1);                  q.put( 1, 0);                  d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  printf("d((%6.3lf,%6.3lf),(%6.3lf,%6.3lf))=%13.10lf ~= %13.10lf  ",p.getx(),p.gety(),q.getx(),q.gety(),d,dN);if(dN<d*0.999||dN>d*1.001){printf("FAIL\n");fails++;}else printf("ok\n");
  p.put( 0,1);                  q.put(-1, 0);                  d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  printf("d((%6.3lf,%6.3lf),(%6.3lf,%6.3lf))=%13.10lf ~= %13.10lf  ",p.getx(),p.gety(),q.getx(),q.gety(),d,dN);if(dN<d*0.999||dN>d*1.001){printf("FAIL\n");fails++;}else printf("ok\n");
  p.put( 0,1);                  q.put( 0,-1);                  d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  printf("d((%6.3lf,%6.3lf),(%6.3lf,%6.3lf))=%13.10lf ~= %13.10lf  ",p.getx(),p.gety(),q.getx(),q.gety(),d,dN);if(dN<d*0.999||dN>d*1.001){printf("FAIL\n");fails++;}else printf("ok\n");
  p.put( 1,0);                  q.put( 0,-1);                  d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  printf("d((%6.3lf,%6.3lf),(%6.3lf,%6.3lf))=%13.10lf ~= %13.10lf  ",p.getx(),p.gety(),q.getx(),q.gety(),d,dN);if(dN<d*0.999||dN>d*1.001){printf("FAIL\n");fails++;}else printf("ok\n");
  p.put( 1,0);                  q.put(-1, 0);                  d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  printf("d((%6.3lf,%6.3lf),(%6.3lf,%6.3lf))=%13.10lf ~= %13.10lf  ",p.getx(),p.gety(),q.getx(),q.gety(),d,dN);if(dN<d*0.999||dN>d*1.001){printf("FAIL\n");fails++;}else printf("ok\n");
  p.put(-1,0);                  q.put( 0,-1);                  d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  printf("d((%6.3lf,%6.3lf),(%6.3lf,%6.3lf))=%13.10lf ~= %13.10lf  ",p.getx(),p.gety(),q.getx(),q.gety(),d,dN);if(dN<d*0.999||dN>d*1.001){printf("FAIL\n");fails++;}else printf("ok\n");
  p.put( 0,1);                  q.put( 2, 0);                  d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  printf("d((%6.3lf,%6.3lf),(%6.3lf,%6.3lf))=%13.10lf ~= %13.10lf  ",p.getx(),p.gety(),q.getx(),q.gety(),d,dN);if(dN<d*0.999||dN>d*1.001){printf("FAIL\n");fails++;}else printf("ok\n");
  p.put( 0,1);                  q.put( 0,-2);                  d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  printf("d((%6.3lf,%6.3lf),(%6.3lf,%6.3lf))=%13.10lf ~= %13.10lf  ",p.getx(),p.gety(),q.getx(),q.gety(),d,dN);if(dN<d*0.999||dN>d*1.001){printf("FAIL\n");fails++;}else printf("ok\n");
  p.put( 0,1);                  q.put(-2, 1);                  d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  printf("d((%6.3lf,%6.3lf),(%6.3lf,%6.3lf))=%13.10lf ~= %13.10lf  ",p.getx(),p.gety(),q.getx(),q.gety(),d,dN);if(dN<d*0.999||dN>d*1.001){printf("FAIL\n");fails++;}else printf("ok\n");
  p.put( 0.000001,0);              q.put(PI, 0);                  d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  printf("d((%6.3lf,%6.3lf),(%6.3lf,%6.3lf))=%13.10lf ~= %13.10lf  ",p.getx(),p.gety(),q.getx(),q.gety(),d,dN);if(dN<d*0.999||dN>d*1.001){printf("FAIL\n");fails++;}else printf("ok\n");
  printf("number of fails = %d\n",fails);
  return fails;
  }

/*-------------------------------------------------------------------------
 * test Lagrange arc metric in R^3
 *-------------------------------------------------------------------------*/
int test_larc_metric_R3(void){
  vectR3 p,q;
  double d,de,s=1.0/PI;
  printf("Lagrange arc metric tests in R3\n");
  printf("--------------------\n");
  p.put(0,2,0);                   q.put(2,0,0);                    d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(0,1,0);                   q.put(1,0,0);                    d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(0,1,0);                   q.put( cos(PI/4), sin(PI/4),0);  d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(0,1,0);                   q.put(-cos(PI/4), sin(PI/4),0);  d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(0,1,0);                   q.put(-cos(PI/4),-sin(PI/4),0);  d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(0,1,0);                   q.put( cos(PI/4),-sin(PI/4),0);  d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(1,0,0);                   q.put(-0.5,-0.5,0);              d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(1,0,0);                   q.put(-2,-2,0);                  d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(1,0,0);                   q.put(0,2,0);                    d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( cos(PI/4),-sin(PI/4),0); q.put(-2,1,0);                   d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( cos(PI/4),-sin(PI/4),0); q.put(-1.63,1.33,0);             d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( cos(PI/4), sin(PI/4),0); q.put( 1, 1,0);                  d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(1,0,0);                   q.put(2,0,0);                    d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(1,0,0);                   q.put(-1,0,0);                   d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(1,1,0);                   q.put(0,0,0);                    d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(1,0,0);                   q.put(0,-1,0);                   d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  printf("\n");
  p.put( 0, 1, 0);                q.put( 0, 1, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 1, 0);                q.put( 0, 0, 1);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 1, 0);                q.put(-1, 0, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 1, 0);                q.put( 0,-1, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 1, 0);                q.put( 0, 0,-1);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 2, 0);                q.put( 1, 0, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 2, 0);                q.put( 0, 1, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 2, 0);                q.put( 0, 0, 1);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 2, 0);                q.put(-1, 0, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 2, 0);                q.put( 0,-1, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 2, 0);                q.put( 0, 0,-1);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 1, 0);                q.put( 1, 1, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 1, 0);                q.put( 0, 1, 1);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 1, 0);                q.put( 1, 0, 1);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 1, 0);                q.put(-1,-1, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 1, 0);                q.put( 0,-1,-1);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 1, 0);                q.put(-1, 0,-1);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1, 1, 1);                q.put(-0.5, 0.25,-2);            d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1,-1,-2);                q.put( 1, 0, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1,-1,-2);                q.put( 0, 1, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1,-1,-2);                q.put( 0, 0, 1);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1,-1,-2);                q.put( 0, 0,-1);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1,-1,-2);                q.put( 0,-1, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1,-1,-2);                q.put(-1, 0, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 0, 0);                q.put( 1, 0, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 0, 0);                q.put( 0, 1, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 0, 0);                q.put( 0, 0, 1);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 0, 0);                q.put( 0, 0,-1);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 0, 0);                q.put( 0,-1, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 0, 0);                q.put(-1, 0, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1, 1,-2);                q.put( 1, 0, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1, 1,-2);                q.put( 0, 1, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1, 1,-2);                q.put( 0, 0, 1);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1, 1,-2);                q.put( 0, 0,-1);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1, 1,-2);                q.put( 0,-1, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1, 1,-2);                q.put(-1, 0, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1, 1, 0);                q.put( 1, 0, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1, 1, 0);                q.put( 0, 1, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1, 1, 0);                q.put( 0, 0, 1);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1, 1, 0);                q.put( 0, 0,-1);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1, 1, 0);                q.put( 0,-1, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1, 1, 0);                q.put(-1, 0, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 2, 2, 0);                q.put( 1, 0, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 2, 2, 0);                q.put( 0, 1, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 2, 2, 0);                q.put( 0, 0, 1);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 2, 2, 0);                q.put( 0, 0,-1);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 2, 2, 0);                q.put( 0,-1, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 2, 2, 0);                q.put(-1, 0, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1, 2, 1);                q.put( 1, 0, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1, 2, 1);                q.put( 0, 1, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1, 2, 1);                q.put( 0, 0, 1);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1, 2, 1);                q.put( 0, 0,-1);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1, 2, 1);                q.put( 0,-1, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1, 2, 1);                q.put(-1, 0, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 3, 1);                q.put( 1, 0, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 3, 1);                q.put( 0, 1, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 3, 1);                q.put( 0, 0, 1);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 3, 1);                q.put( 0, 0,-1);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 3, 1);                q.put( 0,-1, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 3, 1);                q.put(-1, 0, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 2, 0);                q.put( 1, 0, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 2, 0);                q.put( 0, 1, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 2, 0);                q.put( 0, 0, 1);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 2, 0);                q.put( 0, 0,-1);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 2, 0);                q.put( 0,-1, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 2, 0);                q.put(-1, 0, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(-1, 2, 1);                q.put( 1, 0, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(-1, 2, 1);                q.put( 0, 1, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(-1, 2, 1);                q.put( 0, 0, 1);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(-1, 2, 1);                q.put( 0, 0,-1);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(-1, 2, 1);                q.put( 0,-1, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(-1, 2, 1);                q.put(-1, 0, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);

  p.put(-1, 0, 0);                q.put( 1, 0, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(-1.01, 0, 0);                q.put(0.99, 0, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(-1.001, 0, 0);                q.put(0.999, 0, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(-1.0001, 0, 0);                q.put(0.9999, 0, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);

  q.put(-1, 0, 0);                p.put( 1, 0, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  q.put(-1.01, 0, 0);             p.put(0.99, 0, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  q.put(-1.001, 0, 0);            p.put(0.999, 0, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  q.put(-1.0001, 0, 0);           p.put(0.9999, 0, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);

  p.put(0,1,0);                   q.put( cos(PI/4), sin(PI/4),0);  d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(0,1,0);                   q.put(-cos(PI/4), sin(PI/4),0);  d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(0,1,0);                   q.put( cos(PI/6), sin(PI/6),0);  d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(0,1,0);                   q.put(-cos(PI/6), sin(PI/6),0);  d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( cos(PI/4), sin(PI/4),0);   q.put( 0.5*cos(PI/8), 0.5*sin(PI/8),0);  d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(-cos(PI/4), sin(PI/4),0);   q.put(-0.5*cos(3*PI/8), 0.5*sin(3*PI/8),0);  d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0.5*cos(PI/4), 0.5*sin(PI/4),0);   q.put( cos(PI/8),   sin(PI/8),0);  d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(-0.5*cos(PI/4), 0.5*sin(PI/4),0);   q.put(-cos(3*PI/8), sin(3*PI/8),0);  d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  printf("Example 3.13: Lagrange arc distance in R^3\n");
  p.put( 0, 1, 0);                q.put( 2, 0, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 1, 0);                q.put( 0,-2, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 1, 0);                q.put(-2, 1, 0);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 1, 0);                q.put(-1, 0,-1);                 d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1, 1, 1);                q.put(-0.5,0.25,-2);             d=larc_metric(p,q); de=ae_metric(s,p,q); printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);

  return 0;
  }

/*-------------------------------------------------------------------------
 * test Lagrange arc metric in R^6
 *-------------------------------------------------------------------------*/
int test_larc_metric_R6(void){
  vectR6 p,q;
  double d;
  printf("Lagrange arc metric tests in R^6\n");
  printf("--------------------------------\n");
  p.put( 1, 0, 0, 0, 0, 0);  q.put( 0, 0, 0, 0, 0, 1); d=larc_metric(p,q);  printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put( 0, 1, 0, 0, 0, 0);  q.put( 0, 0, 0, 0, 1, 0); d=larc_metric(p,q);  printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put( 0, 0, 1, 0, 0, 0);  q.put( 0, 0, 0, 1, 0, 0); d=larc_metric(p,q);  printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put( 0, 0, 0, 1, 0, 0);  q.put( 0, 0, 1, 0, 0, 0); d=larc_metric(p,q);  printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put( 0, 0, 0, 0, 1, 0);  q.put( 0, 1, 0, 0, 0, 0); d=larc_metric(p,q);  printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put( 0, 0, 0, 0, 0, 1);  q.put( 1, 0, 0, 0, 0, 0); d=larc_metric(p,q);  printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  putchar('\n');
  p.put( 1, 0, 0, 0, 0, 0);  q.put(-1, 0, 0, 0, 0, 0); d=larc_metric(p,q);  printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put( 0, 1, 0, 0, 0, 0);  q.put( 0,-1, 0, 0, 0, 0); d=larc_metric(p,q);  printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put( 0, 0, 1, 0, 0, 0);  q.put( 0, 0,-1, 0, 0, 0); d=larc_metric(p,q);  printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put( 0, 0, 0, 1, 0, 0);  q.put( 0, 0, 0,-1, 0, 0); d=larc_metric(p,q);  printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put( 0, 0, 0, 0, 1, 0);  q.put( 0, 0, 0, 0,-1, 0); d=larc_metric(p,q);  printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put( 0, 0, 0, 0, 0, 1);  q.put( 0, 0, 0, 0, 0,-1); d=larc_metric(p,q);  printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  putchar('\n');
  p.put( 1, 0, 0, 0, 0, 0);  q.put( 0, 0, 0, 0, 0, 1); d=larc_metric(p,q);  printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put( 1, 0, 0, 0, 0, 0);  q.put( 0, 0, 0, 0, 1, 0); d=larc_metric(p,q);  printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put( 1, 0, 0, 0, 0, 0);  q.put( 0, 0, 0, 1, 0, 0); d=larc_metric(p,q);  printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put( 1, 0, 0, 0, 0, 0);  q.put( 0, 0, 1, 0, 0, 0); d=larc_metric(p,q);  printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put( 1, 0, 0, 0, 0, 0);  q.put( 0, 1, 0, 0, 0, 0); d=larc_metric(p,q);  printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put( 1, 0, 0, 0, 0, 0);  q.put( 1, 0, 0, 0, 0, 0); d=larc_metric(p,q);  printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  putchar('\n');
  p.put(-1, 2, 1,-1, 2, 1);  q.put( 1, 0, 0, 1, 0, 0); d=larc_metric(p,q);  printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put(-1, 2, 1,-1, 2, 1);  q.put( 0, 1, 0, 0, 1, 0); d=larc_metric(p,q);  printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put(-1, 2, 1,-1, 2, 1);  q.put( 0, 0, 1, 0, 0, 1); d=larc_metric(p,q);  printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put(-1, 2, 1,-1, 2, 1);  q.put( 0, 0,-1, 0, 0,-1); d=larc_metric(p,q);  printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put(-1, 2, 1,-1, 2, 1);  q.put( 0,-1, 0, 0,-1, 0); d=larc_metric(p,q);  printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put(-1, 2, 1,-1, 2, 1);  q.put(-1, 0, 0,-1, 0, 0); d=larc_metric(p,q);  printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  putchar('\n');
  p.put(0.458,0.022,0.120,0.234,0.167,0.000);  q.put(1, 0, 0, 0, 0, 0); d=larc_metric(p,q);  printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put(0.458,0.022,0.120,0.234,0.167,0.000);  q.put(0, 1, 0, 0, 0, 0); d=larc_metric(p,q);  printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put(0.458,0.022,0.120,0.234,0.167,0.000);  q.put(0, 0, 1, 0, 0, 0); d=larc_metric(p,q);  printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put(0.458,0.022,0.120,0.234,0.167,0.000);  q.put(0, 0, 0, 1, 0, 0); d=larc_metric(p,q);  printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put(0.458,0.022,0.120,0.234,0.167,0.000);  q.put(0, 0, 0, 0, 1, 0); d=larc_metric(p,q);  printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put(0.458,0.022,0.120,0.234,0.167,0.000);  q.put(0, 0, 0, 0, 0, 1); d=larc_metric(p,q);  printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  return 0;
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
  long N=160;
  spinseq x(N);
  seqR1 cseq(2*N+1);

  printf("Test spinner routines\n");
  printf("---------------------\n");
  printf("%3.1lf ",spin_metric('0','0'));
  for(m=0;m<6;m++)printf("%3.1lf ",spin_metric('0',(char)('A'+m)));
  putchar('\n');
  for(n=0;n<6;n++){
    printf("%3.1lf ",spin_metric('0','A'+n));
    for(m=0;m<6;m++)printf("%3.1lf ",spin_metric((char)('A'+n),(char)('A'+m)));
    putchar('\n');
    }

  x.randomize(0x5EED);
  //cseq=spin_correlation(x,x,':'); //auto-correlation for spinner sequence

  x.list();putchar('\n');
  cseq.list(0,N);putchar('\n');
  cseq.list(N+1,2*N);putchar('\n');

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
int test_dieC1(void){
  //fdie_dietoC1('A').list();
  printf("\nTest die to C^1:\n");
  printf(  "---------------\n");
  die_dietoC1c('A').list("die A --> ","\n");
  die_dietoC1c('B').list("die B --> ","\n");
  die_dietoC1c('C').list("die C --> ","\n");
  die_dietoC1c('D').list("die D --> ","\n");
  die_dietoC1c('E').list("die E --> ","\n");
  die_dietoC1c('F').list("die F --> ","\n");
  return 0;
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