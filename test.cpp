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

//-----------------------------------------------------------------------------
//! \brief test complex operations
//-----------------------------------------------------------------------------
TEST( TestSuiteSeqR2, max )
{
  seqR2 x(6);

  x.fill(2);
  ASSERT_EQ( x.getx(0),  2 );  ASSERT_EQ( x.gety(0), 2 );
  ASSERT_EQ( x.getx(1),  2 );  ASSERT_EQ( x.gety(1), 2 );

  x.inc(2,-5,1,-5);
  ASSERT_EQ( x.getx(0),  2 );  ASSERT_EQ( x.gety(0),  -5 );
  ASSERT_EQ( x.getx(1),  3 );  ASSERT_EQ( x.gety(1), -10 );
  ASSERT_EQ( x.getx(2),  4 );  ASSERT_EQ( x.gety(2), -15 );
  ASSERT_EQ( x.getx(3),  5 );  ASSERT_EQ( x.gety(3), -20 );
  ASSERT_EQ( x.getx(4),  6 );  ASSERT_EQ( x.gety(4), -25 );
  ASSERT_EQ( x.getx(5),  7 );  ASSERT_EQ( x.gety(5), -30 );

  const vectR2 maxVal = x.max();
  ASSERT_EQ( maxVal.getx(),   7 );  
  ASSERT_EQ( maxVal.gety(), -30 );

  x.clear();
  for( int n=0; n<6; n++ )
  {
    ASSERT_EQ( x.getx(n), 0 );  
    ASSERT_EQ( x.gety(n), 0 );
  }
}

//-----------------------------------------------------------------------------
//! \brief Test Lagrange arc metric in R^2
//-----------------------------------------------------------------------------
TEST( TestSuiteVectRn, theta )
{
  vectR2 p,q;
  double tr,td,dc;
  const double err = 1e-8;
  p.put(0,0);                       q.put(2,0);                    tr=pqtheta(p,q); td=tr*180/PI; dc= -1*180/PI;  ASSERT_NEAR( td, dc, err * fabs(dc) );  //printf("p=(%6.3lf,%6.3lf) q=(%6.3lf,%6.3lf) theta=%6.2lf radians %6.2lf degrees ",p.getx(),p.gety(),q.getx(),q.gety(),tr,td);
  p.put(0,1);                       q.put(0,0);                    tr=pqtheta(p,q); td=tr*180/PI; dc= -2*180/PI;  ASSERT_NEAR( td, dc, err * fabs(dc) );  //printf("p=(%6.3lf,%6.3lf) q=(%6.3lf,%6.3lf) theta=%6.2lf radians %6.2lf degrees ",p.getx(),p.gety(),q.getx(),q.gety(),tr,td);
  p.put(0,2);                       q.put(2,0);                    tr=pqtheta(p,q); td=tr*180/PI; dc= 90;         ASSERT_NEAR( td, dc, err * fabs(dc) );  //printf("p=(%6.3lf,%6.3lf) q=(%6.3lf,%6.3lf) theta=%6.2lf radians %6.2lf degrees ",p.getx(),p.gety(),q.getx(),q.gety(),tr,td);
  p.put(0,1);                       q.put(1,0);                    tr=pqtheta(p,q); td=tr*180/PI; dc= 90;         ASSERT_NEAR( td, dc, err * fabs(dc) );  //printf("p=(%6.3lf,%6.3lf) q=(%6.3lf,%6.3lf) theta=%6.2lf radians %6.2lf degrees ",p.getx(),p.gety(),q.getx(),q.gety(),tr,td);
  p.put(1,3);                       q.put(-1,-3);                  tr=pqtheta(p,q); td=tr*180/PI; dc=180;         ASSERT_NEAR( td, dc, err * fabs(dc) );  //printf("p=(%6.3lf,%6.3lf) q=(%6.3lf,%6.3lf) theta=%6.2lf radians %6.2lf degrees ",p.getx(),p.gety(),q.getx(),q.gety(),tr,td);
  p.put(-2,5);                      q.put(2,-5);                   tr=pqtheta(p,q); td=tr*180/PI; dc=180;         ASSERT_NEAR( td, dc, err * fabs(dc) );  //printf("p=(%6.3lf,%6.3lf) q=(%6.3lf,%6.3lf) theta=%6.2lf radians %6.2lf degrees ",p.getx(),p.gety(),q.getx(),q.gety(),tr,td);
  p.put(0,1);                       q.put( cos(PI/4), sin(PI/4));  tr=pqtheta(p,q); td=tr*180/PI; dc= 45;         ASSERT_NEAR( td, dc, err * fabs(dc) );  //printf("p=(%6.3lf,%6.3lf) q=(%6.3lf,%6.3lf) theta=%6.2lf radians %6.2lf degrees ",p.getx(),p.gety(),q.getx(),q.gety(),tr,td);
  p.put(0,1);                       q.put( cos(PI/4),-sin(PI/4));  tr=pqtheta(p,q); td=tr*180/PI; dc=135;         ASSERT_NEAR( td, dc, err * fabs(dc) );  //printf("p=(%6.3lf,%6.3lf) q=(%6.3lf,%6.3lf) theta=%6.2lf radians %6.2lf degrees ",p.getx(),p.gety(),q.getx(),q.gety(),tr,td);
  p.put(0,1);                       q.put(-cos(PI/4), sin(PI/4));  tr=pqtheta(p,q); td=tr*180/PI; dc= 45;         ASSERT_NEAR( td, dc, err * fabs(dc) );  //printf("p=(%6.3lf,%6.3lf) q=(%6.3lf,%6.3lf) theta=%6.2lf radians %6.2lf degrees ",p.getx(),p.gety(),q.getx(),q.gety(),tr,td);
  p.put(0,1);                       q.put(-cos(PI/4),-sin(PI/4));  tr=pqtheta(p,q); td=tr*180/PI; dc=135;         ASSERT_NEAR( td, dc, err * fabs(dc) );  //printf("p=(%6.3lf,%6.3lf) q=(%6.3lf,%6.3lf) theta=%6.2lf radians %6.2lf degrees ",p.getx(),p.gety(),q.getx(),q.gety(),tr,td);
  p.put( 7*cos(PI/4), 7*sin(PI/4)); q.put(0,1);                    tr=pqtheta(p,q); td=tr*180/PI; dc= 45;         ASSERT_NEAR( td, dc, err * fabs(dc) );  //printf("p=(%6.3lf,%6.3lf) q=(%6.3lf,%6.3lf) theta=%6.2lf radians %6.2lf degrees ",p.getx(),p.gety(),q.getx(),q.gety(),tr,td);
  p.put( 7*cos(PI/4),-7*sin(PI/4)); q.put(0,1);                    tr=pqtheta(p,q); td=tr*180/PI; dc=135;         ASSERT_NEAR( td, dc, err * fabs(dc) );  //printf("p=(%6.3lf,%6.3lf) q=(%6.3lf,%6.3lf) theta=%6.2lf radians %6.2lf degrees ",p.getx(),p.gety(),q.getx(),q.gety(),tr,td);
  p.put(-7*cos(PI/4), 7*sin(PI/4)); q.put(0,1);                    tr=pqtheta(p,q); td=tr*180/PI; dc= 45;         ASSERT_NEAR( td, dc, err * fabs(dc) );  //printf("p=(%6.3lf,%6.3lf) q=(%6.3lf,%6.3lf) theta=%6.2lf radians %6.2lf degrees ",p.getx(),p.gety(),q.getx(),q.gety(),tr,td);
  p.put(-7*cos(PI/4),-7*sin(PI/4)); q.put(0,1);                    tr=pqtheta(p,q); td=tr*180/PI; dc=135;         ASSERT_NEAR( td, dc, err * fabs(dc) );  //printf("p=(%6.3lf,%6.3lf) q=(%6.3lf,%6.3lf) theta=%6.2lf radians %6.2lf degrees ",p.getx(),p.gety(),q.getx(),q.gety(),tr,td);
}

//-----------------------------------------------------------------------------
//! \brief Half circle test
//-----------------------------------------------------------------------------
TEST( TestSuiteEllipse, circle )
{
  ellipsec circle(1,1);
  vectR2 q;
  vectR2 (ellipsec::*memfptr)(double t);
  
  memfptr = &ellipsec::xy;
  q=(circle.*memfptr)(1); printf("q=(%lf,%lf)\n",q.getx(),q.gety());

  const long N = 100000;
  const double length = circle.pathlength(M_PI/2,-M_PI/2, N);
  const double err = 1e-10;
  ASSERT_NEAR( length, M_PI, err * M_PI );
}

//-----------------------------------------------------------------------------
//! \brief Estimate perimeter length of unit circle test
//-----------------------------------------------------------------------------
TEST( TestSuiteCircle, circle )
{
  ellipsec circle(1,1);
  double length;
  const long N=1000;
  //-----------------------------------
  // Circle centered at (0,0) and with radius 1
  //-----------------------------------
  length  = circle.pathlength(0,2*M_PI,N); // estimate path length around circle
  ASSERT_NEAR( length, 2*M_PI, 1e-5 * (2*M_PI) ); // relative error
  //-----------------------------------
  // Circle centered at (0,0) and with radius 1
  //-----------------------------------
  circle.set(1,1,0,0,0);
  length  = circle.pathlength(0,2*M_PI,N);
  ASSERT_NEAR( length, 2*M_PI, 1e-5 * (2*M_PI) );
  //-----------------------------------
  // Circle centered at (1,-2) and with radius 1
  //-----------------------------------
  circle.set(1,1,0,1,-2);
  length  = circle.pathlength(0,2*M_PI,N);
  ASSERT_NEAR( length, 2*M_PI, 1e-5 * (2*M_PI) );
  //-----------------------------------
  // Circle centered at (-3,-5) and with radius 1
  // measuring half-way around
  //-----------------------------------
  circle.set(1,1,0,-3,-5);
  length  = circle.pathlength(M_PI,2*M_PI,N);
  ASSERT_NEAR( length, M_PI, 1e-5 * M_PI );
}

//-----------------------------------------------------------------------------
//! \brief Find points on unit circle that are distance 1 
//!        from the point at pi/4 radians on unit circle
//-----------------------------------------------------------------------------
TEST( TestSuiteEllipse, circled1 )
{
  ellipsec circle(1,1);
  double phi1, phi2;
  double err1, err2;
  const long N = 1000;
  const double err = 10.0 / N;

  circle.findt_dfroms( M_PI/4, 1, +1, N, &phi1, &err1);
  circle.findt_dfroms( M_PI/4, 1, -1, N, &phi2, &err2);
  ASSERT_NEAR( phi1, M_PI/4 + 1, err * (M_PI/4 + 1) );
  ASSERT_NEAR( phi2, M_PI/4 - 1, err * (1 - M_PI/4) );
}

//-----------------------------------------------------------------------------
//! \brief Test elliptic normalize
//-----------------------------------------------------------------------------
TEST( TestSuiteVectR2, rotate )
{
  vectR2 p,q;
  const double phi = M_PI / 2;
  const double err = 1e-15;
  p.put(0,2); q= p &   phi ;  ASSERT_NEAR( q.getx(), -2, err ); ASSERT_NEAR( q.gety(),  0, err * 2 );
  p.put(0,2); q= p & (-phi);  ASSERT_NEAR( q.getx(),  2, err ); ASSERT_NEAR( q.gety(),  0, err * 2 );
  p.put(2,0); q= p &   phi ;  ASSERT_NEAR( q.getx(),  0, err ); ASSERT_NEAR( q.gety(),  2, err * 2 );
  p.put(2,0); q= p & (-phi);  ASSERT_NEAR( q.getx(),  0, err ); ASSERT_NEAR( q.gety(), -2, err * 2 );
  p.put(2,0);    p &=  phi ;  ASSERT_NEAR( p.getx(),  0, err ); ASSERT_NEAR( p.gety(),  2, err * 2 );
}

//-----------------------------------------------------------------------------
//! \brief Find points on ellipse(0.8,1) that are distance 1 
//!        from the point at pi/2 radians on the ellipse
//-----------------------------------------------------------------------------
TEST( TestSuiteEllipse, distance1 )
{
  ellipsec ellipse(0.8,1);
  const double errMax = 0.0015;
  const long N = 1000;
  double phi1, phi2;
  double err1, err2;

  ellipse.findt_dfroms( M_PI/2, 1, +1, N, &phi1, &err1 );
  ellipse.findt_dfroms( M_PI/2, 1, -1, N, &phi2, &err2 );
//const double perror1 = fabs(100*err1);
//const double perror2 = fabs(100*err2);
//printf("(x,y) = (%9.6lf, %9.6lf) phi = %9.6lf error = %9.6lf (%lf%%)\n", ellipse.x(phi1), ellipse.y(phi1), phi1, err1, perror1 );
//printf("(x,y) = (%9.6lf, %9.6lf) phi = %9.6lf error = %9.6lf (%lf%%)\n", ellipse.x(phi2), ellipse.y(phi2), phi2, err2, perror2 );
  ASSERT_LE( err1, errMax );
  ASSERT_LE( err2, errMax );
}

//-----------------------------------------------------------------------------
//! \brief Find points on ellipse(0.8,1) that are distance 1 
//!        from the point at pi/2 radians on the ellipse
//-----------------------------------------------------------------------------
TEST( TestSuiteEllipse, distance1b )
{
  vectR2 p(cos(PI/4),sin(PI/4));
  ellipsec ellipse;
  double t;
  if( ellipse.setab_givenxyb(p,2) )
  {
    ASSERT_NE( ellipse.tgivenxy(p,&t), 0 );
    printf("t=%lf (x,y)=(%lf,%lf) ellipse(a,b)=(%lf,%lf)\n",t, p.getx(),p.gety(),ellipse.geta(),ellipse.getb());
  }
}

//-----------------------------------------------------------------------------
//! \brief Find the perimeter of an ellipse
//-----------------------------------------------------------------------------
TEST( TestSuiteEllipse, perimeter )
{
  ellipsec ellipse(0.8,2);
  const long   N     = 1000;
  const double err   = 1e-4;
  const double perim = ellipse.perimeter(N);
  const double est   = ellipse.estimate();
  ASSERT_NEAR( est, perim, err * perim );
//printf("ellipse(%lf,%lf) perimeter=%9.6lf estimate=%9.6lf\n",ellipse.geta(),ellipse.getb(),perim,est);
}

//-----------------------------------------------------------------------------
//! \brief Test DNA metric
//-----------------------------------------------------------------------------
TEST( TestSuiteDNA, metric ) 
{
  const char domain[]="ATCG";
  const int N=strlen(domain);
  int n, m;
  double d;

//for(m=0;m<N;m++)printf("     %c",domain[m]);//print column headers: 0 A ...
//putchar('\n');
  for( n=0; n<N; n++ )
  {
  //printf("%c  ",domain[n]);// print row headers
    for(m=0;m<N;m++)
    {
      d=dnan_metric(domain[n],domain[m]);
      if( m==n ) ASSERT_EQ( d, 0 );
      else       ASSERT_EQ( d, 1 );
    //if(d==(double)abs(d)) printf("%5.0lf ",d);//if d is integer
    //else                  printf("%5.3lf ",d);//if d is non-integer
    }
  //putchar('\n');
  }
}

/*-------------------------------------------------------------------------
 * test DNA-N metric
 *-------------------------------------------------------------------------*/
TEST( TestSuiteDNAN, metric ) 
{
  const char domain[]="ATCGN";
  const int N=strlen(domain);
  int n,m;
  double d;

//printf("dnan_metric(a,b) test over domain of N=%d symbols\n",N);
//printf("-------------------------------------------------\n");
//printf("  ");
//for( m=0; m<N; m++ ) printf("     %c",domain[m]);//print column headers: 0 A ...
//putchar('\n');
  for( n=0; n<N; n++ )
  {
  //printf("%c  ", domain[n] );// print row headers
    for( m=0; m<N; m++ )
    {
      d = dnan_metric( domain[n], domain[m] );
      if( m==n ) ASSERT_EQ( d, 0 );
      else       ASSERT_EQ( d, 1 );
    //if(d==(double)abs(d)) printf("%5.0lf ",d);//if d is integer
    //else                  printf("%5.3lf ",d);//if d is non-integer
    }
//putchar('\n');
  }
}

//-----------------------------------------------------------------------------
//! \brief Test balloon metric
//-----------------------------------------------------------------------------
TEST( TestSuiteMetric, balloon )
{
  vectR2 p,q;
  double d;
  const double err = 1e-6;

  printf("Balloon metric tests\n");
  printf("--------------------\n");
  p.put(0,2);                   q.put(2,0);                    d=metric_balloon(p,q);  ASSERT_NEAR( d, 2.000000, err * 2.000000 ); printf("d((%9.6lf,%9.6lf),(%9.6lf,%9.6lf))=%9.6lf\n",p.getx(),p.gety(),q.getx(),q.gety(),d);
  p.put(0,1);                   q.put(1,0);                    d=metric_balloon(p,q);  ASSERT_NEAR( d, 1.000000, err * 1.000000 ); printf("d((%9.6lf,%9.6lf),(%9.6lf,%9.6lf))=%9.6lf\n",p.getx(),p.gety(),q.getx(),q.gety(),d);
  p.put(0,1);                   q.put( cos(PI/4), sin(PI/4));  d=metric_balloon(p,q);  ASSERT_NEAR( d, 0.500000, err * 0.500000 ); printf("d((%9.6lf,%9.6lf),(%9.6lf,%9.6lf))=%9.6lf\n",p.getx(),p.gety(),q.getx(),q.gety(),d);
  p.put(0,1);                   q.put(-cos(PI/4), sin(PI/4));  d=metric_balloon(p,q);  ASSERT_NEAR( d, 0.500000, err * 0.500000 ); printf("d((%9.6lf,%9.6lf),(%9.6lf,%9.6lf))=%9.6lf\n",p.getx(),p.gety(),q.getx(),q.gety(),d);
  p.put(0,1);                   q.put(-cos(PI/4),-sin(PI/4));  d=metric_balloon(p,q);  ASSERT_NEAR( d, 1.500000, err * 1.500000 ); printf("d((%9.6lf,%9.6lf),(%9.6lf,%9.6lf))=%9.6lf\n",p.getx(),p.gety(),q.getx(),q.gety(),d);
  p.put(0,1);                   q.put( cos(PI/4),-sin(PI/4));  d=metric_balloon(p,q);  ASSERT_NEAR( d, 1.500000, err * 1.500000 ); printf("d((%9.6lf,%9.6lf),(%9.6lf,%9.6lf))=%9.6lf\n",p.getx(),p.gety(),q.getx(),q.gety(),d);
  p.put(1,0);                   q.put(-0.5,-0.5);              d=metric_balloon(p,q);  ASSERT_NEAR( d, 1.126357, err * 1.126357 ); printf("d((%9.6lf,%9.6lf),(%9.6lf,%9.6lf))=%9.6lf\n",p.getx(),p.gety(),q.getx(),q.gety(),d);
  p.put(1,0);                   q.put(-2,-2);                  d=metric_balloon(p,q);  ASSERT_NEAR( d, 2.388169, err * 2.388169 ); printf("d((%9.6lf,%9.6lf),(%9.6lf,%9.6lf))=%9.6lf\n",p.getx(),p.gety(),q.getx(),q.gety(),d);
  p.put(1,0);                   q.put(0,2);                    d=metric_balloon(p,q);  ASSERT_NEAR( d, 1.541964, err * 1.541964 ); printf("d((%9.6lf,%9.6lf),(%9.6lf,%9.6lf))=%9.6lf\n",p.getx(),p.gety(),q.getx(),q.gety(),d);
  p.put( cos(PI/4),-sin(PI/4)); q.put(-2,1);                   d=metric_balloon(p,q);  ASSERT_NEAR( d, 2.075940, err * 2.075940 ); printf("d((%9.6lf,%9.6lf),(%9.6lf,%9.6lf))=%9.6lf\n",p.getx(),p.gety(),q.getx(),q.gety(),d);
  p.put( cos(PI/4),-sin(PI/4)); q.put(-1.63,1.33);             d=metric_balloon(p,q);  ASSERT_NEAR( d, 1.980283, err * 1.980283 ); printf("d((%9.6lf,%9.6lf),(%9.6lf,%9.6lf))=%9.6lf\n",p.getx(),p.gety(),q.getx(),q.gety(),d);
}

//-----------------------------------------------------------------------------
//! \brief Test Mean Cicular Arc metric in R^2
//-----------------------------------------------------------------------------
TEST( TestSuiteMetric, mca )
{
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
}

//-----------------------------------------------------------------------------
//! \brief Test conjugate operation
//-----------------------------------------------------------------------------
TEST( TestSuiteGeneral, conj )
{
  complex z,zc;

  z.put(1.0,1.0);   
  zc=z.conj();  
  ASSERT_EQ(  z.getx(),  zc.getx() );
  ASSERT_EQ(  z.gety(), -zc.gety() );

  z.put(-3.2,4.7);  
  zc = z.conj();  
  ASSERT_EQ(  z.getx(),  zc.getx() );
  ASSERT_EQ(  z.gety(), -zc.gety() );
}

//-----------------------------------------------------------------------------
//! \brief Test spinner
//-----------------------------------------------------------------------------
TEST( TestSuiteMetric, spinner )
{
  const long N=160;
  spinseq x(N);
  seqR1 cseq(2*N+1);

  ASSERT_EQ( spin_metric( 'A'+0, 'A'+0 ), 0 ); ASSERT_EQ( spin_metric( 'A'+0, 'A'+1 ), 1 );  ASSERT_EQ( spin_metric( 'A'+0, 'A'+2 ), 2 ); ASSERT_EQ( spin_metric( 'A'+0, 'A'+3 ), 3 ); ASSERT_EQ( spin_metric( 'A'+0, 'A'+4 ), 2 ); ASSERT_EQ( spin_metric( 'A'+0, 'A'+5 ), 2 );
  ASSERT_EQ( spin_metric( 'A'+0, 'A'+0 ), 0 ); ASSERT_EQ( spin_metric( 'A'+0, 'A'+1 ), 1 );  ASSERT_EQ( spin_metric( 'A'+0, 'A'+2 ), 2 ); ASSERT_EQ( spin_metric( 'A'+0, 'A'+3 ), 3 ); ASSERT_EQ( spin_metric( 'A'+0, 'A'+4 ), 2 ); ASSERT_EQ( spin_metric( 'A'+0, 'A'+5 ), 2 );
  ASSERT_EQ( spin_metric( 'A'+1, 'A'+0 ), 1 ); ASSERT_EQ( spin_metric( 'A'+1, 'A'+1 ), 0 );  ASSERT_EQ( spin_metric( 'A'+1, 'A'+2 ), 1 ); ASSERT_EQ( spin_metric( 'A'+1, 'A'+3 ), 2 ); ASSERT_EQ( spin_metric( 'A'+1, 'A'+4 ), 3 ); ASSERT_EQ( spin_metric( 'A'+1, 'A'+5 ), 2 );
  ASSERT_EQ( spin_metric( 'A'+2, 'A'+0 ), 2 ); ASSERT_EQ( spin_metric( 'A'+2, 'A'+1 ), 1 );  ASSERT_EQ( spin_metric( 'A'+2, 'A'+2 ), 0 ); ASSERT_EQ( spin_metric( 'A'+2, 'A'+3 ), 1 ); ASSERT_EQ( spin_metric( 'A'+2, 'A'+4 ), 2 ); ASSERT_EQ( spin_metric( 'A'+2, 'A'+5 ), 3 );
  ASSERT_EQ( spin_metric( 'A'+3, 'A'+0 ), 3 ); ASSERT_EQ( spin_metric( 'A'+3, 'A'+1 ), 2 );  ASSERT_EQ( spin_metric( 'A'+3, 'A'+2 ), 1 ); ASSERT_EQ( spin_metric( 'A'+3, 'A'+3 ), 0 ); ASSERT_EQ( spin_metric( 'A'+3, 'A'+4 ), 1 ); ASSERT_EQ( spin_metric( 'A'+3, 'A'+5 ), 2 );
  ASSERT_EQ( spin_metric( 'A'+4, 'A'+0 ), 2 ); ASSERT_EQ( spin_metric( 'A'+4, 'A'+1 ), 3 );  ASSERT_EQ( spin_metric( 'A'+4, 'A'+2 ), 2 ); ASSERT_EQ( spin_metric( 'A'+4, 'A'+3 ), 1 ); ASSERT_EQ( spin_metric( 'A'+4, 'A'+4 ), 0 ); ASSERT_EQ( spin_metric( 'A'+4, 'A'+5 ), 1 );
  ASSERT_EQ( spin_metric( 'A'+5, 'A'+0 ), 2 ); ASSERT_EQ( spin_metric( 'A'+5, 'A'+1 ), 2 );  ASSERT_EQ( spin_metric( 'A'+5, 'A'+2 ), 3 ); ASSERT_EQ( spin_metric( 'A'+5, 'A'+3 ), 2 ); ASSERT_EQ( spin_metric( 'A'+5, 'A'+4 ), 1 ); ASSERT_EQ( spin_metric( 'A'+5, 'A'+5 ), 0 );

  x.randomize(0x5EED);
  x.list(          ); putchar('\n');
  cseq.list(0,  1*N); putchar('\n');
  cseq.list(N+1,2*N); putchar('\n');
  }

//-------------------------------------------------------------------------
//! \brief Test die seqR1 functions
//-------------------------------------------------------------------------*/
TEST( TestSuiteGeneral, die )
{
  const long N=100;
  dieseq x(N);

  printf("Test die seqR1 functions\n");
  printf("---------------------------\n");
  x.randomize(0x5EED);
  x.list();putchar('\n');
}

/*-------------------------------------------------------------------------
 * test real die seqR1 functions
 *-------------------------------------------------------------------------*/
TEST( TestSuiteDie, rdie )
{
  long m;
  const long N=100;
  rdieseq x(N);
  seqR1 Rxx(2*N+1);
  seqR1 histo(8);

  printf("Test real die seqR1 operations\n");
  x.randomize(0x5EED);
  x.list();putchar('\n');

  printf("Test real die histogram operation\n");
  histo=x.histogram();
  histo.list();putchar('\n');

  char a,b;
  for(a='A';a<='F';a++)
  {
    for(b='A';b<='F';b++)
    {
      if( a==b ) ASSERT_EQ( fdie_metric(a,b), 0.0 );
      else       ASSERT_EQ( fdie_metric(a,b), 1.0 );
    }
  }
  //x.metrictbl();

  //printf("Test auto-correlation Rxx\n");
  ASSERT_EQ( x.Rxx(-N), -2*N );
  ASSERT_EQ( x.Rxx( 0),    0 );
  ASSERT_EQ( x.Rxx( N), -2*N );
  for( m=-N; m< 0; m++ ) ASSERT_LT( x.Rxx(m), -(N-15) );
  for( m= 1; m<=N; m++ ) ASSERT_LT( x.Rxx(m), -(N-15) );

  if(x.Rxx(&Rxx,1))//auto-correlation seqR1 of truncated xR3hl
  // |      |   |____________switch to turn on counting display
  // |      |________________pointer to output correlation sequence
  // |_______________________input real die sequence
  {
    FAIL(); 
  } 

  //Rxx.list();

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

//-----------------------------------------------------------------------------
//! \brief Test expi function
//-----------------------------------------------------------------------------
TEST( TestSuite, expi )
{
  const long N=17;
  long n;
  complex y;
  double theta;
  for(n=0; n<=N; n++)
  {
    theta = 2.0 * M_PI*(double)n/(double)N;
    y = expi(theta);
  //printf("n=%2ld, theta=%lf (%3.0lf), expi(theta)=(%+9.6lf,%+9.6lf) mag=%lf\n",n,theta,theta/PI*180,y.real(),y.imag(), y.mag());
    ASSERT_DOUBLE_EQ( y.real(), cos(theta) );
    ASSERT_DOUBLE_EQ( y.imag(), sin(theta) );
    ASSERT_DOUBLE_EQ( y.mag(),  1.0        );
  }
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
  cos((2.0*M_PI/10.0),&y);
  dft(&y,&Dy);
  mag(&Dy,&mDy);
  printf("\n(y_n)=cos(2 pi n/10):\n");  y.list();
  printf("\nDFT(y_n):\n");              Dy.list();
  printf("\n|DFT(y_n)|:\n");            mDy.list();
  printf("\n");                         mDy.max('p');

  printf("\ntest DFT operations with sin(2pi n/10):\n");
  printf("---------------------------------------\n");
  sin((2.0*M_PI/10.0),&y);
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