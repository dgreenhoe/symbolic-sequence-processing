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
#include "larc.h"
#include "euclid.h"
#include "testlarc.h"
#include "gtest/gtest.h"  // https://github.com/google/googletest/blob/master/googletest/docs/primer.md

//-----------------------------------------------------------------------------
//! \brief Test Lagrange arc metric in R^2
//! \details See Example 3.12: Lagrange arc distance in R^2
//-----------------------------------------------------------------------------
TEST( TestSuiteLarc, R2 )
{
  vectR2 p,q;
  double d,dN;
  const double err = 1e-6;

  p.put(1,0);                   q.put(-1,0);                   d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, err * d );
  q.put(1,0);                   p.put(-1,0);                   d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, err * d );
  p.put(0,2);                   q.put(2,0);                    d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, err * d );
  p.put(0,1);                   q.put(1,0);                    d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, err * d );
  p.put(0,1);                   q.put( cos(PI/4), sin(PI/4));  d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, err * d );
  p.put(0,1);                   q.put(-cos(PI/4), sin(PI/4));  d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, err * d );
  p.put(0,1);                   q.put(-cos(PI/4),-sin(PI/4));  d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, err * d );
  p.put(0,1);                   q.put( cos(PI/4),-sin(PI/4));  d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, err * d );
  p.put(1,0);                   q.put(-0.5,-0.5);              d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, err * d );
  p.put(1,0);                   q.put(-2,-2);                  d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, err * d );
  p.put(1,0);                   q.put(0,2);                    d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, err * d );
  p.put( cos(PI/4),-sin(PI/4)); q.put(-2,1);                   d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, err * d );
  p.put( cos(PI/4),-sin(PI/4)); q.put(-1.63,1.33);             d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, err * d );
  p.put( cos(PI/4), sin(PI/4)); q.put( 1, 1);                  d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, err * d );
  p.put(1,0);                   q.put(2,0);                    d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, err * d );
  p.put(1,0);                   q.put(-1,0);                   d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, err * d );
  p.put(1,1);                   q.put(0,0);                    d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, err * d );
  p.put(1,0);                   q.put(0,-1);                   d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, err * d );
  p.put(0,1);                   q.put(-cos(PI/4), sin(PI/4));  d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, err * d );
  p.put(0,1);                   q.put( cos(PI/4), sin(PI/4));  d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, err * d );
  p.put(1,0);                   q.put(-2,0);                   d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, err * d );
  q.put(1,0);                   p.put(-2,0);                   d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, err * d );
  p.put(1,0);                   q.put(-0.5,0);                 d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, err * d );
//printf("Theorem 3.9 (1): triangle inequality fails\n");
  p.put(1,0);                   q.put(-0.5,0);                 d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, err * d );
  p.put(1,0);                   q.put(-0.5,0.2);               d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, err * d );
  p.put(-0.5,0.2);              q.put(-0.5,0);                 d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, err * d );
//printf("Theorem 3.9 (2): translation invariance fails\n");
  p.put(1,0.5);                 q.put(0.5,1);                  d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, err * d );
  p.put(0.5,0);                 q.put(0,0.5);                  d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, err * d );
//printf("Theorem 3.9 (5): balls are not convex\n");
  p.put(0,1);                   q.put(-0.70,-1.12);            d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, err * d );
  q.put(0,1);                   p.put(-0.70,-1.12);            d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, err * d );
  p.put(0,1);                   q.put( 0.70,-1.12);            d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, err * d );
  q.put(0,1);                   p.put( 0.70,-1.12);            d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, err * d );
  p.put(0,1);                   q.put( 0,   -1.12);            d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, err * d );
  q.put(0,1);                   p.put( 0,   -1.12);            d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, err * d );
//printf("Theorem 3.10: Lagrange arc distance versus Euclidean metric\n");
  p.put(1,0);                   q.put(-0.50,0);                d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, err * d );
  q.put(1,0);                   p.put(-0.50,0);                d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, err * d );
  p.put(1,0);                   q.put(-0.50,0.75);             d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, err * d );
  q.put(1,0);                   p.put(-0.50,0.75);             d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, err * d );
//printof(lptr,"Example 3.12: Lagrange arc distance in R^2 (?)\n");
  p.put( 0,1);                  q.put( 1, 0);                  d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, err * d );
  p.put( 0,1);                  q.put(-1, 0);                  d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, err * d );
  p.put( 0,1);                  q.put( 0,-1);                  d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, err * d );
  p.put( 1,0);                  q.put( 0,-1);                  d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, err * d );
  p.put( 1,0);                  q.put(-1, 0);                  d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, err * d );
  p.put(-1,0);                  q.put( 0,-1);                  d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, err * d );
  p.put( 0,1);                  q.put( 2, 0);                  d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, err * d );
  p.put( 0,1);                  q.put( 0,-2);                  d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, err * d );
  p.put( 0,1);                  q.put(-2, 1);                  d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, err * d );
  p.put( 0.000001,0);           q.put(PI, 0);                  d=larc_metric(p,q);  dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, err * d );
}

//-----------------------------------------------------------------------------
//! \brief Test Lagrange arc metric in R^3
//! \todo Try to tighten error bounds
//-----------------------------------------------------------------------------
TEST( TestSuiteLarc, R3 )
{
  vectR3 p,q;
  double d,de,s=1.0/PI;

  p.put(0,2,0);                   q.put(2,0,0);                    d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(0,1,0);                   q.put(1,0,0);                    d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(0,1,0);                   q.put( cos(PI/4), sin(PI/4),0);  d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(0,1,0);                   q.put(-cos(PI/4), sin(PI/4),0);  d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(0,1,0);                   q.put(-cos(PI/4),-sin(PI/4),0);  d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(0,1,0);                   q.put( cos(PI/4),-sin(PI/4),0);  d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(1,0,0);                   q.put(-0.5,-0.5,0);              d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(1,0,0);                   q.put(-2,-2,0);                  d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(1,0,0);                   q.put(0,2,0);                    d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( cos(PI/4),-sin(PI/4),0); q.put(-2,1,0);                   d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( cos(PI/4),-sin(PI/4),0); q.put(-1.63,1.33,0);             d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 1.0 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( cos(PI/4), sin(PI/4),0); q.put( 1, 1,0);                  d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 1.0 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(1,0,0);                   q.put(2,0,0);                    d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 1.0 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(1,0,0);                   q.put(-1,0,0);                   d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 1.0 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(1,1,0);                   q.put(0,0,0);                    d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 1.0 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(1,0,0);                   q.put(0,-1,0);                   d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 1.0 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);

  p.put( 0, 1, 0);                q.put( 0, 1, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 1, 0);                q.put( 0, 0, 1);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 1, 0);                q.put(-1, 0, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 1, 0);                q.put( 0,-1, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 1, 0);                q.put( 0, 0,-1);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 2, 0);                q.put( 1, 0, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 2, 0);                q.put( 0, 1, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 2, 0);                q.put( 0, 0, 1);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 2, 0);                q.put(-1, 0, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 2, 0);                q.put( 0,-1, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 1.0 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 2, 0);                q.put( 0, 0,-1);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 1.0 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 1, 0);                q.put( 1, 1, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 1.0 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 1, 0);                q.put( 0, 1, 1);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 1.0 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 1, 0);                q.put( 1, 0, 1);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 1.0 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 1, 0);                q.put(-1,-1, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 1.0 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 1, 0);                q.put( 0,-1,-1);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 1.0 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 1, 0);                q.put(-1, 0,-1);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 1.0 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1, 1, 1);                q.put(-0.5, 0.25,-2); d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 1.0 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1,-1,-2);                q.put( 1, 0, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 1.0 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1,-1,-2);                q.put( 0, 1, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 1.0 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1,-1,-2);                q.put( 0, 0, 1);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 1.0 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1,-1,-2);                q.put( 0, 0,-1);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 1.0 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1,-1,-2);                q.put( 0,-1, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 1.0 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1,-1,-2);                q.put(-1, 0, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 1.0 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 0, 0);                q.put( 1, 0, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 1.0 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 0, 0);                q.put( 0, 1, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 1.0 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 0, 0);                q.put( 0, 0, 1);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 1.0 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 0, 0);                q.put( 0, 0,-1);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 1.0 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 0, 0);                q.put( 0,-1, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 1.0 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 0, 0);                q.put(-1, 0, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 1.0 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1, 1,-2);                q.put( 1, 0, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 1.0 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1, 1,-2);                q.put( 0, 1, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 1.0 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1, 1,-2);                q.put( 0, 0, 1);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1, 1,-2);                q.put( 0, 0,-1);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1, 1,-2);                q.put( 0,-1, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1, 1,-2);                q.put(-1, 0, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1, 1, 0);                q.put( 1, 0, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1, 1, 0);                q.put( 0, 1, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1, 1, 0);                q.put( 0, 0, 1);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1, 1, 0);                q.put( 0, 0,-1);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1, 1, 0);                q.put( 0,-1, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1, 1, 0);                q.put(-1, 0, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 2, 2, 0);                q.put( 1, 0, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 2, 2, 0);                q.put( 0, 1, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 2, 2, 0);                q.put( 0, 0, 1);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 2, 2, 0);                q.put( 0, 0,-1);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 2, 2, 0);                q.put( 0,-1, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 2, 2, 0);                q.put(-1, 0, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1, 2, 1);                q.put( 1, 0, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1, 2, 1);                q.put( 0, 1, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1, 2, 1);                q.put( 0, 0, 1);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1, 2, 1);                q.put( 0, 0,-1);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1, 2, 1);                q.put( 0,-1, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1, 2, 1);                q.put(-1, 0, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 3, 1);                q.put( 1, 0, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 3, 1);                q.put( 0, 1, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 3, 1);                q.put( 0, 0, 1);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 3, 1);                q.put( 0, 0,-1);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 3, 1);                q.put( 0,-1, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 1.0 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 3, 1);                q.put(-1, 0, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 1.0 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 2, 0);                q.put( 1, 0, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 1.0 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 2, 0);                q.put( 0, 1, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 1.0 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 2, 0);                q.put( 0, 0, 1);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 1.0 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 2, 0);                q.put( 0, 0,-1);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 1.0 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 2, 0);                q.put( 0,-1, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 1.0 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 2, 0);                q.put(-1, 0, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 1.0 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(-1, 2, 1);                q.put( 1, 0, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 1.0 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(-1, 2, 1);                q.put( 0, 1, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 1.0 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(-1, 2, 1);                q.put( 0, 0, 1);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 1.0 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(-1, 2, 1);                q.put( 0, 0,-1);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 1.0 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(-1, 2, 1);                q.put( 0,-1, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 1.0 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(-1, 2, 1);                q.put(-1, 0, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 1.0 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);

  p.put(-1, 0, 0);                q.put( 1, 0, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(-1.01, 0, 0);             q.put(0.99, 0, 0);    d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(-1.001, 0, 0);            q.put(0.999, 0, 0);   d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(-1.0001, 0, 0);           q.put(0.9999, 0, 0);  d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);

  q.put(-1, 0, 0);                p.put( 1, 0, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  q.put(-1.01, 0, 0);             p.put(0.99, 0, 0);    d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  q.put(-1.001, 0, 0);            p.put(0.999, 0, 0);   d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  q.put(-1.0001, 0, 0);           p.put(0.9999, 0, 0);  d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);

  p.put(0,1,0);                             q.put( cos(PI/4),           sin(PI/4),0);    d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(0,1,0);                             q.put(-cos(PI/4),           sin(PI/4),0);    d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(0,1,0);                             q.put( cos(PI/6),           sin(PI/6),0);    d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(0,1,0);                             q.put(-cos(PI/6),           sin(PI/6),0);    d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( cos(PI/4),         sin(PI/4),0);   q.put( 0.5*cos(PI/8),   0.5*sin(PI/8),0);    d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(-cos(PI/4),         sin(PI/4),0);   q.put(-0.5*cos(3*PI/8), 0.5*sin(3*PI/8),0);  d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0.5*cos(PI/4), 0.5*sin(PI/4),0);   q.put( cos(PI/8),           sin(PI/8),0);    d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put(-0.5*cos(PI/4), 0.5*sin(PI/4),0);   q.put(-cos(3*PI/8),         sin(3*PI/8),0);  d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);

  p.put( 0, 1, 0);  q.put( 2, 0, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 0.5 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 1, 0);  q.put( 0,-2, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 1.0 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 1, 0);  q.put(-2, 1, 0);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 1.0 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 0, 1, 0);  q.put(-1, 0,-1);      d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 1.0 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
  p.put( 1, 1, 1);  q.put(-0.5,0.25,-2);  d=larc_metric(p,q); de=ae_metric(s,p,q); ASSERT_NEAR( d, de, 1.0 ); //printf("d((%6.3lf,%6.3lf,%6.3lf),(%6.3lf,%6.3lf,%6.3lf))=%12.9lf ae=%12.9lf\n",p.getx(),p.gety(),p.getz(),q.getx(),q.gety(),q.getz(),d,de);
}

//-----------------------------------------------------------------------------
//! \brief Test Lagrange arc metric in R^6
//-----------------------------------------------------------------------------
TEST( TestSuiteLarc, R6 )
{
  vectR6 p,q;
  double d;
  const double err = 1e-3;

  p.put( 1, 0, 0, 0, 0, 0);  q.put( 0, 0, 0, 0, 0, 1); d=larc_metric(p,q);  ASSERT_EQ( d, 0.5 ); ASSERT_EQ( pqtheta(p,q)*180/PI,  90 ); //printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put( 0, 1, 0, 0, 0, 0);  q.put( 0, 0, 0, 0, 1, 0); d=larc_metric(p,q);  ASSERT_EQ( d, 0.5 ); ASSERT_EQ( pqtheta(p,q)*180/PI,  90 ); //printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put( 0, 0, 1, 0, 0, 0);  q.put( 0, 0, 0, 1, 0, 0); d=larc_metric(p,q);  ASSERT_EQ( d, 0.5 ); ASSERT_EQ( pqtheta(p,q)*180/PI,  90 ); //printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put( 0, 0, 0, 1, 0, 0);  q.put( 0, 0, 1, 0, 0, 0); d=larc_metric(p,q);  ASSERT_EQ( d, 0.5 ); ASSERT_EQ( pqtheta(p,q)*180/PI,  90 ); //printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put( 0, 0, 0, 0, 1, 0);  q.put( 0, 1, 0, 0, 0, 0); d=larc_metric(p,q);  ASSERT_EQ( d, 0.5 ); ASSERT_EQ( pqtheta(p,q)*180/PI,  90 ); //printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put( 0, 0, 0, 0, 0, 1);  q.put( 1, 0, 0, 0, 0, 0); d=larc_metric(p,q);  ASSERT_EQ( d, 0.5 ); ASSERT_EQ( pqtheta(p,q)*180/PI,  90 ); //printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  
  p.put( 1, 0, 0, 0, 0, 0);  q.put(-1, 0, 0, 0, 0, 0); d=larc_metric(p,q);  ASSERT_EQ( d, 1.0 ); ASSERT_EQ( pqtheta(p,q)*180/PI, 180 ); //printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put( 0, 1, 0, 0, 0, 0);  q.put( 0,-1, 0, 0, 0, 0); d=larc_metric(p,q);  ASSERT_EQ( d, 1.0 ); ASSERT_EQ( pqtheta(p,q)*180/PI, 180 ); //printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put( 0, 0, 1, 0, 0, 0);  q.put( 0, 0,-1, 0, 0, 0); d=larc_metric(p,q);  ASSERT_EQ( d, 1.0 ); ASSERT_EQ( pqtheta(p,q)*180/PI, 180 ); //printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put( 0, 0, 0, 1, 0, 0);  q.put( 0, 0, 0,-1, 0, 0); d=larc_metric(p,q);  ASSERT_EQ( d, 1.0 ); ASSERT_EQ( pqtheta(p,q)*180/PI, 180 ); //printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put( 0, 0, 0, 0, 1, 0);  q.put( 0, 0, 0, 0,-1, 0); d=larc_metric(p,q);  ASSERT_EQ( d, 1.0 ); ASSERT_EQ( pqtheta(p,q)*180/PI, 180 ); //printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put( 0, 0, 0, 0, 0, 1);  q.put( 0, 0, 0, 0, 0,-1); d=larc_metric(p,q);  ASSERT_EQ( d, 1.0 ); ASSERT_EQ( pqtheta(p,q)*180/PI, 180 ); //printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  
  p.put( 1, 0, 0, 0, 0, 0);  q.put( 0, 0, 0, 0, 0, 1); d=larc_metric(p,q);  ASSERT_EQ( d, 0.5 ); ASSERT_EQ( pqtheta(p,q)*180/PI,  90 ); //printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put( 1, 0, 0, 0, 0, 0);  q.put( 0, 0, 0, 0, 1, 0); d=larc_metric(p,q);  ASSERT_EQ( d, 0.5 ); ASSERT_EQ( pqtheta(p,q)*180/PI,  90 ); //printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put( 1, 0, 0, 0, 0, 0);  q.put( 0, 0, 0, 1, 0, 0); d=larc_metric(p,q);  ASSERT_EQ( d, 0.5 ); ASSERT_EQ( pqtheta(p,q)*180/PI,  90 ); //printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put( 1, 0, 0, 0, 0, 0);  q.put( 0, 0, 1, 0, 0, 0); d=larc_metric(p,q);  ASSERT_EQ( d, 0.5 ); ASSERT_EQ( pqtheta(p,q)*180/PI,  90 ); //printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put( 1, 0, 0, 0, 0, 0);  q.put( 0, 1, 0, 0, 0, 0); d=larc_metric(p,q);  ASSERT_EQ( d, 0.5 ); ASSERT_EQ( pqtheta(p,q)*180/PI,  90 ); //printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put( 1, 0, 0, 0, 0, 0);  q.put( 1, 0, 0, 0, 0, 0); d=larc_metric(p,q);  ASSERT_EQ( d, 0.0 ); ASSERT_EQ( pqtheta(p,q)*180/PI,   0 ); //printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  //! \todo Assertion targets correct?
  p.put(-1, 2, 1,-1, 2, 1);  q.put( 1, 0, 0, 1, 0, 0); d=larc_metric(p,q);  ASSERT_NEAR( d, 1.684968, err*d ); ASSERT_NEAR( pqtheta(p,q)*180/PI, 114.09, err*114.09 ); //printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put(-1, 2, 1,-1, 2, 1);  q.put( 0, 1, 0, 0, 1, 0); d=larc_metric(p,q);  ASSERT_NEAR( d, 0.814214, err*d ); ASSERT_NEAR( pqtheta(p,q)*180/PI, 35.26 , err*35.26  ); //printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put(-1, 2, 1,-1, 2, 1);  q.put( 0, 0, 1, 0, 0, 1); d=larc_metric(p,q);  ASSERT_NEAR( d, 1.113722, err*d ); ASSERT_NEAR( pqtheta(p,q)*180/PI, 65.91 , err*65.91  ); //printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put(-1, 2, 1,-1, 2, 1);  q.put( 0, 0,-1, 0, 0,-1); d=larc_metric(p,q);  ASSERT_NEAR( d, 1.684968, err*d ); ASSERT_NEAR( pqtheta(p,q)*180/PI, 114.09, err*114.09 ); //printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put(-1, 2, 1,-1, 2, 1);  q.put( 0,-1, 0, 0,-1, 0); d=larc_metric(p,q);  ASSERT_NEAR( d, 2.072943, err*d ); ASSERT_NEAR( pqtheta(p,q)*180/PI, 144.74, err*144.74 ); //printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put(-1, 2, 1,-1, 2, 1);  q.put(-1, 0, 0,-1, 0, 0); d=larc_metric(p,q);  ASSERT_NEAR( d, 1.113722, err*d ); ASSERT_NEAR( pqtheta(p,q)*180/PI, 65.91 , err*65.91  ); //printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  //! \todo Assertion targets correct?
  p.put(0.458,0.022,0.120,0.234,0.167,0.000);  q.put(1, 0, 0, 0, 0, 0); d=larc_metric(p,q); ASSERT_NEAR( d, 0.205748, err*d ); ASSERT_NEAR( pqtheta(p,q)*180/PI, 34.29, err*34.29 ); //printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put(0.458,0.022,0.120,0.234,0.167,0.000);  q.put(0, 1, 0, 0, 0, 0); d=larc_metric(p,q); ASSERT_NEAR( d, 0.405079, err*d ); ASSERT_NEAR( pqtheta(p,q)*180/PI, 87.73, err*87.73 ); //printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put(0.458,0.022,0.120,0.234,0.167,0.000);  q.put(0, 0, 1, 0, 0, 0); d=larc_metric(p,q); ASSERT_NEAR( d, 0.364101, err*d ); ASSERT_NEAR( pqtheta(p,q)*180/PI, 77.50, err*77.50 ); //printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put(0.458,0.022,0.120,0.234,0.167,0.000);  q.put(0, 0, 0, 1, 0, 0); d=larc_metric(p,q); ASSERT_NEAR( d, 0.315300, err*d ); ASSERT_NEAR( pqtheta(p,q)*180/PI, 65.03, err*65.03 ); //printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put(0.458,0.022,0.120,0.234,0.167,0.000);  q.put(0, 0, 0, 0, 1, 0); d=larc_metric(p,q); ASSERT_NEAR( d, 0.344228, err*d ); ASSERT_NEAR( pqtheta(p,q)*180/PI, 72.47, err*72.47 ); //printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
  p.put(0.458,0.022,0.120,0.234,0.167,0.000);  q.put(0, 0, 0, 0, 0, 1); d=larc_metric(p,q); ASSERT_NEAR( d, 0.414280, err*d ); ASSERT_NEAR( pqtheta(p,q)*180/PI, 90.00, err*90.00 ); //printf("d((%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf),(%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf,%3.1lf))=%9.6lf t=%6.2lf\n",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6(),q.get1(),q.get2(),q.get3(),q.get4(),q.get5(),q.get6(),d,pqtheta(p,q)*180/PI);
}
