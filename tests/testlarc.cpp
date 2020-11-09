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
  printf("number of fails = %d\n",fails);
  return fails;
  }

//-----------------------------------------------------------------------------
//! \brief Test Lagrange arc metric in R^2
//-----------------------------------------------------------------------------
TEST( TestSuiteLarc, R2 )
{
  vectR2 p,q;
  double d,dN;

  p.put( 0,1);        q.put( 1, 0); d=larc_metric(p,q); dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, 1e-6 * d );
  p.put( 0,1);        q.put(-1, 0); d=larc_metric(p,q); dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, 1e-6 * d );
  p.put( 0,1);        q.put( 0,-1); d=larc_metric(p,q); dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, 1e-6 * d );
  p.put( 1,0);        q.put( 0,-1); d=larc_metric(p,q); dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, 1e-6 * d );
  p.put( 1,0);        q.put(-1, 0); d=larc_metric(p,q); dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, 1e-6 * d );
  p.put(-1,0);        q.put( 0,-1); d=larc_metric(p,q); dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, 1e-6 * d );
  p.put( 0,1);        q.put( 2, 0); d=larc_metric(p,q); dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, 1e-6 * d );
  p.put( 0,1);        q.put( 0,-2); d=larc_metric(p,q); dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, 1e-6 * d );
  p.put( 0,1);        q.put(-2, 1); d=larc_metric(p,q); dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, 1e-6 * d );
  p.put( 0.000001,0); q.put(PI, 0); d=larc_metric(p,q); dN=larc_metric(p,q,1000);  ASSERT_NEAR( dN, d, 1e-6 * d );
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
