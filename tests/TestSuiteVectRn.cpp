//=============================================================================
// Daniel J. Greenhoe
// Test routines
//=============================================================================
#include <stdio.h>
#include <math.h>
#include "r1.h"
#include "r2.h"
#include "r3.h"
#include "r4.h"
#include "r6.h"
#include "gtest/gtest.h"

//-----------------------------------------------------------------------------
//! \brief Test R2 operations
//-----------------------------------------------------------------------------
TEST( TestSuiteVectRn, R2 )
{
  vectR2 p(3,-4);
  vectR2 q,s;
  ASSERT_EQ( p.getx(),  3 );
  ASSERT_EQ( p.gety(), -4 );
  ASSERT_EQ( p.norm(),  5 );
  p.put(sqrt(3),1);
  ASSERT_DOUBLE_EQ( p.theta()/M_PI*180, 30.0 );
  q=-p;
  ASSERT_EQ( q.getx(), -sqrt(3) );
  ASSERT_EQ( q.gety(), -1       );
  p.put(3,-5); q.put(-2,7); s=p+q;
  ASSERT_EQ( s.getx(), 1 );
  ASSERT_EQ( s.gety(), 2 );
  p.put(3,-5); q.put(-2,7); s=p-q;
  ASSERT_EQ( s.getx(),  5 );
  ASSERT_EQ( s.gety(),-12 );
  p.put(2,-3); q.put(-5,7);
  ASSERT_DOUBLE_EQ( p^q, -31 );
  q=p; q&=(M_PI/2);
  ASSERT_DOUBLE_EQ( q.getx(),  3 );
  ASSERT_DOUBLE_EQ( q.gety(),  2 );
  q=p; q&=(M_PI);
  ASSERT_DOUBLE_EQ( q.getx(), -2 );
  ASSERT_DOUBLE_EQ( q.gety(),  3 );
  p.clear();
  ASSERT_DOUBLE_EQ( p.getx(),  0 );
  ASSERT_DOUBLE_EQ( p.gety(),  0 );
}

//-----------------------------------------------------------------------------
//! \brief Test R6 operations
//-----------------------------------------------------------------------------
TEST( TestSuiteVectRn, R6 )
{
  vectR6 x(1,-2,3,-5,7,-11);
  vectR6 y,z,v;
  ASSERT_EQ( x.get1(),  1 );
  ASSERT_EQ( x.get2(), -2 );
  ASSERT_EQ( x.get3(),  3 );
  ASSERT_EQ( x.get4(), -5 );
  ASSERT_EQ( x.get5(),  7 );
  ASSERT_EQ( x.get6(),-11 );

  y=x;
  ASSERT_EQ( y.get1(), x.get1() );
  ASSERT_EQ( y.get2(), x.get2() );
  ASSERT_EQ( y.get3(), x.get3() );
  ASSERT_EQ( y.get4(), x.get4() );
  ASSERT_EQ( y.get5(), x.get5() );
  ASSERT_EQ( y.get6(), x.get6() );

  z.put(y);
  ASSERT_EQ( z.get1(), x.get1() );
  ASSERT_EQ( z.get2(), x.get2() );
  ASSERT_EQ( z.get3(), x.get3() );
  ASSERT_EQ( z.get4(), x.get4() );
  ASSERT_EQ( z.get5(), x.get5() );
  ASSERT_EQ( z.get6(), x.get6() );

  x.clear();
  ASSERT_EQ( x.get1(), 0 );
  ASSERT_EQ( x.get2(), 0 );
  ASSERT_EQ( x.get3(), 0 );
  ASSERT_EQ( x.get4(), 0 );
  ASSERT_EQ( x.get5(), 0 );
  ASSERT_EQ( x.get6(), 0 );

  x.put(1,2,3,4,5,6); 
  ASSERT_EQ( x.get1(), 1 );
  ASSERT_EQ( x.get2(), 2 );
  ASSERT_EQ( x.get3(), 3 );
  ASSERT_EQ( x.get4(), 4 );
  ASSERT_EQ( x.get5(), 5 );
  ASSERT_EQ( x.get6(), 6 );

  y.put(1,1./2.,1./3.,1./4.,1./5.,1./6.);
  ASSERT_EQ( y.get1(), 1./1. );
  ASSERT_EQ( y.get2(), 1./2. );
  ASSERT_EQ( y.get3(), 1./3. );
  ASSERT_EQ( y.get4(), 1./4. );
  ASSERT_EQ( y.get5(), 1./5. );
  ASSERT_EQ( y.get6(), 1./6. );

  z=x+y;   //z.list("z=x+y:       ");
  ASSERT_EQ( z.get1(), 1 + 1./1. );
  ASSERT_EQ( z.get2(), 2 + 1./2. );
  ASSERT_EQ( z.get3(), 3 + 1./3. );
  ASSERT_EQ( z.get4(), 4 + 1./4. );
  ASSERT_EQ( z.get5(), 5 + 1./5. );
  ASSERT_EQ( z.get6(), 6 + 1./6. );

  v=z-y;
  ASSERT_EQ( v.get1(), 1 + 1./1. - 1./1. );
  ASSERT_EQ( v.get2(), 2 + 1./2. - 1./2. );
  ASSERT_EQ( v.get3(), 3 + 1./3. - 1./3. );
  ASSERT_EQ( v.get4(), 4 + 1./4. - 1./4. );
  ASSERT_EQ( v.get5(), 5 + 1./5. - 1./5. );
  ASSERT_EQ( v.get6(), 6 + 1./6. - 1./6. );

  v+=y;
  ASSERT_EQ( v.get1(), 1 + 1./1. - 1./1. + 1./1. );
  ASSERT_EQ( v.get2(), 2 + 1./2. - 1./2. + 1./2. );
  ASSERT_EQ( v.get3(), 3 + 1./3. - 1./3. + 1./3. );
  ASSERT_EQ( v.get4(), 4 + 1./4. - 1./4. + 1./4. );
  ASSERT_EQ( v.get5(), 5 + 1./5. - 1./5. + 1./5. );
  ASSERT_EQ( v.get6(), 6 + 1./6. - 1./6. + 1./6. );

  v-=y;
  ASSERT_EQ( v.get1(), 1 + 1./1. - 1./1. + 1./1.  - 1./1. );
  ASSERT_EQ( v.get2(), 2 + 1./2. - 1./2. + 1./2.  - 1./2. );
  ASSERT_EQ( v.get3(), 3 + 1./3. - 1./3. + 1./3.  - 1./3. );
  ASSERT_EQ( v.get4(), 4 + 1./4. - 1./4. + 1./4.  - 1./4. );
  ASSERT_EQ( v.get5(), 5 + 1./5. - 1./5. + 1./5.  - 1./5. );
  ASSERT_EQ( v.get6(), 6 + 1./6. - 1./6. + 1./6.  - 1./6. );

  y=2*v;   //y.list("y=2*v:       ");
  ASSERT_EQ( y.get1(), 2 * (1 + 1./1. - 1./1. + 1./1.  - 1./1.) );
  ASSERT_EQ( y.get2(), 2 * (2 + 1./2. - 1./2. + 1./2.  - 1./2.) );
  ASSERT_EQ( y.get3(), 2 * (3 + 1./3. - 1./3. + 1./3.  - 1./3.) );
  ASSERT_EQ( y.get4(), 2 * (4 + 1./4. - 1./4. + 1./4.  - 1./4.) );
  ASSERT_EQ( y.get5(), 2 * (5 + 1./5. - 1./5. + 1./5.  - 1./5.) );
  ASSERT_EQ( y.get6(), 2 * (6 + 1./6. - 1./6. + 1./6.  - 1./6.) );

  v*=-3;   //v.list("v*=-3:       ");
  ASSERT_EQ( v.get1(), -3 * (1 + 1./1. - 1./1. + 1./1.  - 1./1.) );
  ASSERT_EQ( v.get2(), -3 * (2 + 1./2. - 1./2. + 1./2.  - 1./2.) );
  ASSERT_EQ( v.get3(), -3 * (3 + 1./3. - 1./3. + 1./3.  - 1./3.) );
  ASSERT_EQ( v.get4(), -3 * (4 + 1./4. - 1./4. + 1./4.  - 1./4.) );
  ASSERT_EQ( v.get5(), -3 * (5 + 1./5. - 1./5. + 1./5.  - 1./5.) );
  ASSERT_EQ( v.get6(), -3 * (6 + 1./6. - 1./6. + 1./6.  - 1./6.) );
}

