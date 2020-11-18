//=============================================================================
// Daniel J. Greenhoe
// Test routines
//=============================================================================
#include <stdio.h>
#include <math.h>
#include "r1.h"
#include "r2.h"
#include "r3.h"
#include "gtest/gtest.h"  // https://github.com/google/googletest/blob/master/googletest/docs/primer.md

//-----------------------------------------------------------------------------
// \brief Test opair Set operation
//-----------------------------------------------------------------------------
TEST( TestSuiteR2, set  )
{
  opair ab(3,-4);
  ASSERT_EQ( ab.getx(),  3 );
  ASSERT_EQ( ab.gety(), -4 );
}

//-----------------------------------------------------------------------------
// \brief Test opair Clear operation
//-----------------------------------------------------------------------------
TEST( TestSuiteR2, clear  )
{
  opair ab(3,-4);
  ab.clear();
  ASSERT_EQ( ab.getx(),  0 );
  ASSERT_EQ( ab.gety(),  0 );
}

//-----------------------------------------------------------------------------
//! \brief Test R2 operations
//-----------------------------------------------------------------------------
TEST( TestSuiteR2, misc )
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

