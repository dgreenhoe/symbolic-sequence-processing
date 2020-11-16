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
#include "gtest/gtest.h"  // https://github.com/google/googletest/blob/master/googletest/docs/primer.md

//-----------------------------------------------------------------------------
// \brief Test opair Set operation
//-----------------------------------------------------------------------------
TEST( TestSuiteTuples, opairSet  )
{
  opair ab(3,-4);
  ASSERT_EQ( ab.getx(),  3 );
  ASSERT_EQ( ab.gety(), -4 );
}

//-----------------------------------------------------------------------------
// \brief Test opair Clear operation
//-----------------------------------------------------------------------------
TEST( TestSuiteTuples, opairClear  )
{
  opair ab(3,-4);
  ab.clear();
  ASSERT_EQ( ab.getx(),  0 );
  ASSERT_EQ( ab.gety(),  0 );
}

//-----------------------------------------------------------------------------
// \brief Test otriple clear operation
//-----------------------------------------------------------------------------
TEST( TestSuiteTuples, otripleSet  )
{
  otriple abc( 3, -4, 5 );
  ASSERT_EQ( abc.getx(),  3 );
  ASSERT_EQ( abc.gety(), -4 );
  ASSERT_EQ( abc.getz(),  5 );
}

//-----------------------------------------------------------------------------
// \brief Test otriple clear operation
//-----------------------------------------------------------------------------
TEST( TestSuiteTuples, otripleClear  )
{
  otriple abc(3,-4,5);
  abc.clear();
  ASSERT_EQ( abc.getx(), 0 );
  ASSERT_EQ( abc.gety(), 0 );
  ASSERT_EQ( abc.getz(), 0 );
}

//-----------------------------------------------------------------------------
// \brief Test osix Clear operation
//-----------------------------------------------------------------------------
TEST( TestSuiteTuples, osixSet )
{
  osix x(1,-2,3,-5,7,-11);
  osix y,z;
  y=x;
  z.put(y);
  ASSERT_EQ( x.get1(),  1 );
  ASSERT_EQ( x.get2(), -2 );
  ASSERT_EQ( x.get3(),  3 );
  ASSERT_EQ( x.get4(), -5 );
  ASSERT_EQ( x.get5(),  7 );
  ASSERT_EQ( x.get6(),-11 );
  ASSERT_EQ( x.get1(), y.get1() );
  ASSERT_EQ( x.get2(), y.get2() );
  ASSERT_EQ( x.get3(), y.get3() );
  ASSERT_EQ( x.get4(), y.get4() );
  ASSERT_EQ( x.get5(), y.get5() );
  ASSERT_EQ( x.get6(), y.get6() );
  ASSERT_EQ( x.get1(), z.get1() );
  ASSERT_EQ( x.get2(), z.get2() );
  ASSERT_EQ( x.get3(), z.get3() );
  ASSERT_EQ( x.get4(), z.get4() );
  ASSERT_EQ( x.get5(), z.get5() );
  ASSERT_EQ( x.get6(), z.get6() );
}

//-----------------------------------------------------------------------------
//! \brief Test osix Clear operation
//-----------------------------------------------------------------------------
TEST( TestSuiteTuples, osixClear  )
{
  osix x(1,-2,3,-5,7,-11);
  x.clear();
  ASSERT_EQ( x.get1(),  0 );
  ASSERT_EQ( x.get2(),  0 );
  ASSERT_EQ( x.get3(),  0 );
  ASSERT_EQ( x.get4(),  0 );
  ASSERT_EQ( x.get5(),  0 );
  ASSERT_EQ( x.get6(),  0 );
}

//-----------------------------------------------------------------------------
//! \brief Test otriple min operation
//-----------------------------------------------------------------------------
TEST( TestSuiteR3, min  )
{
  const otriple abc(1, -5, -3 );
  const double min = abc.min();
  ASSERT_EQ( min, -5 );
}

//-----------------------------------------------------------------------------
//! \brief Test otriple max operation
//-----------------------------------------------------------------------------
TEST( TestSuiteR3, max  )
{
  const otriple abc(1, -5, -3 );
  const double max = abc.max();
  ASSERT_EQ( max, 1 );
}

