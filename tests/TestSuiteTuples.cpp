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
//! \brief Test otriple angle operation
//-----------------------------------------------------------------------------
TEST( TestSuiteR3, angle )
{
  vectR3 p, q;
  p.put( 1, 0, 0 ); q.put( 0, 3, 0 ); ASSERT_DOUBLE_EQ( pqtheta(p,q)/M_PI*180., 90.0 );
  p.put( 1, 0, 0 ); q.put( 0, 0, 5 ); ASSERT_DOUBLE_EQ( pqtheta(p,q)/M_PI*180., 90.0 );
  p.put( 0, 3, 0 ); q.put( 1, 0, 0 ); ASSERT_DOUBLE_EQ( pqtheta(p,q)/M_PI*180., 90.0 );
  p.put( 0, 3, 0 ); q.put( 0, 0, 7 ); ASSERT_DOUBLE_EQ( pqtheta(p,q)/M_PI*180., 90.0 );
  p.put( 0, 0, 7 ); q.put( 2, 0, 0 ); ASSERT_DOUBLE_EQ( pqtheta(p,q)/M_PI*180., 90.0 );
  p.put( 0, 0, 7 ); q.put( 0, 5, 0 ); ASSERT_DOUBLE_EQ( pqtheta(p,q)/M_PI*180., 90.0 );
  p.put( 1, 2, 0 ); q.put(-2, 1, 0 ); ASSERT_DOUBLE_EQ( pqtheta(p,q)/M_PI*180., 90.0 );
  p.put( 0,-3,-5 ); q.put( 0, 5,-3 ); ASSERT_DOUBLE_EQ( pqtheta(p,q)/M_PI*180., 90.0 );
}

//-----------------------------------------------------------------------------
//! \brief Test otriple angle operation
//-----------------------------------------------------------------------------
TEST( TestSuiteR3, add )
{
  vectR3 p, q, z;
  p.put( 1, 0, 1  ); q.put( 0, 3, 2   ); z = p + q; ASSERT_EQ( z.getx(), 1+( 0  ) ); ASSERT_EQ( z.gety(),  0+(3  ) ); ASSERT_EQ( z.getz(),  1  +( 2  ) );
  p.put( 1, 0, 3  ); q.put( 0, 0, 5   ); z = p + q; ASSERT_EQ( z.getx(), 1+( 0  ) ); ASSERT_EQ( z.gety(),  0+(0  ) ); ASSERT_EQ( z.getz(),  3  +( 5  ) );
  p.put( 0, 3, 4  ); q.put( 1, 0, 9   ); z = p + q; ASSERT_EQ( z.getx(), 0+( 1  ) ); ASSERT_EQ( z.gety(),  3+(0  ) ); ASSERT_EQ( z.getz(),  4  +( 9  ) );
  p.put( 0, 3, 8  ); q.put( 0, 0, 7   ); z = p + q; ASSERT_EQ( z.getx(), 0+( 0  ) ); ASSERT_EQ( z.gety(),  3+(0  ) ); ASSERT_EQ( z.getz(),  8  +( 7  ) );
  p.put( 0, 0, 7  ); q.put( 2, 0, 0   ); z = p + q; ASSERT_EQ( z.getx(), 0+( 2  ) ); ASSERT_EQ( z.gety(),  0+(0  ) ); ASSERT_EQ( z.getz(),  7  +( 0  ) );
  p.put( 0, 0, 7  ); q.put( 0, 5, 1.1 ); z = p + q; ASSERT_EQ( z.getx(), 0+( 0  ) ); ASSERT_EQ( z.gety(),  0+(5  ) ); ASSERT_EQ( z.getz(),  7  +( 1.1) );
  p.put( 0, 0, 7  ); q.put( 0, 5.1, 1 ); z = p + q; ASSERT_EQ( z.getx(), 0+( 0  ) ); ASSERT_EQ( z.gety(),  0+(5.1) ); ASSERT_EQ( z.getz(),  7  +(1   ) );
  p.put( 0, 0, 7  ); q.put( 0.3, 5, 1 ); z = p + q; ASSERT_EQ( z.getx(), 0+( 0.3) ); ASSERT_EQ( z.gety(),  0+(5  ) ); ASSERT_EQ( z.getz(),  7  +(1   ) );
  p.put( 1, 2, 2.3); q.put(-2, 1, 5.6 ); z = p + q; ASSERT_EQ( z.getx(), 1+(-2  ) ); ASSERT_EQ( z.gety(),  2+(1  ) ); ASSERT_EQ( z.getz(),  2.3+( 5.6) );
  p.put( 0,-3,-5  ); q.put( 0, 5,-3   ); z = p + q; ASSERT_EQ( z.getx(), 0+( 0  ) ); ASSERT_EQ( z.gety(), -3+(5  ) ); ASSERT_EQ( z.getz(), -5  +(-3  ) );
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

