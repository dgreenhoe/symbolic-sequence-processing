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
// \brief Test otriple clear operation
//-----------------------------------------------------------------------------
TEST( TestSuiteR3, set  )
{
  otriple abc( 3, -4, 5 );
  ASSERT_EQ( abc.getx(),  3 );
  ASSERT_EQ( abc.gety(), -4 );
  ASSERT_EQ( abc.getz(),  5 );
}

//-----------------------------------------------------------------------------
// \brief Test otriple clear operation
//-----------------------------------------------------------------------------
TEST( TestSuiteR3, clear  )
{
  otriple abc(3,-4,5);
  abc.clear();
  ASSERT_EQ( abc.getx(), 0 );
  ASSERT_EQ( abc.gety(), 0 );
  ASSERT_EQ( abc.getz(), 0 );
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
//! \brief Test VectR3 add operation
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
//! \brief Test VectR3 subtract operation
//-----------------------------------------------------------------------------
TEST( TestSuiteR3, sub )
{
  vectR3 p, q, z;
  p.put( 1, 0, 1  ); q.put( 0, 3, 2   ); z = p - q; ASSERT_EQ( z.getx(), 1-( 0  ) ); ASSERT_EQ( z.gety(),  0-(3  ) ); ASSERT_EQ( z.getz(),  1  -( 2  ) );
  p.put( 1, 0, 3  ); q.put( 0, 0, 5   ); z = p - q; ASSERT_EQ( z.getx(), 1-( 0  ) ); ASSERT_EQ( z.gety(),  0-(0  ) ); ASSERT_EQ( z.getz(),  3  -( 5  ) );
  p.put( 0, 3, 4  ); q.put( 1, 0, 9   ); z = p - q; ASSERT_EQ( z.getx(), 0-( 1  ) ); ASSERT_EQ( z.gety(),  3-(0  ) ); ASSERT_EQ( z.getz(),  4  -( 9  ) );
  p.put( 0, 3, 8  ); q.put( 0, 0, 7   ); z = p - q; ASSERT_EQ( z.getx(), 0-( 0  ) ); ASSERT_EQ( z.gety(),  3-(0  ) ); ASSERT_EQ( z.getz(),  8  -( 7  ) );
  p.put( 0, 0, 7  ); q.put( 2, 0, 0   ); z = p - q; ASSERT_EQ( z.getx(), 0-( 2  ) ); ASSERT_EQ( z.gety(),  0-(0  ) ); ASSERT_EQ( z.getz(),  7  -( 0  ) );
  p.put( 0, 0, 7  ); q.put( 0, 5, 1.1 ); z = p - q; ASSERT_EQ( z.getx(), 0-( 0  ) ); ASSERT_EQ( z.gety(),  0-(5  ) ); ASSERT_EQ( z.getz(),  7  -( 1.1) );
  p.put( 0, 0, 7  ); q.put( 0, 5.1, 1 ); z = p - q; ASSERT_EQ( z.getx(), 0-( 0  ) ); ASSERT_EQ( z.gety(),  0-(5.1) ); ASSERT_EQ( z.getz(),  7  -(1   ) );
  p.put( 0, 0, 7  ); q.put( 0.3, 5, 1 ); z = p - q; ASSERT_EQ( z.getx(), 0-( 0.3) ); ASSERT_EQ( z.gety(),  0-(5  ) ); ASSERT_EQ( z.getz(),  7  -(1   ) );
  p.put( 1, 2, 2.3); q.put(-2, 1, 5.6 ); z = p - q; ASSERT_EQ( z.getx(), 1-(-2  ) ); ASSERT_EQ( z.gety(),  2-(1  ) ); ASSERT_EQ( z.getz(),  2.3-( 5.6) );
  p.put( 0,-3,-5  ); q.put( 0, 5,-3   ); z = p - q; ASSERT_EQ( z.getx(), 0-( 0  ) ); ASSERT_EQ( z.gety(), -3-(5  ) ); ASSERT_EQ( z.getz(), -5  -(-3  ) );
}

//-----------------------------------------------------------------------------
//! \brief Test VectR3 plus equals operation
//-----------------------------------------------------------------------------
TEST( TestSuiteR3, plusEquals )
{
  vectR3 p, q;
  p.put( 1, 0, 1  ); q.put( 0, 3, 2   ); p += q; ASSERT_EQ( p.getx(), 1+( 0  ) ); ASSERT_EQ( p.gety(),  0+(3  ) ); ASSERT_EQ( p.getz(),  1  +( 2  ) );
  p.put( 1, 0, 3  ); q.put( 0, 0, 5   ); p += q; ASSERT_EQ( p.getx(), 1+( 0  ) ); ASSERT_EQ( p.gety(),  0+(0  ) ); ASSERT_EQ( p.getz(),  3  +( 5  ) );
  p.put( 0, 3, 4  ); q.put( 1, 0, 9   ); p += q; ASSERT_EQ( p.getx(), 0+( 1  ) ); ASSERT_EQ( p.gety(),  3+(0  ) ); ASSERT_EQ( p.getz(),  4  +( 9  ) );
  p.put( 0, 3, 8  ); q.put( 0, 0, 7   ); p += q; ASSERT_EQ( p.getx(), 0+( 0  ) ); ASSERT_EQ( p.gety(),  3+(0  ) ); ASSERT_EQ( p.getz(),  8  +( 7  ) );
  p.put( 0, 0, 7  ); q.put( 2, 0, 0   ); p += q; ASSERT_EQ( p.getx(), 0+( 2  ) ); ASSERT_EQ( p.gety(),  0+(0  ) ); ASSERT_EQ( p.getz(),  7  +( 0  ) );
  p.put( 0, 0, 7  ); q.put( 0, 5, 1.1 ); p += q; ASSERT_EQ( p.getx(), 0+( 0  ) ); ASSERT_EQ( p.gety(),  0+(5  ) ); ASSERT_EQ( p.getz(),  7  +( 1.1) );
  p.put( 0, 0, 7  ); q.put( 0, 5.1, 1 ); p += q; ASSERT_EQ( p.getx(), 0+( 0  ) ); ASSERT_EQ( p.gety(),  0+(5.1) ); ASSERT_EQ( p.getz(),  7  +(1   ) );
  p.put( 0, 0, 7  ); q.put( 0.3, 5, 1 ); p += q; ASSERT_EQ( p.getx(), 0+( 0.3) ); ASSERT_EQ( p.gety(),  0+(5  ) ); ASSERT_EQ( p.getz(),  7  +(1   ) );
  p.put( 1, 2, 2.3); q.put(-2, 1, 5.6 ); p += q; ASSERT_EQ( p.getx(), 1+(-2  ) ); ASSERT_EQ( p.gety(),  2+(1  ) ); ASSERT_EQ( p.getz(),  2.3+( 5.6) );
  p.put( 0,-3,-5  ); q.put( 0, 5,-3   ); p += q; ASSERT_EQ( p.getx(), 0+( 0  ) ); ASSERT_EQ( p.gety(), -3+(5  ) ); ASSERT_EQ( p.getz(), -5  +(-3  ) );
}

//-----------------------------------------------------------------------------
//! \brief Test VectR3 plus equals operation
//-----------------------------------------------------------------------------
TEST( TestSuiteR3, subEquals )
{
  vectR3 p, q;
  p.put( 1, 0, 1  ); q.put( 0, 3, 2   ); p -= q; ASSERT_EQ( p.getx(), 1-( 0  ) ); ASSERT_EQ( p.gety(),  0-(3  ) ); ASSERT_EQ( p.getz(),  1  -( 2  ) );
  p.put( 1, 0, 3  ); q.put( 0, 0, 5   ); p -= q; ASSERT_EQ( p.getx(), 1-( 0  ) ); ASSERT_EQ( p.gety(),  0-(0  ) ); ASSERT_EQ( p.getz(),  3  -( 5  ) );
  p.put( 0, 3, 4  ); q.put( 1, 0, 9   ); p -= q; ASSERT_EQ( p.getx(), 0-( 1  ) ); ASSERT_EQ( p.gety(),  3-(0  ) ); ASSERT_EQ( p.getz(),  4  -( 9  ) );
  p.put( 0, 3, 8  ); q.put( 0, 0, 7   ); p -= q; ASSERT_EQ( p.getx(), 0-( 0  ) ); ASSERT_EQ( p.gety(),  3-(0  ) ); ASSERT_EQ( p.getz(),  8  -( 7  ) );
  p.put( 0, 0, 7  ); q.put( 2, 0, 0   ); p -= q; ASSERT_EQ( p.getx(), 0-( 2  ) ); ASSERT_EQ( p.gety(),  0-(0  ) ); ASSERT_EQ( p.getz(),  7  -( 0  ) );
  p.put( 0, 0, 7  ); q.put( 0, 5, 1.1 ); p -= q; ASSERT_EQ( p.getx(), 0-( 0  ) ); ASSERT_EQ( p.gety(),  0-(5  ) ); ASSERT_EQ( p.getz(),  7  -( 1.1) );
  p.put( 0, 0, 7  ); q.put( 0, 5.1, 1 ); p -= q; ASSERT_EQ( p.getx(), 0-( 0  ) ); ASSERT_EQ( p.gety(),  0-(5.1) ); ASSERT_EQ( p.getz(),  7  -(1   ) );
  p.put( 0, 0, 7  ); q.put( 0.3, 5, 1 ); p -= q; ASSERT_EQ( p.getx(), 0-( 0.3) ); ASSERT_EQ( p.gety(),  0-(5  ) ); ASSERT_EQ( p.getz(),  7  -(1   ) );
  p.put( 1, 2, 2.3); q.put(-2, 1, 5.6 ); p -= q; ASSERT_EQ( p.getx(), 1-(-2  ) ); ASSERT_EQ( p.gety(),  2-(1  ) ); ASSERT_EQ( p.getz(),  2.3-( 5.6) );
  p.put( 0,-3,-5  ); q.put( 0, 5,-3   ); p -= q; ASSERT_EQ( p.getx(), 0-( 0  ) ); ASSERT_EQ( p.gety(), -3-(5  ) ); ASSERT_EQ( p.getz(), -5  -(-3  ) );
}

//-----------------------------------------------------------------------------
//! \brief Test VectR3 plus equals operation
//-----------------------------------------------------------------------------
TEST( TestSuiteR3, negate )
{
  vectR3 p, q;
  p.put( 1, 0, 1  );  q = -p; ASSERT_EQ( q.getx(), -1 ); ASSERT_EQ( q.gety(), -0 ); ASSERT_EQ( q.getz(), -1   );
  p.put( 1, 0, 3  );  q = -p; ASSERT_EQ( q.getx(), -1 ); ASSERT_EQ( q.gety(), -0 ); ASSERT_EQ( q.getz(), -3   );
  p.put( 0, 3, 4  );  q = -p; ASSERT_EQ( q.getx(), -0 ); ASSERT_EQ( q.gety(), -3 ); ASSERT_EQ( q.getz(), -4   );
  p.put( 0, 3, 8  );  q = -p; ASSERT_EQ( q.getx(), -0 ); ASSERT_EQ( q.gety(), -3 ); ASSERT_EQ( q.getz(), -8   );
  p.put( 0, 0, 7  );  q = -p; ASSERT_EQ( q.getx(), -0 ); ASSERT_EQ( q.gety(), -0 ); ASSERT_EQ( q.getz(), -7   );
  p.put( 0, 0, 7  );  q = -p; ASSERT_EQ( q.getx(), -0 ); ASSERT_EQ( q.gety(), -0 ); ASSERT_EQ( q.getz(), -7   );
  p.put( 0, 0, 7  );  q = -p; ASSERT_EQ( q.getx(), -0 ); ASSERT_EQ( q.gety(), -0 ); ASSERT_EQ( q.getz(), -7   );
  p.put( 0, 0, 7  );  q = -p; ASSERT_EQ( q.getx(), -0 ); ASSERT_EQ( q.gety(), -0 ); ASSERT_EQ( q.getz(), -7   );
  p.put( 1, 2, 2.3);  q = -p; ASSERT_EQ( q.getx(), -1 ); ASSERT_EQ( q.gety(), -2 ); ASSERT_EQ( q.getz(), -2.3 );
  p.put( 0,-3,-5  );  q = -p; ASSERT_EQ( q.getx(), -0 ); ASSERT_EQ( q.gety(),  3 ); ASSERT_EQ( q.getz(),  5   );
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

