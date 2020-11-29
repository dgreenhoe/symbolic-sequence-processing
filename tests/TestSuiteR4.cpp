//=============================================================================
// Daniel J. Greenhoe
// Test routines
//=============================================================================
#include <stdio.h>
#include <math.h>
#include "r1.h"
#include "r4.h"
#include "gtest/gtest.h"

//-----------------------------------------------------------------------------
//! \brief Test R4 set
//-----------------------------------------------------------------------------
TEST( TestSuiteR4, set )
{
  const vectR4 x( 1, -2, 3, -5 );
  ASSERT_EQ( x.get1(),  1 );
  ASSERT_EQ( x.get2(), -2 );
  ASSERT_EQ( x.get3(),  3 );
  ASSERT_EQ( x.get4(), -5 );
}

//-----------------------------------------------------------------------------
//! \brief Test set-equal operation
//-----------------------------------------------------------------------------
TEST( TestSuiteR4, seto )
{
  oquad x( 1, -2, 3, -5 );
  oquad y,z;
  y=x;
  z.put(y);
  ASSERT_EQ( x.get1(),  1 );
  ASSERT_EQ( x.get2(), -2 );
  ASSERT_EQ( x.get3(),  3 );
  ASSERT_EQ( x.get4(), -5 );
  ASSERT_EQ( x.get1(), y.get1() );
  ASSERT_EQ( x.get2(), y.get2() );
  ASSERT_EQ( x.get3(), y.get3() );
  ASSERT_EQ( x.get4(), y.get4() );
  ASSERT_EQ( x.get1(), z.get1() );
  ASSERT_EQ( x.get2(), z.get2() );
  ASSERT_EQ( x.get3(), z.get3() );
  ASSERT_EQ( x.get4(), z.get4() );
}

//-----------------------------------------------------------------------------
//! \brief Test clear operation
//-----------------------------------------------------------------------------
TEST( TestSuiteR4, clear  )
{
  oquad x( 1, -2, 3, -5 );
  x.clear();
  ASSERT_EQ( x.get1(),  0 );
  ASSERT_EQ( x.get2(),  0 );
  ASSERT_EQ( x.get3(),  0 );
  ASSERT_EQ( x.get4(),  0 );
}

//-----------------------------------------------------------------------------
//! \brief Test R4 equals operations
//-----------------------------------------------------------------------------
TEST( TestSuiteR4, equal )
{
  vectR4 y;
  vectR4 x( 1, -2, 3, -5 );
  y=x;
  ASSERT_EQ( y.get1(), x.get1() );
  ASSERT_EQ( y.get2(), x.get2() );
  ASSERT_EQ( y.get3(), x.get3() );
  ASSERT_EQ( y.get4(), x.get4() );
}

//-----------------------------------------------------------------------------
//! \brief Test R4 equals operations
//-----------------------------------------------------------------------------
TEST( TestSuiteR4, put )
{
  vectR4 z;
  const vectR4 x( 1, -2, 3, -5 );
  const vectR4 y( 1, -2, 3, -5 );
  z.put(y);
  ASSERT_EQ( z.get1(), x.get1() );
  ASSERT_EQ( z.get2(), x.get2() );
  ASSERT_EQ( z.get3(), x.get3() );
  ASSERT_EQ( z.get4(), x.get4() );
}

//-----------------------------------------------------------------------------
//! \brief Test R4 equal
//-----------------------------------------------------------------------------
TEST( TestSuiteR4, equal2 )
{
  vectR4 x;
  x.clear();
  ASSERT_EQ( x.get1(), 0 );
  ASSERT_EQ( x.get2(), 0 );
  ASSERT_EQ( x.get3(), 0 );
  ASSERT_EQ( x.get4(), 0 );
}

//-----------------------------------------------------------------------------
//! \brief Test R4 put operation
//-----------------------------------------------------------------------------
TEST( TestSuiteR4, put2 )
{
  vectR4 x,y;
  x.put( 1, 2, 3, 4 );
  ASSERT_EQ( x.get1(), 1 );
  ASSERT_EQ( x.get2(), 2 );
  ASSERT_EQ( x.get3(), 3 );
  ASSERT_EQ( x.get4(), 4 );

  y.put( 1, 1./2., 1./3., 1./4. );
  ASSERT_EQ( y.get1(), 1./1. );
  ASSERT_EQ( y.get2(), 1./2. );
  ASSERT_EQ( y.get3(), 1./3. );
  ASSERT_EQ( y.get4(), 1./4. );
}

//-----------------------------------------------------------------------------
//! \brief Test R4 equals operations
//-----------------------------------------------------------------------------
TEST( TestSuiteR4, min )
{
  vectR4 x;
  x.put( -1   ,   -2  , -3   , -4     ); ASSERT_EQ( x.min(),  -4   );
  x.put( -1   ,   -2  , -13.0, -4     ); ASSERT_EQ( x.min(), -13   );
  x.put( -11.0,   -2.0, -10.0, -4.1   ); ASSERT_EQ( x.min(), -11   );
  x.put( -11.0,   -2.0, -10.0, -4.1   ); ASSERT_EQ( x.min(), -11.0 );
  x.put( -11.0,   -2.0, -10.0, -411.1 ); ASSERT_EQ( x.min(),-411.1 );
  x.put( -11.0, -512.0, -10.0,  811.1 ); ASSERT_EQ( x.min(),-512.0 );
}

//-----------------------------------------------------------------------------
//! \brief Test R4 equals operations
//-----------------------------------------------------------------------------
TEST( TestSuiteR4, max )
{
  vectR4 x;
  x.put( 1   ,   2  , 3   , 5     ); ASSERT_EQ( x.max(),   5   );
  x.put( 1   ,   2  , 13.0, 5     ); ASSERT_EQ( x.max(),  13   );
  x.put( 11.0,   2.0, 10.0, 4.1   ); ASSERT_EQ( x.max(),  11   );
  x.put( 11.0,   2.0, 10.0, 4.1   ); ASSERT_EQ( x.max(),  11.0 );
  x.put( 11.0,   2.0, 10.0, 411.1 ); ASSERT_EQ( x.max(), 411.1 );
  x.put( 11.0, 512.0, 10.0,-811.1 ); ASSERT_EQ( x.max(), 512.0 );
}

//-----------------------------------------------------------------------------
//! \brief Test R4 equals operations
//-----------------------------------------------------------------------------
TEST( TestSuiteR4, mag )
{
  const vectR4 x( sqrt(2), sqrt(3), sqrt(5), sqrt(7) );
  const vectR4 y(-sqrt(2),-sqrt(3),-sqrt(5),-sqrt(7) );
  ASSERT_DOUBLE_EQ( x.mag()*x.mag(), 2+3+5+7 );
  ASSERT_DOUBLE_EQ( y.mag()*y.mag(), 2+3+5+7 );
}

//-----------------------------------------------------------------------------
//! \brief Test R4 equals operations
//-----------------------------------------------------------------------------
TEST( TestSuiteR4, mpy )
{
  vectR4 x( sqrt(2.0*2.0/2.0), sqrt(3.0*3.0/2.0), sqrt(5.0*5.0/2.0), sqrt(7.0*7.0/2.0) );
  vectR4 y = x.mpy(sqrt(2));
  x.list();
  y.list();
  ASSERT_DOUBLE_EQ( y.get(0),  2.0 );
  ASSERT_DOUBLE_EQ( y.get(1),  3.0 );
  ASSERT_DOUBLE_EQ( y.get(2),  5.0 );
  ASSERT_DOUBLE_EQ( y.get(3),  7.0 );
}

//-----------------------------------------------------------------------------
//! \brief Test R4 add operation
//-----------------------------------------------------------------------------
TEST( TestSuiteR4, add )
{
  vectR4 x,y,z;
  x.put( 1, 2    , 3    , 4     );
  y.put( 1, 1./2., 1./3., 1./4. );
  z = x + y;
  ASSERT_EQ( z.get1(), 1 + 1./1. );
  ASSERT_EQ( z.get2(), 2 + 1./2. );
  ASSERT_EQ( z.get3(), 3 + 1./3. );
  ASSERT_EQ( z.get4(), 4 + 1./4. );
}

//-----------------------------------------------------------------------------
//! \brief Test R4 subtract operation
//-----------------------------------------------------------------------------
TEST( TestSuiteR4, sub )
{
  vectR4 x,y,v;
  x.put(1,2,3,4);
  y.put(1,1./2.,1./3.,1./4.);
  v = x + y - y;
  ASSERT_EQ( v.get1(), 1 + 1./1. - 1./1. );
  ASSERT_EQ( v.get2(), 2 + 1./2. - 1./2. );
  ASSERT_EQ( v.get3(), 3 + 1./3. - 1./3. );
  ASSERT_EQ( v.get4(), 4 + 1./4. - 1./4. );
}

//-----------------------------------------------------------------------------
//! \brief Test R4 negation operation
//-----------------------------------------------------------------------------
TEST( TestSuiteR4, neg )
{
  vectR4 x,y;
  x.put(1,1./2.,1./3.,1./4.);
  y = -x;
  ASSERT_EQ( y.get1(), -1./1. );
  ASSERT_EQ( y.get2(), -1./2. );
  ASSERT_EQ( y.get3(), -1./3. );
  ASSERT_EQ( y.get4(), -1./4. );
}

//-----------------------------------------------------------------------------
//! \brief Test R4 plus equals operation
//-----------------------------------------------------------------------------
TEST( TestSuiteR4, plusequal )
{
  vectR4 x,y,v;
  x.put( 1, 2    , 3    , 4     );
  y.put( 1, 1./2., 1./3., 1./4. );
  v  = x;
  v += y;
  v -= y;
  v += y;
  ASSERT_EQ( v.get1(), 1 + 1./1. - 1./1. + 1./1. );
  ASSERT_EQ( v.get2(), 2 + 1./2. - 1./2. + 1./2. );
  ASSERT_EQ( v.get3(), 3 + 1./3. - 1./3. + 1./3. );
  ASSERT_EQ( v.get4(), 4 + 1./4. - 1./4. + 1./4. );
}

//-----------------------------------------------------------------------------
//! \brief Test R4 subtract equals operation
//-----------------------------------------------------------------------------
TEST( TestSuiteR4, subequal )
{
  const vectR4 x( 1, 2    , 3    , 4     );
  const vectR4 y( 1, 1./2., 1./3., 1./4. );
  vectR4 v;
  v  = x;
  v -= y;
  ASSERT_EQ( v.get1(), 1 - 1./1. );
  ASSERT_EQ( v.get2(), 2 - 1./2. );
  ASSERT_EQ( v.get3(), 3 - 1./3. );
  ASSERT_EQ( v.get4(), 4 - 1./4. );
}

//-----------------------------------------------------------------------------
//! \brief Test R4 subtract equals operations
//-----------------------------------------------------------------------------
TEST( TestSuiteR4, equal4 )
{
  vectR4 x,v,y;
  x.put( 1, 2, 3, 4 );
  y.put(1,1./2.,1./3.,1./4. );
  v = x;
  v += y;
  v -= y;
  v += y;
  v -= y;
  ASSERT_EQ( v.get1(), 1 + 1./1. - 1./1. + 1./1.  - 1./1. );
  ASSERT_EQ( v.get2(), 2 + 1./2. - 1./2. + 1./2.  - 1./2. );
  ASSERT_EQ( v.get3(), 3 + 1./3. - 1./3. + 1./3.  - 1./3. );
  ASSERT_EQ( v.get4(), 4 + 1./4. - 1./4. + 1./4.  - 1./4. );
}

//-----------------------------------------------------------------------------
//! \brief Test R4 multiplication operation
//-----------------------------------------------------------------------------
TEST( TestSuiteR4, mult )
{
  vectR4 x,v,y;
  x.put( 1, 2    , 3    , 4     );
  y.put( 1, 1./2., 1./3., 1./4. );

  v = x;
  v += y;
  v -= y;
  v += y;
  v -= y;
  y = 2*v;

  ASSERT_EQ( y.get1(),  2 * ( 1 + 1./1. - 1./1. + 1./1.  - 1./1. ) );
  ASSERT_EQ( y.get2(),  2 * ( 2 + 1./2. - 1./2. + 1./2.  - 1./2. ) );
  ASSERT_EQ( y.get3(),  2 * ( 3 + 1./3. - 1./3. + 1./3.  - 1./3. ) );
  ASSERT_EQ( y.get4(),  2 * ( 4 + 1./4. - 1./4. + 1./4.  - 1./4. ) );

  v*=-M_PI;
  ASSERT_DOUBLE_EQ( v.get1(), -M_PI * ( 1 + 1./1. - 1./1. + 1./1.  - 1./1. ) );
  ASSERT_DOUBLE_EQ( v.get2(), -M_PI * ( 2 + 1./2. - 1./2. + 1./2.  - 1./2. ) );
  ASSERT_DOUBLE_EQ( v.get3(), -M_PI * ( 3 + 1./3. - 1./3. + 1./3.  - 1./3. ) );
  ASSERT_DOUBLE_EQ( v.get4(), -M_PI * ( 4 + 1./4. - 1./4. + 1./4.  - 1./4. ) );
}

//-----------------------------------------------------------------------------
//! \brief Test R4 dot product operation
//-----------------------------------------------------------------------------
TEST( TestSuiteR4, dotprod )
{
  vectR4 x;
  vectR4 y;
  double dotProd;

  x.put( 1, 2,     3,     4     );
  y.put( 1, 1./2., 1./3., 1./4. );
  dotProd = x^y;
  ASSERT_DOUBLE_EQ( dotProd, 1.0 + 1.0 + 1.0 + 1.0 );

  x.put( sqrt(2), sqrt(3), sqrt(5), sqrt(7) );
  y.put( sqrt(2), sqrt(3), sqrt(5), sqrt(7) );
  dotProd = x^y;
  ASSERT_DOUBLE_EQ( dotProd, 2 + 3 + 5 + 7 );

  x.put( sqrt(2), sqrt(3), sqrt(5), sqrt(7) );
  y.put(-sqrt(2),-sqrt(3),-sqrt(5),-sqrt(7) );
  dotProd = x^y;
  ASSERT_DOUBLE_EQ( dotProd, -(2 + 3 + 5 + 7 ) );
}

//-----------------------------------------------------------------------------
//! \brief Test otriple angle operation
//-----------------------------------------------------------------------------
TEST( TestSuiteR4, angle )
{
  vectR4 x, y;
  x.put( 1, 0, 0, 0 ); y.put( 0, 3, 0, 0 ); ASSERT_DOUBLE_EQ( pqtheta(x,y)/M_PI*180., 90.0 );
  x.put( 1, 0, 0, 0 ); y.put( 0, 0, 5, 0 ); ASSERT_DOUBLE_EQ( pqtheta(x,y)/M_PI*180., 90.0 );
  x.put( 0, 3, 0, 0 ); y.put( 1, 0, 0, 0 ); ASSERT_DOUBLE_EQ( pqtheta(x,y)/M_PI*180., 90.0 );
  x.put( 0, 3, 0, 0 ); y.put( 0, 0, 7, 0 ); ASSERT_DOUBLE_EQ( pqtheta(x,y)/M_PI*180., 90.0 );
  x.put( 0, 0, 7, 0 ); y.put( 2, 0, 0, 0 ); ASSERT_DOUBLE_EQ( pqtheta(x,y)/M_PI*180., 90.0 );
  x.put( 0, 0, 7, 0 ); y.put( 0, 5, 0, 0 ); ASSERT_DOUBLE_EQ( pqtheta(x,y)/M_PI*180., 90.0 );
  x.put( 1, 2, 0, 0 ); y.put(-2, 1, 0, 0 ); ASSERT_DOUBLE_EQ( pqtheta(x,y)/M_PI*180., 90.0 );
  x.put( 0,-3,-5, 0 ); y.put( 0, 5,-3, 0 ); ASSERT_DOUBLE_EQ( pqtheta(x,y)/M_PI*180., 90.0 );
}
