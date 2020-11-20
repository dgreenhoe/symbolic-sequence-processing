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
//! \brief Test R6 set
//-----------------------------------------------------------------------------
TEST( TestSuiteR6, set )
{
  const vectR6 x(1,-2,3,-5,7,-11);
  ASSERT_EQ( x.get1(),  1 );
  ASSERT_EQ( x.get2(), -2 );
  ASSERT_EQ( x.get3(),  3 );
  ASSERT_EQ( x.get4(), -5 );
  ASSERT_EQ( x.get5(),  7 );
  ASSERT_EQ( x.get6(),-11 );
}

//-----------------------------------------------------------------------------
// \brief Test osix Clear operation
//-----------------------------------------------------------------------------
TEST( TestSuiteR6, seto )
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
TEST( TestSuiteR6, clear  )
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
//! \brief Test R6 equals operations
//-----------------------------------------------------------------------------
TEST( TestSuiteR6, equal )
{
  vectR6 y;
  vectR6 x(1,-2,3,-5,7,-11);
  y=x;
  ASSERT_EQ( y.get1(), x.get1() );
  ASSERT_EQ( y.get2(), x.get2() );
  ASSERT_EQ( y.get3(), x.get3() );
  ASSERT_EQ( y.get4(), x.get4() );
  ASSERT_EQ( y.get5(), x.get5() );
  ASSERT_EQ( y.get6(), x.get6() );
}

//-----------------------------------------------------------------------------
//! \brief Test R6 equals operations
//-----------------------------------------------------------------------------
TEST( TestSuiteR6, put )
{
  vectR6 z;
  const vectR6 x(1,-2,3,-5,7,-11);
  const vectR6 y(1,-2,3,-5,7,-11);
  z.put(y);
  ASSERT_EQ( z.get1(), x.get1() );
  ASSERT_EQ( z.get2(), x.get2() );
  ASSERT_EQ( z.get3(), x.get3() );
  ASSERT_EQ( z.get4(), x.get4() );
  ASSERT_EQ( z.get5(), x.get5() );
  ASSERT_EQ( z.get6(), x.get6() );
}

//-----------------------------------------------------------------------------
//! \brief Test R6 equal
//-----------------------------------------------------------------------------
TEST( TestSuiteR6, equal2 )
{
  vectR6 x;
  x.clear();
  ASSERT_EQ( x.get1(), 0 );
  ASSERT_EQ( x.get2(), 0 );
  ASSERT_EQ( x.get3(), 0 );
  ASSERT_EQ( x.get4(), 0 );
  ASSERT_EQ( x.get5(), 0 );
  ASSERT_EQ( x.get6(), 0 );
}

//-----------------------------------------------------------------------------
//! \brief Test R6 put operation
//-----------------------------------------------------------------------------
TEST( TestSuiteR6, put2 )
{
  vectR6 x,y;
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
}

//-----------------------------------------------------------------------------
//! \brief Test R6 equals operations
//-----------------------------------------------------------------------------
TEST( TestSuiteR6, min )
{
  vectR6 x;
  x.put( -1   ,   -2  , -3   , -4    , -5  , -6    ); ASSERT_EQ( x.min(),  -6   );
  x.put( -1   ,   -2  , -13.0, -4    , -5  , -6    ); ASSERT_EQ( x.min(), -13   );
  x.put( -11.0,   -2.0, -10.0, -4.1  , -5.2, -6.3  ); ASSERT_EQ( x.min(), -11   );
  x.put( -11.0,   -2.0, -10.0, -4.1  , -5.2, -67.3 ); ASSERT_EQ( x.min(), -67.3 );
  x.put( -11.0,   -2.0, -10.0, -411.1, -5.2, -67.3 ); ASSERT_EQ( x.min(),-411.1 );
  x.put( -11.0, -512.0, -10.0,  811.1, -5.2, -67.3 ); ASSERT_EQ( x.min(),-512.0 );
}

//-----------------------------------------------------------------------------
//! \brief Test R6 equals operations
//-----------------------------------------------------------------------------
TEST( TestSuiteR6, max )
{
  vectR6 x;
  x.put( 1   ,   2  , 3   , 4    , 5  , 6    ); ASSERT_EQ( x.max(),   6   );
  x.put( 1   ,   2  , 13.0, 4    , 5  , 6    ); ASSERT_EQ( x.max(),  13   );
  x.put( 11.0,   2.0, 10.0, 4.1  , 5.2, 6.3  ); ASSERT_EQ( x.max(),  11   );
  x.put( 11.0,   2.0, 10.0, 4.1  , 5.2, 67.3 ); ASSERT_EQ( x.max(),  67.3 );
  x.put( 11.0,   2.0, 10.0, 411.1, 5.2, 67.3 ); ASSERT_EQ( x.max(), 411.1 );
  x.put( 11.0, 512.0, 10.0,-811.1, 5.2, 67.3 ); ASSERT_EQ( x.max(), 512.0 );
}

//-----------------------------------------------------------------------------
//! \brief Test R6 equals operations
//-----------------------------------------------------------------------------
TEST( TestSuiteR6, mag )
{
  const vectR6 x( sqrt(2), sqrt(3), sqrt(5), sqrt(7), sqrt(11), sqrt(13) ); 
  const vectR6 y(-sqrt(2),-sqrt(3),-sqrt(5),-sqrt(7),-sqrt(11),-sqrt(13) ); 
  ASSERT_DOUBLE_EQ( x.mag()*x.mag(), 2+3+5+7+11+13 );
  ASSERT_DOUBLE_EQ( y.mag()*y.mag(), 2+3+5+7+11+13 );
}

//-----------------------------------------------------------------------------
//! \brief Test R6 equals operations
//-----------------------------------------------------------------------------
TEST( TestSuiteR6, mpy )
{
  vectR6 x( sqrt(2.0*2.0/2.0), sqrt(3.0*3.0/2.0), sqrt(5.0*5.0/2.0), sqrt(7.0*7.0/2.0), sqrt(11.0*11.0/2.0), sqrt(13.0*13.0/2.0) ); 
  vectR6 y = x.mpy(sqrt(2));
  x.list();
  y.list();
  ASSERT_DOUBLE_EQ( y.get(0),  2.0 );
  ASSERT_DOUBLE_EQ( y.get(1),  3.0 );
  ASSERT_DOUBLE_EQ( y.get(2),  5.0 );
  ASSERT_DOUBLE_EQ( y.get(3),  7.0 );
  ASSERT_DOUBLE_EQ( y.get(4), 11.0 );
  ASSERT_DOUBLE_EQ( y.get(5), 13.0 );
}

//-----------------------------------------------------------------------------
//! \brief Test R6 add operation
//-----------------------------------------------------------------------------
TEST( TestSuiteR6, add )
{
  vectR6 x,y,z;
  x.put( 1, 2    , 3    , 4     , 5    , 6     ); 
  y.put( 1, 1./2., 1./3., 1./4. , 1./5., 1./6. );
  z = x + y;
  ASSERT_EQ( z.get1(), 1 + 1./1. );
  ASSERT_EQ( z.get2(), 2 + 1./2. );
  ASSERT_EQ( z.get3(), 3 + 1./3. );
  ASSERT_EQ( z.get4(), 4 + 1./4. );
  ASSERT_EQ( z.get5(), 5 + 1./5. );
  ASSERT_EQ( z.get6(), 6 + 1./6. );
}

//-----------------------------------------------------------------------------
//! \brief Test R6 subtract operation
//-----------------------------------------------------------------------------
TEST( TestSuiteR6, sub )
{
  vectR6 x,y,v;
  x.put(1,2,3,4,5,6); 
  y.put(1,1./2.,1./3.,1./4.,1./5.,1./6.);
  v = x + y - y;
  ASSERT_EQ( v.get1(), 1 + 1./1. - 1./1. );
  ASSERT_EQ( v.get2(), 2 + 1./2. - 1./2. );
  ASSERT_EQ( v.get3(), 3 + 1./3. - 1./3. );
  ASSERT_EQ( v.get4(), 4 + 1./4. - 1./4. );
  ASSERT_EQ( v.get5(), 5 + 1./5. - 1./5. );
  ASSERT_EQ( v.get6(), 6 + 1./6. - 1./6. );
}

//-----------------------------------------------------------------------------
//! \brief Test R6 negation operation
//-----------------------------------------------------------------------------
TEST( TestSuiteR6, neg )
{
  vectR6 x,y;
  x.put(1,1./2.,1./3.,1./4.,1./5.,1./6.);
  y = -x;
  ASSERT_EQ( y.get1(), -1./1. );
  ASSERT_EQ( y.get2(), -1./2. );
  ASSERT_EQ( y.get3(), -1./3. );
  ASSERT_EQ( y.get4(), -1./4. );
  ASSERT_EQ( y.get5(), -1./5. );
  ASSERT_EQ( y.get6(), -1./6. );
}

//-----------------------------------------------------------------------------
//! \brief Test R6 plus equals operation
//-----------------------------------------------------------------------------
TEST( TestSuiteR6, plusequal )
{
  vectR6 x,y,v;
  x.put( 1, 2, 3, 4, 5, 6 ); 
  y.put( 1, 1./2., 1./3., 1./4., 1./5., 1./6. );
  v  = x;
  v += y;
  v -= y;
  v += y;
  ASSERT_EQ( v.get1(), 1 + 1./1. - 1./1. + 1./1. );
  ASSERT_EQ( v.get2(), 2 + 1./2. - 1./2. + 1./2. );
  ASSERT_EQ( v.get3(), 3 + 1./3. - 1./3. + 1./3. );
  ASSERT_EQ( v.get4(), 4 + 1./4. - 1./4. + 1./4. );
  ASSERT_EQ( v.get5(), 5 + 1./5. - 1./5. + 1./5. );
  ASSERT_EQ( v.get6(), 6 + 1./6. - 1./6. + 1./6. );
}

//-----------------------------------------------------------------------------
//! \brief Test R6 subtract equals operation
//-----------------------------------------------------------------------------
TEST( TestSuiteR6, subequal )
{
  const vectR6 x( 1, 2,     3,     4,     5,     6     ); 
  const vectR6 y( 1, 1./2., 1./3., 1./4., 1./5., 1./6. );
  vectR6 v;
  v  = x;
  v -= y;
  ASSERT_EQ( v.get1(), 1 - 1./1. );
  ASSERT_EQ( v.get2(), 2 - 1./2. );
  ASSERT_EQ( v.get3(), 3 - 1./3. );
  ASSERT_EQ( v.get4(), 4 - 1./4. );
  ASSERT_EQ( v.get5(), 5 - 1./5. );
  ASSERT_EQ( v.get6(), 6 - 1./6. );
}

//-----------------------------------------------------------------------------
//! \brief Test R6 subtract equals operations
//-----------------------------------------------------------------------------
TEST( TestSuiteR6, equal4 )
{
  vectR6 x,v,y;
  x.put( 1, 2, 3, 4, 5, 6 );
  y.put(1,1./2.,1./3.,1./4.,1./5.,1./6.);
  v = x;
  v += y;
  v -= y;
  v += y;
  v -= y;
  ASSERT_EQ( v.get1(), 1 + 1./1. - 1./1. + 1./1.  - 1./1. );
  ASSERT_EQ( v.get2(), 2 + 1./2. - 1./2. + 1./2.  - 1./2. );
  ASSERT_EQ( v.get3(), 3 + 1./3. - 1./3. + 1./3.  - 1./3. );
  ASSERT_EQ( v.get4(), 4 + 1./4. - 1./4. + 1./4.  - 1./4. );
  ASSERT_EQ( v.get5(), 5 + 1./5. - 1./5. + 1./5.  - 1./5. );
  ASSERT_EQ( v.get6(), 6 + 1./6. - 1./6. + 1./6.  - 1./6. );
}

//-----------------------------------------------------------------------------
//! \brief Test R6 multiplication operation
//-----------------------------------------------------------------------------
TEST( TestSuiteR6, mult )
{
  vectR6 x,v,y;
  x.put( 1, 2, 3, 4, 5, 6 );
  y.put(1,1./2.,1./3.,1./4.,1./5.,1./6.);

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
  ASSERT_EQ( y.get5(),  2 * ( 5 + 1./5. - 1./5. + 1./5.  - 1./5. ) );
  ASSERT_EQ( y.get6(),  2 * ( 6 + 1./6. - 1./6. + 1./6.  - 1./6. ) );

  v*=-M_PI;
  ASSERT_DOUBLE_EQ( v.get1(), -M_PI * ( 1 + 1./1. - 1./1. + 1./1.  - 1./1. ) );
  ASSERT_DOUBLE_EQ( v.get2(), -M_PI * ( 2 + 1./2. - 1./2. + 1./2.  - 1./2. ) );
  ASSERT_DOUBLE_EQ( v.get3(), -M_PI * ( 3 + 1./3. - 1./3. + 1./3.  - 1./3. ) );
  ASSERT_DOUBLE_EQ( v.get4(), -M_PI * ( 4 + 1./4. - 1./4. + 1./4.  - 1./4. ) );
  ASSERT_DOUBLE_EQ( v.get5(), -M_PI * ( 5 + 1./5. - 1./5. + 1./5.  - 1./5. ) );
  ASSERT_DOUBLE_EQ( v.get6(), -M_PI * ( 6 + 1./6. - 1./6. + 1./6.  - 1./6. ) );
}

//-----------------------------------------------------------------------------
//! \brief Test R6 dot product operation
//-----------------------------------------------------------------------------
TEST( TestSuiteR6, dotprod )
{
  vectR6 x; 
  vectR6 y;
  double dotProd;

  x.put( 1, 2,     3,     4,     5,     6     ); 
  y.put( 1, 1./2., 1./3., 1./4., 1./5., 1./6. );
  dotProd = x^y;
  ASSERT_DOUBLE_EQ( dotProd, 1.0 + 1.0 + 1.0 + 1.0 + 1.0 + 1.0 );
  
  x.put( sqrt(2), sqrt(3), sqrt(5), sqrt(7), sqrt(11), sqrt(13) ); 
  y.put( sqrt(2), sqrt(3), sqrt(5), sqrt(7), sqrt(11), sqrt(13) ); 
  dotProd = x^y;
  ASSERT_DOUBLE_EQ( dotProd, 2 + 3 + 5 + 7 + 11 + 13 );
  
  x.put( sqrt(2), sqrt(3), sqrt(5), sqrt(7), sqrt(11), sqrt(13) ); 
  y.put(-sqrt(2),-sqrt(3),-sqrt(5),-sqrt(7),-sqrt(11),-sqrt(13) ); 
  dotProd = x^y;
  ASSERT_DOUBLE_EQ( dotProd, -(2 + 3 + 5 + 7 + 11 + 13) );
}

//-----------------------------------------------------------------------------
//! \brief Test otriple angle operation
//-----------------------------------------------------------------------------
TEST( TestSuiteR6, angle )
{
  vectR6 x, y;
  x.put( 1, 0, 0, 0, 0, 0 ); y.put( 0, 3, 0, 0, 0, 0 ); ASSERT_DOUBLE_EQ( pqtheta(x,y)/M_PI*180., 90.0 );
  x.put( 1, 0, 0, 0, 0, 0 ); y.put( 0, 0, 5, 0, 0, 0 ); ASSERT_DOUBLE_EQ( pqtheta(x,y)/M_PI*180., 90.0 );
  x.put( 0, 3, 0, 0, 0, 0 ); y.put( 1, 0, 0, 0, 0, 0 ); ASSERT_DOUBLE_EQ( pqtheta(x,y)/M_PI*180., 90.0 );
  x.put( 0, 3, 0, 0, 0, 0 ); y.put( 0, 0, 7, 0, 0, 0 ); ASSERT_DOUBLE_EQ( pqtheta(x,y)/M_PI*180., 90.0 );
  x.put( 0, 0, 7, 0, 0, 0 ); y.put( 2, 0, 0, 0, 0, 0 ); ASSERT_DOUBLE_EQ( pqtheta(x,y)/M_PI*180., 90.0 );
  x.put( 0, 0, 7, 0, 0, 0 ); y.put( 0, 5, 0, 0, 0, 0 ); ASSERT_DOUBLE_EQ( pqtheta(x,y)/M_PI*180., 90.0 );
  x.put( 1, 2, 0, 0, 0, 0 ); y.put(-2, 1, 0, 0, 0, 0 ); ASSERT_DOUBLE_EQ( pqtheta(x,y)/M_PI*180., 90.0 );
  x.put( 0,-3,-5, 0, 0, 0 ); y.put( 0, 5,-3, 0, 0, 0 ); ASSERT_DOUBLE_EQ( pqtheta(x,y)/M_PI*180., 90.0 );
}
