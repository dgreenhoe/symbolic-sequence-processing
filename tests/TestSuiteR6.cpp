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
//! \brief Test R6 add operation
//-----------------------------------------------------------------------------
TEST( TestSuiteR6, add )
{
  vectR6 x,y,z;
  x.put(1,2,3,4,5,6); 
  y.put(1,1./2.,1./3.,1./4.,1./5.,1./6.);
  z=x+y;
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
TEST( TestSuiteR6, equal3 )
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

  v*=-3;
  ASSERT_EQ( v.get1(), -3 * ( 1 + 1./1. - 1./1. + 1./1.  - 1./1. ) );
  ASSERT_EQ( v.get2(), -3 * ( 2 + 1./2. - 1./2. + 1./2.  - 1./2. ) );
  ASSERT_EQ( v.get3(), -3 * ( 3 + 1./3. - 1./3. + 1./3.  - 1./3. ) );
  ASSERT_EQ( v.get4(), -3 * ( 4 + 1./4. - 1./4. + 1./4.  - 1./4. ) );
  ASSERT_EQ( v.get5(), -3 * ( 5 + 1./5. - 1./5. + 1./5.  - 1./5. ) );
  ASSERT_EQ( v.get6(), -3 * ( 6 + 1./6. - 1./6. + 1./6.  - 1./6. ) );
}

