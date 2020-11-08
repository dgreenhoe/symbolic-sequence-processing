/*============================================================================
 * Daniel J. Greenhoe
 *============================================================================*/
/*=====================================
 * headers
 *=====================================*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gtest/gtest.h"  // https://github.com/google/googletest/blob/master/googletest/docs/primer.md
//#include "test.h"

TEST(TestSuite, Test1)
{
  printf("TestSuite Test1\n");
  ASSERT_EQ(cos(0), 1);
}

TEST(TestSuite, Test2)
{
  printf("TestSuite Test2\n");
  ASSERT_EQ(sin(0), 0);
}


