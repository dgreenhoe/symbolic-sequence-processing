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

TEST(TestSuiteB, TestB1)
{
  printf("TestSuiteB TestB1\n");
  ASSERT_EQ(cos(0), 1);
}

TEST(TestSuiteB, TestB2)
{
  printf("TestSuiteB TestB2\n");
  ASSERT_EQ(0, 1);
}


