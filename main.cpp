//=============================================================================
// Daniel J. Greenhoe
//=============================================================================
//=====================================
// headers
//=====================================
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "main.h"
#include "bsplines.h"
#include "lab2015ssp.h"
#include "lab2015larc.h"
#include "gtest/gtest.h"  // https://github.com/google/googletest/blob/master/googletest/docs/primer.md

//-------------------------------------
// main
// https://stackoverflow.com/questions/12657596/
//-------------------------------------
int main(int argc, char *argv[])
{
  printf("This is opseq.exe by Daniel J. Greenhoe \n");
  //printf("  for support of version 0.50 (2016 July 04 Monday) of the text \n");
  //printf("  \"A book concerning symbolic sequence processing\" \n");
  //printf("   by Daniel J. Greenhoe \n");

  //bspline_Sdat();
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

