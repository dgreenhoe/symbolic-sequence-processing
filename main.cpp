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
  const int N = 8;
  const int M = 2*N-1;
  const int p = 3;
  const double x[p][N]  = { -1, -1,  2,  3,  0, -1, -1, -1,
                            -1, -1, -1, -1, -1,  2,  3,  0,
                            -1,  2,  3,  0, -1, -1, -1, -1  };
  const double y[p][N]  = { -1, -1,  2,  3,  0, -1, -1, -1,
                            -1, -1, -1, -1, -1,  2,  3,  0,
                            -1,  2,  3,  0, -1, -1, -1, -1  };
  const double z[p][N]  = { -1, -1, -1,  2,  3,  0, -1, -1, 
                             0, -1, -1, -1, -1, -1,  2,  3, 
                            -1, -1,  2,  3,  0, -1, -1, -1   };
  double xval;
  double Rxy[p][M];
  int row, col, index, m, m1, m2, m3, n1, n2, n3;

  printf("\nsequence x:");
  for(row=0; row<p; row++) {
    printf("\n");
    for(col=0; col<N; col++){
      printf("%2.0lf  ",x[row][col]);
    }}

  printf("\nsequence y:");
  for(row=0; row<p; row++) {
    printf("\n");
    for(col=0; col<N; col++){
      printf("%2.0lf  ",y[row][col]);
    }}

  for(row=0; row<p; row++){
    for(col=0; col<M; col++){
      Rxy[row][col] = 0;
    }}  

  for(m=0; m<M; m++){
    for(n1=0; n1<N; n1++) {
      for(n2=0; n2<N; n2++) {
        for(n3=0; n3<N; n3++) {
          index = m+n1-N;
          if(index<0)          xval = 0;
          else if(index>(N-1)) xval = 0;
          else                 xval = x[0][index];
          Rxy[0][m] += xval * y[0][n1]; /* dimension-1 */

          index = m+n2-N;
          if(index<0)          xval = 0;
          else if(index>(N-1)) xval = 0;
          else                 xval = x[0][index];
          Rxy[1][m] += xval *y[1][n2]; /* dimension-2 */

          index = m+n3-N;
          if(index<0)          xval = 0;
          else if(index>(N-1)) xval = 0;
          else                 xval = x[0][index];
          Rxy[2][m] += xval *y[2][n3]; /* dimension-3 */
  }}}}

  printf("\nestimated Rxy:");
  for(row=0; row<p; row++) {
    printf("\n");
    for(col=0; col<M; col++){
      printf("%5.0lf  ",Rxy[row][col]);
    }}
}

//int main(int argc, char *argv[])
//{
//  printf("This is opseq.exe by Daniel J. Greenhoe \n");
//  //printf("  for support of version 0.50 (2016 July 04 Monday) of the text \n");
//  //printf("  \"A book concerning symbolic sequence processing\" \n");
//  //printf("   by Daniel J. Greenhoe \n");
//
//  //bspline_Sdat();
//  testing::InitGoogleTest(&argc, argv);
//  return RUN_ALL_TESTS();
//}

