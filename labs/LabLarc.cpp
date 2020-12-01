//=============================================================================
// Daniel J. Greenhoe
// Lagrange Arc Distance analysis
//=============================================================================
#include <math.h>
#include "lab2015larc.h"
#include "gtest/gtest.h"

//---------------------------------------------------------------------------
//! \brief Calculate data for the paper 2015sphx:
//! "An extension to the spherical metric using Lagrange interpolation"
//---------------------------------------------------------------------------
TEST( LabLarc, ballR2 )
{
  lab_larc_distances_R2("../data/lab_larc_distances_R2");     //Thm 3.13, Rem 3.14, Exm 3.15

  lab_larc_ball_R2(0,   0,   1, 0.5,  4, 0.001,  2000, 1000, "../data/larc_ball(0_0)");    //Example 3.17: larc in R^2
  lab_larc_ball_R2(0.5, 0.5, 1, 0.5,  4, 0.001,  2000, 1000, "../data/larc_ball(05_05)");  //Example 3.17: larc in R^2
  lab_larc_ball_R2(1,   1,   1, 0.5,  4, 0.001,  2000, 1000, "../data/larc_ball(1_1)");    //Example 3.17: larc in R^2
  lab_larc_ball_R2(2,   2,   1, 0.5,  4, 0.001,  8000, 1000, "../data/larc_ball(2_2)");    //Example 3.17: larc in R^2
  lab_larc_ball_R2(M_PI, M_PI,  1, 0.5,  4, 0.001,  2000, 1000, "../data/larc_ball(pi_pi)");  //Example 3.17: larc in R^2
  lab_larc_ball_R2(4,   4,   1, 0.5,  4, 0.001,  2000, 1000, "../data/larc_ball(4_4)");    //Example 3.17: larc in R^2
    //             |    |    |  |     |  |          |     |   |___base filename
    //             |    |    |  |     |  |          |     |_______number of plot data points
    //             |    |    |  |     |  |          |_____________number of search iterations per plot point
    //             |    |    |  |     |  |________________________maximum acceptable error
    //             |    |    |  |     |___________________________maximum search distance from point p
    //             |    |    |  |_________________________________minimum search distance from point p
    //             |    |    |____________________________________radius of ball
    //             |    |_________________________________________x of point p=(x,y)
    //             |______________________________________________y of point p=(x,y)
}

//---------------------------------------------------------------------------
//! \brief Calculate data for the paper 2015sphx:
//! "An extension to the spherical metric using Lagrange interpolation"
//---------------------------------------------------------------------------
TEST( LabLarc, ballR3 )
{
  lab_larc_distances_R3("data/lab_larc_distances_R3");     //Example 3.16: larc in R^3
    
  lab_larc_ball_R3(1, 1,  1, 1, "../data/larc_ball(1_1_1)");  //Example 3.18: larc in R^3
  lab_larc_ball_R3(1,-1,  1, 1, "../data/larc_ball(1_1_1)");  //Example 3.18: larc in R^3
  lab_larc_ball_R3(0, 0,  1, 1, "../data/larc_ball(0_0_1)");  //Example 3.18: larc in R^3
  lab_larc_ball_R3(0, 0,  0, 1, "../data/larc_ball(0_0_0)");  //Example 3.18: larc in R^3
  lab_larc_ball_R3(0, 0, -1, 1, "../data/larc_ball(0_0_-1)"); //Example 3.18: larc in R^3
  lab_larc_ball_R3(0, 0, -2, 1, "../data/larc_ball(0_0_-2)"); //Example 3.18: larc in R^3
  lab_larc_ball_R3(0, 0, -3, 1, "../data/larc_ball(0_0_-3)"); //Example 3.18: larc in R^3
  lab_larc_ball_R3(0, 0, -5, 1, "../data/larc_ball(0_0_-5)"); //Example 3.18: larc in R^3
  lab_larc_ball_R3(0, 0,-10, 1, "../data/larc_ball(0_0_-10)");//Example 3.18: larc in R^3
}

