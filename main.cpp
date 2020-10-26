/*============================================================================
 * Daniel J. Greenhoe
 *============================================================================*/
/*=====================================
 * headers
 *=====================================*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "main.h"
#include "test.h"
#include "bsplines.h"
#include "lab2015ssp.h"
#include "lab2015larc.h"

/*-------------------------------------
 * prototypes
 *-------------------------------------*/
int perform_tests(void);
int test_2015larc(void);
int make_2015larc_data(void);
int make_2015ssp_texplots(void);

/*-------------------------------------
 * main
 *-------------------------------------*/
int main(int argc, char *argv[]){
  //printf("This is opseq.exe by Daniel J. Greenhoe \n");
  //printf("  for support of version 0.50 (2016 July 04 Monday) of the text \n");
  //printf("  \"A book concerning symbolic sequence processing\" \n");
  //printf("   by Daniel J. Greenhoe \n");

  bspline_Sdat();
  perform_tests();
  make_2015ssp_texplots();
  test_2015larc();
  make_2015larc_data();
  //lab_larc_pti(0.5, 0, -1, 1, 20, "larc_pti_05");
  //lab_larc_pti(0.5,   0, -1, 1, 2, "tmp");
  //lab_larc_pti(2,   0, -1, 1, 2, "tmp");
    //         |_____________________________p

  return 0;
  }

/*---------------------------------------------------------------------------
 * calculate data for the paper 2015sphx:
 * "An extension to the spherical metric using Lagrange interpolation"
 *---------------------------------------------------------------------------*/
int test_2015larc(void){
  test_opair();
  test_vectR2();
  test_complex();
  test_seqR2();
  test_otriple();
  test_osix();
  test_vectR6();
  test_conj();
  test_larc_metric_R2();
  test_larc_metric_R3();
  test_larc_metric_R6();
  return 0;
  }

/*---------------------------------------------------------------------------
 * calculate data for the paper 2015sphx:
 * "An extension to the spherical metric using Lagrange interpolation"
 *---------------------------------------------------------------------------*/
int make_2015larc_data(void){
  lab_larc_distances_R2("data/lab_larc_distances_R2");     //Thm 3.13, Rem 3.14, Exm 3.15
  lab_larc_distances_R3("data/lab_larc_distances_R3");     //Example 3.16: larc in R^3
    
  lab_larc_ball_R2(0,   0,   1, 0.5,  4, 0.001,  2000, 1000, "data/larc_ball(0_0)");    //Example 3.17: larc in R^2
  lab_larc_ball_R2(0.5, 0.5, 1, 0.5,  4, 0.001,  2000, 1000, "data/larc_ball(05_05)");  //Example 3.17: larc in R^2
  lab_larc_ball_R2(1,   1,   1, 0.5,  4, 0.001,  2000, 1000, "data/larc_ball(1_1)");    //Example 3.17: larc in R^2
  lab_larc_ball_R2(2,   2,   1, 0.5,  4, 0.001,  8000, 1000, "data/larc_ball(2_2)");    //Example 3.17: larc in R^2
  lab_larc_ball_R2(PI,  PI,  1, 0.5,  4, 0.001,  2000, 1000, "data/larc_ball(pi_pi)");  //Example 3.17: larc in R^2
  lab_larc_ball_R2(4,   4,   1, 0.5,  4, 0.001,  2000, 1000, "data/larc_ball(4_4)");    //Example 3.17: larc in R^2
    //             |    |    |  |     |  |          |     |   |___base filename
    //             |    |    |  |     |  |          |     |_______number of plot data points
    //             |    |    |  |     |  |          |_____________number of search iterations per plot point
    //             |    |    |  |     |  |________________________maximum acceptable error
    //             |    |    |  |     |___________________________maximum search distance from point p
    //             |    |    |  |_________________________________minimum search distance from point p
    //             |    |    |____________________________________radius of ball
    //             |    |_________________________________________x of point p=(x,y)
    //             |______________________________________________y of point p=(x,y)

  lab_larc_ball_R3(1, 1,  1, 1, "data/larc_ball(1_1_1)");  //Example 3.18: larc in R^3
  lab_larc_ball_R3(1,-1,  1, 1, "data/larc_ball(1_1_1)");  //Example 3.18: larc in R^3
  lab_larc_ball_R3(0, 0,  1, 1, "data/larc_ball(0_0_1)");  //Example 3.18: larc in R^3
  lab_larc_ball_R3(0, 0,  0, 1, "data/larc_ball(0_0_0)");  //Example 3.18: larc in R^3
  lab_larc_ball_R3(0, 0, -1, 1, "data/larc_ball(0_0_-1)"); //Example 3.18: larc in R^3
  lab_larc_ball_R3(0, 0, -2, 1, "data/larc_ball(0_0_-2)"); //Example 3.18: larc in R^3
  lab_larc_ball_R3(0, 0, -3, 1, "data/larc_ball(0_0_-3)"); //Example 3.18: larc in R^3
  lab_larc_ball_R3(0, 0, -5, 1, "data/larc_ball(0_0_-5)"); //Example 3.18: larc in R^3
  lab_larc_ball_R3(0, 0,-10, 1, "data/larc_ball(0_0_-10)");//Example 3.18: larc in R^3

  return 0;
  }

/*---------------------------------------------------------------------------
 * make plot TeX files
 * these files can be compiled using xelatex to make pdf files 
 * for inclusion into main pdf file for paper
 *---------------------------------------------------------------------------*/
int make_2015ssp_texplots(void){
  lab_fdie_ocs (0x5EED, 16002, "plots\\fdie"); // Example 3.4  (fair die sequence), 3 plots
  lab_rdie_ocs (0x5EED, 16002, "plots\\rdie"); // Example 3.5  (real die sequence), 3 plots
  lab_spin_ocs (0x5EED, 16002, "plots\\spin"); // Example 3.6  (spinner  sequence), 3 plots
  lab_wrdie_ocs(0x5EED, 16002, "plots\\wrdie");// Example 3.7  (weighted real die sequence), 3 plots
  lab_wdie_ocs (0x5EED, 16002, "plots\\wdie"); // Example 3.8  (weighted die sequence), 3 plots
  lab_wspin_ocs(0x5EED, 16002, "plots\\wspin");// Example 3.9  (weighted spinner sequence), 3 plots
  lab_dna_ocs  (0x5EED, 16000, "plots\\dna");  // Example 3.10 (artificial DNA  sequence), 3 plots
  lab_dna_ocs  ("..\\..\\common\\symseq\\fasta\\NC004718_sars.dat","plots\\dna_sars");  // Example 3.11 (SARS virus)
  lab_dna_ocs  ("..\\..\\common\\symseq\\fasta\\AF086833_ebola.dat","plots\\dna_ebola"); // Example 3.12 (Ebola virus)
  lab_dna_ocs  ("..\\..\\common\\symseq\\fasta\\gi880815890_MelissococcusPlutonius.dat","plots\\dna_mpbacterium"); // Example 3.13 (Bacterium)
  lab_dnan_ocs ("..\\..\\common\\symseq\\fasta\\gi187567196_papaya1446.dat","plots\\dna_papaya1446");  // Example 3.14 (papaya)
  lab_rdie_lp(0x5EED, 12000, 16, "plots\\rdie_lp");  // Example 4.1 (low pass filtering real die sequence)
  lab_rdie_lp(0x5EED, 12000, 50, "plots\\rdie_lp");  // Example 4.1 (low pass filtering real die sequence)
  lab_spin_lp(0x5EED, 12000, 16, "plots\\spin_lp");  // Example 4.2 (low pass filtering spinner sequence)
  lab_spin_lp(0x5EED, 12000, 50, "plots\\spin_lp");  // Example 4.2 (low pass filtering spinner sequence)
  lab_fdie_lp( 0x5EED, 12000, 16, "plots\\fdie_lp"); // Example 4.3 (low pass filtering fair die sequence)
  lab_fdie_lp( 0x5EED, 12000, 50, "plots\\fdie_lp"); // Example 4.3 (low pass filtering fair die sequence)
  lab_wrdie_hp(0x5EED, 1200, 50, "plots\\wrdie_hp"); // Example 4.4 (high pass filtering weighted real die sequence)
  lab_wspin_hp(0x5EED, 1200, 50, "plots\\wspin_hp"); // Example 4.5 (high pass filtering weighted spinner sequence) 
  lab_wdie_hp(0x5EED, 1200, 16, "plots\\wdie_hp");   // Example 4.6 (high pass filtering weighted die sequence)     
  lab_wdie_hp(0x5EED, 1200, 50, "plots\\wdie_hp");   // Example 4.6 (high pass filtering weighted die sequence)     
  lab_die_nonstat34(0x5eed, 1200,  120,  15, "plots\\diedft"); // Example 4.7 (DFT analyis of die sequence)
  lab_die_nonstat34(0x5eed, 12000, 1200, 16, "plots\\diedft"); //   Example 4.8 (DFT analyis of die sequence)
  lab_die_nonstat34(0x5eed, 12000, 120,  16, "plots\\diedft"); // Example 4.9 (DFT analyis of die sequence)
  lab_dna_nonstatCT(0x5eed, 12000, 1200, 24, -1, 5, "plots\\dnadft");// Example 4.10 (DFT analyis of DNA sequence)
  lab_dna_dft("..\\..\\common\\symseq\\fasta\\AF086833_ebola.dat","plots\\dna_AF086833_ebola_dft");// Example 4.11 (DFT of Ebola DNA sequence)
  lab_die_edge(0x5eed, 12000, 4000, 200, 10, "plots\\diehaar");// Example 4.12 (Wavelet analysis of die sequence)
  lab_dna_edge(0x5eed, 12000, 4000, 200, 17, "plots\\dnahaar");// Example 4.13 (Wavelet analysis of DNA sequence)
  lab_dna_averaging(1600,"..\\..\\common\\symseq\\fasta\\NC001416_phagelambda.dat","plots\\dna_NC001416_phagelambda");// Example 4.14 (sliding histogram of Phage Lambda)
  lab_dna_edge(1600,"..\\..\\common\\symseq\\fasta\\NC001416_phagelambda.dat","plots\\dna_NC001416_phagelambda");// Example 4.14 (sliding histogram of Phage Lambda)
  lab_dna_edge(4000,"..\\..\\common\\symseq\\fasta\\NC001416_phagelambda.dat","plots\\dna_NC001416_phagelambda");// Example 4.14 (sliding histogram of Phage Lambda)
  return 0;
  }

/*---------------------------------------------------------------------------
 * perform tests
 *---------------------------------------------------------------------------*/
int perform_tests(void){
  test_opair();
  test_vectR2();
  test_complex();
  test_seqR2();
  if(test_otriple() !=0) return -1;
  if(test_osix()    !=0) return -1;
  if(test_vectR6()  !=0) return -1;
  if(test_dieC1()  !=0) return -1;
  test_conj();
  if(test_dft_R1()!=0)            return -1;
  if(test_pqtheta()!=0)        return -1;
  if(test_larc_metric_R2()!=0) return -1;
  if(test_larc_metric_R3()!=0) return -1;
  if(test_larc_metric_R6()!=0) return -1;
  if(test_circle()        !=0) return -1;
  if(test_circle_d1()     !=0) return -1;
  if(test_ellipse_d1()    !=0) return -1;
  //test_rdie_metric();
  //test_correlation();
  if(test_halfcircle()    ==0) return -1;
  if(test_findt()         ==0) return -1;
  if(test_perimeter()     ==0) return -1;
  if(test_balloon_metric()==0) return -1;
  if(test_mca_metric()    ==0) return -1;
  if(test_spinner()    ==0) return -1;
  if(test_rdie()    ==0) return -1;
  if(test_dna_metric()==0) return -1;
  if(test_dnan_metric()==0) return -1;
  test_die();
  test_rdie();
  test_dft_R1();
  test_expi();
  return 0;
  }


