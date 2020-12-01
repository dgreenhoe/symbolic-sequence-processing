//=============================================================================
// Daniel J. Greenhoe
// DNA analysis
//=============================================================================
#include <stdio.h>
#include <math.h>
#include "lab2015ssp.h"
#include "gtest/gtest.h"

//---------------------------------------------------------------------------
//! \brief Make plot TeX files
//! \details These files can be compiled using xelatex to make pdf files 
//!          for inclusion into main pdf file for paper
//---------------------------------------------------------------------------
TEST( LabOCS, die )
{
  lab_fdie_ocs (     0x5EED, 16002,                "plots\\fdie"    ); // Example 3.4  (fair die sequence), 3 plots
  lab_rdie_ocs (     0x5EED, 16002,                "plots\\rdie"    ); // Example 3.5  (real die sequence), 3 plots
  lab_spin_ocs (     0x5EED, 16002,                "plots\\spin"    ); // Example 3.6  (spinner  sequence), 3 plots
  lab_wrdie_ocs(     0x5EED, 16002,                "plots\\wrdie"   ); // Example 3.7  (weighted real die sequence), 3 plots
  lab_wdie_ocs (     0x5EED, 16002,                "plots\\wdie"    ); // Example 3.8  (weighted die sequence), 3 plots
  lab_wspin_ocs(     0x5EED, 16002,                "plots\\wspin"   ); // Example 3.9  (weighted spinner sequence), 3 plots
  lab_rdie_lp(       0x5EED, 12000,   16,          "plots\\rdie_lp" ); // Example 4.1  (low pass filtering real die sequence)
  lab_rdie_lp(       0x5EED, 12000,   50,          "plots\\rdie_lp" ); // Example 4.1  (low pass filtering real die sequence)
  lab_spin_lp(       0x5EED, 12000,   16,          "plots\\spin_lp" ); // Example 4.2  (low pass filtering spinner sequence)
  lab_spin_lp(       0x5EED, 12000,   50,          "plots\\spin_lp" ); // Example 4.2  (low pass filtering spinner sequence)
  lab_fdie_lp(       0x5EED, 12000,   16,          "plots\\fdie_lp" ); // Example 4.3  (low pass filtering fair die sequence)
  lab_fdie_lp(       0x5EED, 12000,   50,          "plots\\fdie_lp" ); // Example 4.3  (low pass filtering fair die sequence)
  lab_wrdie_hp(      0x5EED, 1200 ,   50,          "plots\\wrdie_hp"); // Example 4.4  (high pass filtering weighted real die sequence)
  lab_wspin_hp(      0x5EED, 1200 ,   50,          "plots\\wspin_hp"); // Example 4.5  (high pass filtering weighted spinner sequence) 
  lab_wdie_hp(       0x5EED, 1200 ,   16,          "plots\\wdie_hp" ); // Example 4.6  (high pass filtering weighted die sequence)     
  lab_wdie_hp(       0x5EED, 1200 ,   50,          "plots\\wdie_hp" ); // Example 4.6  (high pass filtering weighted die sequence)     
  lab_die_nonstat34( 0x5eed, 1200 ,  120,  15,     "plots\\diedft"  ); // Example 4.7  (DFT analyis of die sequence)
  lab_die_nonstat34( 0x5eed, 12000, 1200,  16,     "plots\\diedft"  ); // Example 4.8  (DFT analyis of die sequence)
  lab_die_nonstat34( 0x5eed, 12000,  120,  16,     "plots\\diedft"  ); // Example 4.9  (DFT analyis of die sequence)
  lab_die_edge(      0x5eed, 12000, 4000, 200, 10, "plots\\diehaar" ); // Example 4.12 (Wavelet analysis of die sequence)
}
