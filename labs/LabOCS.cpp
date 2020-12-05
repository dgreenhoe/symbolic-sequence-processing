//=============================================================================
// Daniel J. Greenhoe
//=============================================================================
#include <stdio.h>
#include <math.h>
#include "lab2015ssp.h"
#include "gtest/gtest.h"

//---------------------------------------------------------------------------
//! \brief Fair die analysis
//---------------------------------------------------------------------------
TEST( LabOCS, fdie )
{
  const unsigned seed = 0x5EED;
  //ASSERT_EQ( lab_fdie_ocs ( seed, 16002,       "../plots/fdie"    ), 0 ); // Example 3.4  (fair die sequence), 3 plots
  //ASSERT_EQ( lab_fdie_lp(   seed, 12000,   16, "../plots/fdie_lp" ), 0 ); // Example 4.3  (low pass filtering fair die sequence)
  ASSERT_EQ( lab_fdie_lp(   seed, 120,   16, "../plots/fdie_lp_tmp" ), 0 ); // Example 4.3  (low pass filtering fair die sequence)
  //ASSERT_EQ( lab_fdie_lp(   seed, 12000,   50, "../plots/fdie_lp" ), 0 ); // Example 4.3  (low pass filtering fair die sequence)
}

//---------------------------------------------------------------------------
//! \brief Real die analysis
//---------------------------------------------------------------------------
TEST( LabOCS, rdie )
{
  const unsigned seed = 0x5EED;
  ASSERT_EQ( lab_rdie_ocs ( seed, 16002,     "../plots/rdie"    ), 0 ); // Example 3.5  (real die sequence), 3 plots
  ASSERT_EQ( lab_rdie_lp(   seed, 12000, 16, "../plots/rdie_lp" ), 0 ); // Example 4.1  (low pass filtering real die sequence)
  ASSERT_EQ( lab_rdie_lp(   seed, 12000, 50, "../plots/rdie_lp" ), 0 ); // Example 4.1  (low pass filtering real die sequence)
}

//---------------------------------------------------------------------------
//! \brief Weighted die analysis
//---------------------------------------------------------------------------
TEST( LabOCS, wdie )
{
  const unsigned seed = 0x5EED;
//  ASSERT_EQ( lab_wrdie_ocs( seed, 16002,       "../plots/wrdie"   ), 0 ); // Example 3.7  (weighted real die sequence), 3 plots
//  ASSERT_EQ( lab_wdie_ocs ( seed, 16002,       "../plots/wdie"    ), 0 ); // Example 3.8  (weighted die sequence), 3 plots
//  ASSERT_EQ( lab_wrdie_hp(  seed, 1200 ,   50, "../plots/wrdie_hp"), 0 ); // Example 4.4  (high pass filtering weighted real die sequence)
  ASSERT_EQ( lab_wdie_hp(   seed, 120 ,   16, 'r', "../plots/wdie_hp" ), 0 ); // Example 4.6  (high pass filtering weighted die sequence)     
//  ASSERT_EQ( lab_wdie_hp(   seed, 1200 ,   16, 'r', "../plots/wdie_hp" ), 0 ); // Example 4.6  (high pass filtering weighted die sequence)     
//  ASSERT_EQ( lab_wdie_hp(   seed, 1200 ,   50, 'r', "../plots/wdie_hp" ), 0 ); // Example 4.6  (high pass filtering weighted die sequence)     
}

//---------------------------------------------------------------------------
//! \brief Non-stationary die analysis
//---------------------------------------------------------------------------
TEST( LabOCS, ndie )
{
  const unsigned     seed = 0x5EED;
  ASSERT_EQ( lab_die_nonstat34( seed, 1200 ,  120,  15,     "../plots/diedft"  ), 0 ); // Example 4.7  (DFT analyis of die sequence)
  ASSERT_EQ( lab_die_nonstat34( seed, 12000, 1200,  16,     "../plots/diedft"  ), 0 ); // Example 4.8  (DFT analyis of die sequence)
  ASSERT_EQ( lab_die_nonstat34( seed, 12000,  120,  16,     "../plots/diedft"  ), 0 ); // Example 4.9  (DFT analyis of die sequence)
  ASSERT_EQ( lab_die_edge(      seed, 12000, 4000, 200, 10, "../plots/diehaar" ), 0 ); // Example 4.12 (Wavelet analysis of die sequence)
}

//---------------------------------------------------------------------------
//! \brief Spinner analysis
//---------------------------------------------------------------------------
TEST( LabOCS, spinner )
{
  const unsigned seed = 0x5EED;
  ASSERT_EQ( lab_spin_ocs ( seed, 16002,     "../plots/spin"    ), 0 ); // Example 3.6  (spinner  sequence), 3 plots
  ASSERT_EQ( lab_wspin_ocs( seed, 16002,     "../plots/wspin"   ), 0 ); // Example 3.9  (weighted spinner sequence), 3 plots
  ASSERT_EQ( lab_spin_lp(   seed, 12000, 16, "../plots/spin_lp" ), 0 ); // Example 4.2  (low pass filtering spinner sequence)
  ASSERT_EQ( lab_spin_lp(   seed, 12000, 50, "../plots/spin_lp" ), 0 ); // Example 4.2  (low pass filtering spinner sequence)
  ASSERT_EQ( lab_wspin_hp(  seed, 1200 , 50, "../plots/wspin_hp"), 0 ); // Example 4.5  (high pass filtering weighted spinner sequence) 
}
