//=============================================================================
// Daniel J. Greenhoe
// DNA analysis
//=============================================================================
#include <stdio.h>
#include <math.h>
#include "lab2015ssp.h"
#include "gtest/gtest.h"

//---------------------------------------------------------------------------
//! \brief Analysis of pseudo-randomly generated DNA sequence
//---------------------------------------------------------------------------
TEST( LabDNA, prng )
{
  const unsigned seed = 0x5EED;
  lab_dna_ocs  (     seed, 16000,                   "../plots/dna"     );
  lab_dna_nonstatCT( seed, 12000, 1200, 24 , -1, 5, "../plots/dnadft"  );
  lab_dna_edge(      seed, 12000, 4000, 200, 17,    "../plots/dnahaar" );
}

//---------------------------------------------------------------------------
//! \brief Melissococcus Plutonius bacterium DNA analysis
//---------------------------------------------------------------------------
TEST( LabDNA, mpbacterium )
{
  lab_dna_ocs  ( "../fasta/gi880815890_MelissococcusPlutonius.fasta", "../plots/dna_mpbacterium" );
}

//---------------------------------------------------------------------------
//! \brief Phage-Lambda DNA analysis
//! \details Sliding histogram of Phage Lambda
//---------------------------------------------------------------------------
TEST( LabDNA, phagelambda )
{
  lab_dna_averaging( 1600, "../fasta/NC001416_phagelambda.fasta", "../plots/dna_NC001416_phagelambda" );
  lab_dna_edge(      1600, "../fasta/NC001416_phagelambda.fasta", "../plots/dna_NC001416_phagelambda" );
  lab_dna_edge(      4000, "../fasta/NC001416_phagelambda.fasta", "../plots/dna_NC001416_phagelambda" );
}

//---------------------------------------------------------------------------
//! \brief Ebola DNA analysis
//---------------------------------------------------------------------------
TEST( LabDNA, ebola )
{
  lab_dna_ocs( "../fasta/AF086833_ebola.fasta", "../plots/dna_ebola"              );
  lab_dna_dft( "../fasta/AF086833_ebola.fasta", "../plots/dna_AF086833_ebola_dft" );
}

//---------------------------------------------------------------------------
//! \brief SARS DNA analysis
//---------------------------------------------------------------------------
TEST( LabDNA, sars )
{
  lab_dna_ocs( "../fasta/NC004718_sars.fasta", "../plots/dna_sars" );
}

//---------------------------------------------------------------------------
//! \brief Papaya DNA-N analysis
//---------------------------------------------------------------------------
TEST( LabDNA, papaya )
{
  lab_dnan_ocs( "../fasta/gi187567196_papaya1446.fasta", "../plots/dna_papaya1446" ); 
}

//-----------------------------------------------------------------------------
//! \brief SARS-CoV-2 DNA analysis
//-----------------------------------------------------------------------------
TEST( LabDNA, covid )
{
  lab_dna_ocs( "../fasta/GenBank_MT072688-1_SARS-COV-2.fasta", "../plots/dna_SARSCoV2"     );
  lab_dna_dft( "../fasta/GenBank_MT072688-1_SARS-COV-2.fasta", "../plots/dna_SARSCoV2_dft" ); 
}

