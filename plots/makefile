#-----------------------------------------------------------------------------
# Project Makefile
# Daniel J. Greenhoe
#-----------------------------------------------------------------------------
#--------------------------------------
# Paths
#--------------------------------------
#DIRPLOTS    = ../plots
DIRPLOTS    = .
DIRSTY      = ../../common/sty

#--------------------------------------
# Files
#--------------------------------------
GRP_DIE =                                                \
  $(DIRPLOTS)/diedftC1_1525_1200m120.pdf                 \
  $(DIRPLOTS)/diedftC1_1620_12000m120.pdf                \
  $(DIRPLOTS)/diedftC1_1620_12000m1200.pdf               \
  $(DIRPLOTS)/diedftR1_1525_1200m120.pdf                 \
  $(DIRPLOTS)/diedftR1_1620_12000m120.pdf                \
  $(DIRPLOTS)/diedftR1_1620_12000m1200.pdf               \
  $(DIRPLOTS)/diedftR6_1525_1200m120.pdf                 \
  $(DIRPLOTS)/diedftR6_1620_12000m120.pdf                \
  $(DIRPLOTS)/diedftR6_1620_12000m1200.pdf               \
  $(DIRPLOTS)/diehaarC1_12000m4000_h200_1050_D10.pdf     \
  $(DIRPLOTS)/diehaarR1_12000m4000_h200_1050_D10.pdf     \
  $(DIRPLOTS)/diehaarR6_12000m4000_h200_1050_D10.pdf     \

GRP_FDIE =                                               \
  $(DIRPLOTS)/fdie_5eed_16002_auto.pdf                   \
  $(DIRPLOTS)/fdie_5eed_16002_histo.pdf                  \
  $(DIRPLOTS)/fdie_5eed_51_seq.pdf                       \
  $(DIRPLOTS)/fdie_lp_12000_R1_hann16_euclid_auto.pdf    \
  $(DIRPLOTS)/fdie_lp_12000_R1_hann16_euclid_histo.pdf   \
  $(DIRPLOTS)/fdie_lp_12000_R1_hann16_euclid_seq.pdf     \
  $(DIRPLOTS)/fdie_lp_12000_R1_hann50_euclid_auto.pdf    \
  $(DIRPLOTS)/fdie_lp_12000_R1_hann50_euclid_histo.pdf   \
  $(DIRPLOTS)/fdie_lp_12000_R1_hann50_euclid_seq.pdf     \
  $(DIRPLOTS)/fdie_lp_12000_R1_rect16_euclid_auto.pdf    \
  $(DIRPLOTS)/fdie_lp_12000_R1_rect16_euclid_histo.pdf   \
  $(DIRPLOTS)/fdie_lp_12000_R1_rect16_euclid_seq.pdf     \
  $(DIRPLOTS)/fdie_lp_12000_R1_rect50_euclid_auto.pdf    \
  $(DIRPLOTS)/fdie_lp_12000_R1_rect50_euclid_histo.pdf   \
  $(DIRPLOTS)/fdie_lp_12000_R1_rect50_euclid_seq.pdf     \
  $(DIRPLOTS)/fdie_lp_12000_R6_hann16_euclid_auto.pdf    \
  $(DIRPLOTS)/fdie_lp_12000_R6_hann16_euclid_histo.pdf   \
  $(DIRPLOTS)/fdie_lp_12000_R6_hann16_euclid_seq.pdf     \
  $(DIRPLOTS)/fdie_lp_12000_R6_hann16_larc_auto.pdf      \
  $(DIRPLOTS)/fdie_lp_12000_R6_hann16_larc_histo.pdf     \
  $(DIRPLOTS)/fdie_lp_12000_R6_hann16_larc_seq.pdf       \
  $(DIRPLOTS)/fdie_lp_12000_R6_hann50_euclid_auto.pdf    \
  $(DIRPLOTS)/fdie_lp_12000_R6_hann50_euclid_histo.pdf   \
  $(DIRPLOTS)/fdie_lp_12000_R6_hann50_euclid_seq.pdf     \
  $(DIRPLOTS)/fdie_lp_12000_R6_hann50_larc_auto.pdf      \
  $(DIRPLOTS)/fdie_lp_12000_R6_hann50_larc_histo.pdf     \
  $(DIRPLOTS)/fdie_lp_12000_R6_hann50_larc_seq.pdf       \
  $(DIRPLOTS)/fdie_lp_12000_R6_rect16_euclid_auto.pdf    \
  $(DIRPLOTS)/fdie_lp_12000_R6_rect16_euclid_histo.pdf   \
  $(DIRPLOTS)/fdie_lp_12000_R6_rect16_euclid_seq.pdf     \
  $(DIRPLOTS)/fdie_lp_12000_R6_rect16_larc_auto.pdf      \
  $(DIRPLOTS)/fdie_lp_12000_R6_rect16_larc_histo.pdf     \
  $(DIRPLOTS)/fdie_lp_12000_R6_rect16_larc_seq.pdf       \
  $(DIRPLOTS)/fdie_lp_12000_R6_rect50_euclid_auto.pdf    \
  $(DIRPLOTS)/fdie_lp_12000_R6_rect50_euclid_histo.pdf   \
  $(DIRPLOTS)/fdie_lp_12000_R6_rect50_euclid_seq.pdf     \
  $(DIRPLOTS)/fdie_lp_12000_R6_rect50_larc_auto.pdf      \
  $(DIRPLOTS)/fdie_lp_12000_R6_rect50_larc_histo.pdf     \
  $(DIRPLOTS)/fdie_lp_12000_R6_rect50_larc_seq.pdf

GRP_RDIE =                                               \
  $(DIRPLOTS)/rdie_5eed_16002_auto.pdf                   \
  $(DIRPLOTS)/rdie_5eed_16002_histo.pdf                  \
  $(DIRPLOTS)/rdie_5eed_51_seq.pdf                       \
  $(DIRPLOTS)/rdie_lp_12000_R1_hann16_euclid_auto.pdf    \
  $(DIRPLOTS)/rdie_lp_12000_R1_hann16_euclid_histo.pdf   \
  $(DIRPLOTS)/rdie_lp_12000_R1_hann16_euclid_seq.pdf     \
  $(DIRPLOTS)/rdie_lp_12000_R1_hann50_euclid_auto.pdf    \
  $(DIRPLOTS)/rdie_lp_12000_R1_hann50_euclid_histo.pdf   \
  $(DIRPLOTS)/rdie_lp_12000_R1_hann50_euclid_seq.pdf     \
  $(DIRPLOTS)/rdie_lp_12000_R1_rect16_euclid_auto.pdf    \
  $(DIRPLOTS)/rdie_lp_12000_R1_rect16_euclid_histo.pdf   \
  $(DIRPLOTS)/rdie_lp_12000_R1_rect16_euclid_seq.pdf     \
  $(DIRPLOTS)/rdie_lp_12000_R1_rect50_euclid_auto.pdf    \
  $(DIRPLOTS)/rdie_lp_12000_R1_rect50_euclid_histo.pdf   \
  $(DIRPLOTS)/rdie_lp_12000_R1_rect50_euclid_seq.pdf     \
  $(DIRPLOTS)/rdie_lp_12000_R3_hann16_euclid_auto.pdf    \
  $(DIRPLOTS)/rdie_lp_12000_R3_hann16_euclid_histo.pdf   \
  $(DIRPLOTS)/rdie_lp_12000_R3_hann16_euclid_seq.pdf     \
  $(DIRPLOTS)/rdie_lp_12000_R3_hann16_larc_auto.pdf      \
  $(DIRPLOTS)/rdie_lp_12000_R3_hann16_larc_histo.pdf     \
  $(DIRPLOTS)/rdie_lp_12000_R3_hann16_larc_seq.pdf       \
  $(DIRPLOTS)/rdie_lp_12000_R3_hann50_euclid_auto.pdf    \
  $(DIRPLOTS)/rdie_lp_12000_R3_hann50_euclid_histo.pdf   \
  $(DIRPLOTS)/rdie_lp_12000_R3_hann50_euclid_seq.pdf     \
  $(DIRPLOTS)/rdie_lp_12000_R3_hann50_larc_auto.pdf      \
  $(DIRPLOTS)/rdie_lp_12000_R3_hann50_larc_histo.pdf     \
  $(DIRPLOTS)/rdie_lp_12000_R3_hann50_larc_seq.pdf       \
  $(DIRPLOTS)/rdie_lp_12000_R3_rect16_euclid_auto.pdf    \
  $(DIRPLOTS)/rdie_lp_12000_R3_rect16_euclid_histo.pdf   \
  $(DIRPLOTS)/rdie_lp_12000_R3_rect16_euclid_seq.pdf     \
  $(DIRPLOTS)/rdie_lp_12000_R3_rect16_larc_auto.pdf      \
  $(DIRPLOTS)/rdie_lp_12000_R3_rect16_larc_histo.pdf     \
  $(DIRPLOTS)/rdie_lp_12000_R3_rect16_larc_seq.pdf       \
  $(DIRPLOTS)/rdie_lp_12000_R3_rect50_euclid_auto.pdf    \
  $(DIRPLOTS)/rdie_lp_12000_R3_rect50_euclid_histo.pdf   \
  $(DIRPLOTS)/rdie_lp_12000_R3_rect50_euclid_seq.pdf     \
  $(DIRPLOTS)/rdie_lp_12000_R3_rect50_larc_auto.pdf      \
  $(DIRPLOTS)/rdie_lp_12000_R3_rect50_larc_histo.pdf     \
  $(DIRPLOTS)/rdie_lp_12000_R3_rect50_larc_seq.pdf

GRP_WDIE =                                               \
  $(DIRPLOTS)/wdie_5eed_16002_auto.pdf                   \
  $(DIRPLOTS)/wdie_5eed_16002_histo.pdf                  \
  $(DIRPLOTS)/wdie_5eed_51_seq.pdf                       \
  $(DIRPLOTS)/wdie_hp_1200_R1_hann16_euclid_auto.pdf     \
  $(DIRPLOTS)/wdie_hp_1200_R1_hann16_euclid_histo.pdf    \
  $(DIRPLOTS)/wdie_hp_1200_R1_hann16_euclid_seq.pdf      \
  $(DIRPLOTS)/wdie_hp_1200_R1_hann50_euclid_auto.pdf     \
  $(DIRPLOTS)/wdie_hp_1200_R1_hann50_euclid_histo.pdf    \
  $(DIRPLOTS)/wdie_hp_1200_R1_hann50_euclid_seq.pdf      \
  $(DIRPLOTS)/wdie_hp_1200_R1_rect16_euclid_auto.pdf     \
  $(DIRPLOTS)/wdie_hp_1200_R1_rect16_euclid_histo.pdf    \
  $(DIRPLOTS)/wdie_hp_1200_R1_rect16_euclid_seq.pdf      \
  $(DIRPLOTS)/wdie_hp_1200_R1_rect50_euclid_auto.pdf     \
  $(DIRPLOTS)/wdie_hp_1200_R1_rect50_euclid_histo.pdf    \
  $(DIRPLOTS)/wdie_hp_1200_R1_rect50_euclid_seq.pdf      \
  $(DIRPLOTS)/wdie_hp_1200_R6_hann16_euclid_auto.pdf     \
  $(DIRPLOTS)/wdie_hp_1200_R6_hann16_euclid_histo.pdf    \
  $(DIRPLOTS)/wdie_hp_1200_R6_hann16_euclid_seq.pdf      \
  $(DIRPLOTS)/wdie_hp_1200_R6_hann16_larc_auto.pdf       \
  $(DIRPLOTS)/wdie_hp_1200_R6_hann16_larc_histo.pdf      \
  $(DIRPLOTS)/wdie_hp_1200_R6_hann16_larc_seq.pdf        \
  $(DIRPLOTS)/wdie_hp_1200_R6_hann50_euclid_auto.pdf     \
  $(DIRPLOTS)/wdie_hp_1200_R6_hann50_euclid_histo.pdf    \
  $(DIRPLOTS)/wdie_hp_1200_R6_hann50_euclid_seq.pdf      \
  $(DIRPLOTS)/wdie_hp_1200_R6_hann50_larc_auto.pdf       \
  $(DIRPLOTS)/wdie_hp_1200_R6_hann50_larc_histo.pdf      \
  $(DIRPLOTS)/wdie_hp_1200_R6_hann50_larc_seq.pdf        \
  $(DIRPLOTS)/wdie_hp_1200_R6_rect16_euclid_auto.pdf     \
  $(DIRPLOTS)/wdie_hp_1200_R6_rect16_euclid_histo.pdf    \
  $(DIRPLOTS)/wdie_hp_1200_R6_rect16_euclid_seq.pdf      \
  $(DIRPLOTS)/wdie_hp_1200_R6_rect16_larc_auto.pdf       \
  $(DIRPLOTS)/wdie_hp_1200_R6_rect16_larc_histo.pdf      \
  $(DIRPLOTS)/wdie_hp_1200_R6_rect16_larc_seq.pdf        \
  $(DIRPLOTS)/wdie_hp_1200_R6_rect50_euclid_auto.pdf     \
  $(DIRPLOTS)/wdie_hp_1200_R6_rect50_euclid_histo.pdf    \
  $(DIRPLOTS)/wdie_hp_1200_R6_rect50_euclid_seq.pdf      \
  $(DIRPLOTS)/wdie_hp_1200_R6_rect50_larc_auto.pdf       \
  $(DIRPLOTS)/wdie_hp_1200_R6_rect50_larc_histo.pdf      \
  $(DIRPLOTS)/wdie_hp_1200_R6_rect50_larc_seq.pdf

GRP_WRDIE =                                              \
  $(DIRPLOTS)/wrdie_5eed_16002_auto.pdf                  \
  $(DIRPLOTS)/wrdie_5eed_16002_histo.pdf                 \
  $(DIRPLOTS)/wrdie_5eed_51_seq.pdf                      \
  $(DIRPLOTS)/wrdie_hp_1200_R1_hann50_euclid_auto.pdf    \
  $(DIRPLOTS)/wrdie_hp_1200_R1_hann50_euclid_histo.pdf   \
  $(DIRPLOTS)/wrdie_hp_1200_R1_hann50_euclid_seq.pdf     \
  $(DIRPLOTS)/wrdie_hp_1200_R1_rect50_euclid_auto.pdf    \
  $(DIRPLOTS)/wrdie_hp_1200_R1_rect50_euclid_histo.pdf   \
  $(DIRPLOTS)/wrdie_hp_1200_R1_rect50_euclid_seq.pdf     \
  $(DIRPLOTS)/wrdie_hp_1200_R3_hann50_euclid_auto.pdf    \
  $(DIRPLOTS)/wrdie_hp_1200_R3_hann50_euclid_histo.pdf   \
  $(DIRPLOTS)/wrdie_hp_1200_R3_hann50_euclid_seq.pdf     \
  $(DIRPLOTS)/wrdie_hp_1200_R3_hann50_larc_auto.pdf      \
  $(DIRPLOTS)/wrdie_hp_1200_R3_hann50_larc_histo.pdf     \
  $(DIRPLOTS)/wrdie_hp_1200_R3_hann50_larc_seq.pdf       \
  $(DIRPLOTS)/wrdie_hp_1200_R3_rect50_euclid_auto.pdf    \
  $(DIRPLOTS)/wrdie_hp_1200_R3_rect50_euclid_histo.pdf   \
  $(DIRPLOTS)/wrdie_hp_1200_R3_rect50_euclid_seq.pdf     \
  $(DIRPLOTS)/wrdie_hp_1200_R3_rect50_larc_auto.pdf      \
  $(DIRPLOTS)/wrdie_hp_1200_R3_rect50_larc_histo.pdf     \
  $(DIRPLOTS)/wrdie_hp_1200_R3_rect50_larc_seq.pdf

GRP_SPIN =                                               \
  $(DIRPLOTS)/spin_5eed_16002_auto.pdf                   \
  $(DIRPLOTS)/spin_5eed_16002_histo.pdf                  \
  $(DIRPLOTS)/spin_5eed_51_seq.pdf                       \
  $(DIRPLOTS)/spin_lp_12000_R1_hann16_euclid_auto.pdf    \
  $(DIRPLOTS)/spin_lp_12000_R1_hann16_euclid_histo.pdf   \
  $(DIRPLOTS)/spin_lp_12000_R1_hann16_euclid_seq.pdf     \
  $(DIRPLOTS)/spin_lp_12000_R1_hann50_euclid_auto.pdf    \
  $(DIRPLOTS)/spin_lp_12000_R1_hann50_euclid_histo.pdf   \
  $(DIRPLOTS)/spin_lp_12000_R1_hann50_euclid_seq.pdf     \
  $(DIRPLOTS)/spin_lp_12000_R1_rect16_euclid_auto.pdf    \
  $(DIRPLOTS)/spin_lp_12000_R1_rect16_euclid_histo.pdf   \
  $(DIRPLOTS)/spin_lp_12000_R1_rect16_euclid_seq.pdf     \
  $(DIRPLOTS)/spin_lp_12000_R1_rect50_euclid_auto.pdf    \
  $(DIRPLOTS)/spin_lp_12000_R1_rect50_euclid_histo.pdf   \
  $(DIRPLOTS)/spin_lp_12000_R1_rect50_euclid_seq.pdf     \
  $(DIRPLOTS)/spin_lp_12000_R2_hann16_euclid_auto.pdf    \
  $(DIRPLOTS)/spin_lp_12000_R2_hann16_euclid_histo.pdf   \
  $(DIRPLOTS)/spin_lp_12000_R2_hann16_euclid_seq.pdf     \
  $(DIRPLOTS)/spin_lp_12000_R2_hann16_larc_auto.pdf      \
  $(DIRPLOTS)/spin_lp_12000_R2_hann16_larc_histo.pdf     \
  $(DIRPLOTS)/spin_lp_12000_R2_hann16_larc_seq.pdf       \
  $(DIRPLOTS)/spin_lp_12000_R2_hann50_euclid_auto.pdf    \
  $(DIRPLOTS)/spin_lp_12000_R2_hann50_euclid_histo.pdf   \
  $(DIRPLOTS)/spin_lp_12000_R2_hann50_euclid_seq.pdf     \
  $(DIRPLOTS)/spin_lp_12000_R2_hann50_larc_auto.pdf      \
  $(DIRPLOTS)/spin_lp_12000_R2_hann50_larc_histo.pdf     \
  $(DIRPLOTS)/spin_lp_12000_R2_hann50_larc_seq.pdf       \
  $(DIRPLOTS)/spin_lp_12000_R2_rect16_euclid_auto.pdf    \
  $(DIRPLOTS)/spin_lp_12000_R2_rect16_euclid_histo.pdf   \
  $(DIRPLOTS)/spin_lp_12000_R2_rect16_euclid_seq.pdf     \
  $(DIRPLOTS)/spin_lp_12000_R2_rect16_larc_auto.pdf      \
  $(DIRPLOTS)/spin_lp_12000_R2_rect16_larc_histo.pdf     \
  $(DIRPLOTS)/spin_lp_12000_R2_rect16_larc_seq.pdf       \
  $(DIRPLOTS)/spin_lp_12000_R2_rect50_euclid_auto.pdf    \
  $(DIRPLOTS)/spin_lp_12000_R2_rect50_euclid_histo.pdf   \
  $(DIRPLOTS)/spin_lp_12000_R2_rect50_euclid_seq.pdf     \
  $(DIRPLOTS)/spin_lp_12000_R2_rect50_larc_auto.pdf      \
  $(DIRPLOTS)/spin_lp_12000_R2_rect50_larc_histo.pdf     \
  $(DIRPLOTS)/spin_lp_12000_R2_rect50_larc_seq.pdf       \
  $(DIRPLOTS)/wspin_5eed_16002_auto.pdf                  \
  $(DIRPLOTS)/wspin_5eed_16002_histo.pdf                 \
  $(DIRPLOTS)/wspin_5eed_51_seq.pdf                      \
  $(DIRPLOTS)/wspin_hp_1200_R1_hann50_euclid_auto.pdf    \
  $(DIRPLOTS)/wspin_hp_1200_R1_hann50_euclid_histo.pdf   \
  $(DIRPLOTS)/wspin_hp_1200_R1_hann50_euclid_seq.pdf     \
  $(DIRPLOTS)/wspin_hp_1200_R1_rect50_euclid_auto.pdf    \
  $(DIRPLOTS)/wspin_hp_1200_R1_rect50_euclid_histo.pdf   \
  $(DIRPLOTS)/wspin_hp_1200_R1_rect50_euclid_seq.pdf     \
  $(DIRPLOTS)/wspin_hp_1200_R2_hann50_euclid_auto.pdf    \
  $(DIRPLOTS)/wspin_hp_1200_R2_hann50_euclid_histo.pdf   \
  $(DIRPLOTS)/wspin_hp_1200_R2_hann50_euclid_seq.pdf     \
  $(DIRPLOTS)/wspin_hp_1200_R2_hann50_larc_auto.pdf      \
  $(DIRPLOTS)/wspin_hp_1200_R2_hann50_larc_histo.pdf     \
  $(DIRPLOTS)/wspin_hp_1200_R2_hann50_larc_seq.pdf       \
  $(DIRPLOTS)/wspin_hp_1200_R2_rect50_euclid_auto.pdf    \
  $(DIRPLOTS)/wspin_hp_1200_R2_rect50_euclid_histo.pdf   \
  $(DIRPLOTS)/wspin_hp_1200_R2_rect50_euclid_seq.pdf     \
  $(DIRPLOTS)/wspin_hp_1200_R2_rect50_larc_auto.pdf      \
  $(DIRPLOTS)/wspin_hp_1200_R2_rect50_larc_histo.pdf     \
  $(DIRPLOTS)/wspin_hp_1200_R2_rect50_larc_seq.pdf

GRP_DNA =                                                \
  $(DIRPLOTS)/dnadftC1.pdf                               \
  $(DIRPLOTS)/dnadftR1.pdf                               \
  $(DIRPLOTS)/dnadftR4.pdf                               \
  $(DIRPLOTS)/dnahaarC1_12000m4000_h200_1749_D5.pdf      \
  $(DIRPLOTS)/dnahaarR1_12000m4000_h200_1749_D5.pdf      \
  $(DIRPLOTS)/dnahaarR4_12000m4000_h200_1749_D5.pdf      \
  $(DIRPLOTS)/dna_5eed_16000_auto.pdf                    \
  $(DIRPLOTS)/dna_5eed_16000_histo.pdf                   \
  $(DIRPLOTS)/dna_5eed_51_seq.pdf                        \
  $(DIRPLOTS)/dna_AF086833_ebola_dftC1.pdf               \
  $(DIRPLOTS)/dna_AF086833_ebola_dftR1.pdf               \
  $(DIRPLOTS)/dna_AF086833_ebola_dftR4.pdf               \
  $(DIRPLOTS)/dna_ebola_51_seq.pdf                       \
  $(DIRPLOTS)/dna_ebola_auto.pdf                         \
  $(DIRPLOTS)/dna_ebola_histo.pdf                        \
  $(DIRPLOTS)/dna_mpbacterium_51_seq.pdf                 \
  $(DIRPLOTS)/dna_mpbacterium_auto.pdf                   \
  $(DIRPLOTS)/dna_mpbacterium_histo.pdf                  \
  $(DIRPLOTS)/dna_NC001416_phagelambdaC1_h1600.pdf       \
  $(DIRPLOTS)/dna_NC001416_phagelambdaC1_h4000.pdf       \
  $(DIRPLOTS)/dna_NC001416_phagelambdaR1b_hs1600.pdf     \
  $(DIRPLOTS)/dna_NC001416_phagelambdaR1_h1600.pdf       \
  $(DIRPLOTS)/dna_NC001416_phagelambdaR1_h4000.pdf       \
  $(DIRPLOTS)/dna_NC001416_phagelambdaR4_h1600.pdf       \
  $(DIRPLOTS)/dna_NC001416_phagelambdaR4_h4000.pdf       \
  $(DIRPLOTS)/dna_papaya1446_51_seq.pdf                  \
  $(DIRPLOTS)/dna_papaya1446_auto.pdf                    \
  $(DIRPLOTS)/dna_papaya1446_histo.pdf                   \
  $(DIRPLOTS)/dna_SARSCoV2_51_seq.pdf                    \
  $(DIRPLOTS)/dna_SARSCoV2_auto.pdf                      \
  $(DIRPLOTS)/dna_SARSCoV2_dftC1.pdf                     \
  $(DIRPLOTS)/dna_SARSCoV2_dftR1.pdf                     \
  $(DIRPLOTS)/dna_SARSCoV2_dftR4.pdf                     \
  $(DIRPLOTS)/dna_SARSCoV2_histo.pdf                     \
  $(DIRPLOTS)/dna_sars_51_seq.pdf                        \
  $(DIRPLOTS)/dna_sars_auto.pdf                          \
  $(DIRPLOTS)/dna_sars_histo.pdf

FILE_STY  =                                \
	$(DIRSTY)/wavelets.sty                   \
	$(DIRSTY)/math.sty                       \
	$(DIRSTY)/dan.sty                        \
	$(DIRSTY)/defaults.sty                   \
	$(DIRSTY)/fonts.sty                      \
	$(DIRSTY)/packages.sty                   \
	$(DIRSTY)/xsd.sty                        \
	$(DIRSTY)/language.sty

#--------------------------------------
# Programs
#--------------------------------------
PRG_TYPESET = xelatex.exe

#--------------------------------------
# master build rule
#--------------------------------------
top: $(GRP_DIE) $(GRP_FDIE) $(GRP_RDIE)  $(GRP_WDIE) $(GRP_WRDIE) $(GRP_SPIN) $(GRP_DNA)
	@echo Done ... graphics pdfs are ready!

#--------------------------------------
# implicit build control
#--------------------------------------
.tex.pdf:
	$(PRG_TYPESET) $(DIRPLOTS)/$&.tex

#--------------------------------------
# commands
#--------------------------------------
clean:
  del *.aux
  del *.log
  del *.out

scrub:
  make clean
  del *.pdf


