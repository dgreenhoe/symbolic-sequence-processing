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
GRP_DIE =                                  \
  $(DIRPLOTS)/fdie_5eed_16002_histo.pdf    \
  $(DIRPLOTS)/fdie_5eed_51_seq.pdf         \
  $(DIRPLOTS)/fdie_5eed_16002_auto.pdf     \
  $(DIRPLOTS)/rdie_5eed_51_seq.pdf         \
  $(DIRPLOTS)/rdie_5eed_16002_histo.pdf    \
  $(DIRPLOTS)/rdie_5eed_16002_auto.pdf     \
  $(DIRPLOTS)/wrdie_5eed_51_seq.pdf        \
  $(DIRPLOTS)/wrdie_5eed_16002_histo.pdf   \
  $(DIRPLOTS)/wrdie_5eed_16002_auto.pdf    \
  $(DIRPLOTS)/wdie_5eed_51_seq.pdf         \
  $(DIRPLOTS)/wdie_5eed_16002_histo.pdf    \
  $(DIRPLOTS)/wdie_5eed_16002_auto.pdf     \
  $(DIRPLOTS)/wspin_5eed_51_seq.pdf        \
  $(DIRPLOTS)/wspin_5eed_16002_histo.pdf   \
  $(DIRPLOTS)/wspin_5eed_16002_auto.pdf

GRP_DNA =                                              \
  $(DIRPLOTS)/dna_5eed_51_seq.pdf                      \
  $(DIRPLOTS)/dna_5eed_16000_histo.pdf                 \
  $(DIRPLOTS)/dna_5eed_16000_auto.pdf                  \
  $(DIRPLOTS)/dnadftR1.pdf                             \
  $(DIRPLOTS)/dnadftC1.pdf                             \
  $(DIRPLOTS)/dnadftR4.pdf                             \
  $(DIRPLOTS)/dnahaarR1_12000m4000_h200_1749_D5.pdf    \
  $(DIRPLOTS)/dnahaarC1_12000m4000_h200_1749_D5.pdf    \
  $(DIRPLOTS)/dnahaarR4_12000m4000_h200_1749_D5.pdf    \
  $(DIRPLOTS)/dna_mpbacterium_51_seq.pdf               \
  $(DIRPLOTS)/dna_mpbacterium_histo.pdf                \
  $(DIRPLOTS)/dna_mpbacterium_auto.pdf                 \
  $(DIRPLOTS)/dna_NC001416_phagelambdaR1b_hs1600.pdf   \
  $(DIRPLOTS)/dna_NC001416_phagelambdaR1_h1600.pdf     \
  $(DIRPLOTS)/dna_NC001416_phagelambdaC1_h1600.pdf     \
  $(DIRPLOTS)/dna_NC001416_phagelambdaR4_h1600.pdf     \
  $(DIRPLOTS)/dna_NC001416_phagelambdaR1_h4000.pdf     \
  $(DIRPLOTS)/dna_NC001416_phagelambdaC1_h4000.pdf     \
  $(DIRPLOTS)/dna_NC001416_phagelambdaR4_h4000.pdf     \
  $(DIRPLOTS)/dna_ebola_histo.pdf                      \
  $(DIRPLOTS)/dna_ebola_51_seq.pdf                     \
  $(DIRPLOTS)/dna_ebola_auto.pdf                       \
  $(DIRPLOTS)/dna_AF086833_ebola_dftR1.pdf             \
  $(DIRPLOTS)/dna_AF086833_ebola_dftC1.pdf             \
  $(DIRPLOTS)/dna_AF086833_ebola_dftR4.pdf             \
  $(DIRPLOTS)/dna_sars_51_seq.pdf                      \
  $(DIRPLOTS)/dna_sars_histo.pdf                       \
  $(DIRPLOTS)/dna_sars_auto.pdf                        \
  $(DIRPLOTS)/dna_papaya1446_51_seq.pdf                \
  $(DIRPLOTS)/dna_papaya1446_histo.pdf                 \
  $(DIRPLOTS)/dna_papaya1446_auto.pdf                  \
  $(DIRPLOTS)/dna_SARSCoV2_histo.pdf                   \
  $(DIRPLOTS)/dna_SARSCoV2_51_seq.pdf                  \
  $(DIRPLOTS)/dna_SARSCoV2_auto.pdf                    \
  $(DIRPLOTS)/dna_SARSCoV2_dftR1.pdf                   \
  $(DIRPLOTS)/dna_SARSCoV2_dftC1.pdf                   \
  $(DIRPLOTS)/dna_SARSCoV2_dftR4.pdf

GRP_SPIN =                                 \
  $(DIRPLOTS)/spin_5eed_16002_auto.pdf     \
  $(DIRPLOTS)/spin_5eed_16002_histo.pdf    \
  $(DIRPLOTS)/spin_5eed_51_seq.pdf         \
  $(DIRPLOTS)/wspin_5eed_16002_auto.pdf    \
  $(DIRPLOTS)/wspin_5eed_16002_histo.pdf   \
  $(DIRPLOTS)/wspin_5eed_51_seq.pdf

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
top: $(GRP_DIE) $(GRP_DNA) $(GRP_SPIN) 
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

