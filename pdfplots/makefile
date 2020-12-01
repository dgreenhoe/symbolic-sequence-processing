#-----------------------------------------------------------------------------
# Project Makefile
# Daniel J. Greenhoe
#-----------------------------------------------------------------------------
#--------------------------------------
# Paths
#--------------------------------------
DIRPLOTS    = ../plots
DIRSTY      = ../../common/sty

#--------------------------------------
# Files
#--------------------------------------
GRP_DIE =                                  \
  $(DIRPLOTS)/fdie_5eed_16002_auto.pdf     \
  $(DIRPLOTS)/fdie_5eed_16002_histo.pdf    \
  $(DIRPLOTS)/fdie_5eed_51_seq.pdf         \
  $(DIRPLOTS)/rdie_5eed_16002_auto.pdf     \
  $(DIRPLOTS)/rdie_5eed_16002_histo.pdf    \
  $(DIRPLOTS)/rdie_5eed_51_seq.pdf         \
  $(DIRPLOTS)/wdie_5eed_16002_auto.pdf     \
  $(DIRPLOTS)/wdie_5eed_16002_histo.pdf    \
  $(DIRPLOTS)/wdie_5eed_51_seq.pdf         \
  $(DIRPLOTS)/wrdie_5eed_16002_auto.pdf    \
  $(DIRPLOTS)/wrdie_5eed_16002_histo.pdf   \
  $(DIRPLOTS)/wrdie_5eed_51_seq.pdf

GRP_DNA =                                  \
  $(DIRPLOTS)/dna_5eed_16000_auto.pdf      \
  $(DIRPLOTS)/dna_5eed_16000_histo.pdf     \
  $(DIRPLOTS)/dna_5eed_51_seq.pdf          \
  $(DIRPLOTS)/dna_papaya1446_51_seq.pdf    \
  $(DIRPLOTS)/dna_papaya1446_auto.pdf      \
  $(DIRPLOTS)/dna_papaya1446_histo.pdf     \
  $(DIRPLOTS)/dna_SARSCoV2_51_seq.pdf      \
  $(DIRPLOTS)/dna_SARSCoV2_auto.pdf        \
  $(DIRPLOTS)/dna_SARSCoV2_histo.pdf

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
{$(DIRPLOTS)}.tex.pdf:
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

