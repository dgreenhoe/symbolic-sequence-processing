#-----------------------------------------------------------------------------
# Project Makefile
# Daniel J. Greenhoe
#-----------------------------------------------------------------------------
#--------------------------------------
# TARGET
#--------------------------------------
#TARGET = symseq
TARGET = ssp
#.RECIPEPREFIX=>

#--------------------------------------
# Files
#--------------------------------------
FILES_OBJ = main.obj c1.obj c4.obj c6.obj die.obj dna.obj dnan.obj dft.obj elliptic.obj euclid.obj fairdie.obj fileplot.obj lab2015larc.obj lab2015ssp.obj larc.obj mca.obj r1.obj r2.obj r3.obj r4.obj r6.obj r1op.obj r2op.obj r3op.obj r4op.obj r6op.obj realdie.obj spinner.obj symseq.obj test.obj
FILES_H   = main.h   c1.h   c4.h   c6.h   die.h   dna.h   dnan.h   dft.h   elliptic.h   euclid.h   fairdie.h   fileplot.h   lab2015larc.h   lab2015ssp.h   larc.h   mca.h   r1.h   r2.h   r3.h   r4.h   r6.h   r1op.h   r2op.h   r3op.h   r4op.h   r6op.h   realdie.h   spinner.h   symseq.h   test.h  

#FILES_OBJ = main.obj c1.obj euclid.obj fileplot.obj lab2015larc.obj larc.obj r1.obj r2.obj r3.obj r4.obj r6.obj r1op.obj r2op.obj r3op.obj r4op.obj r6op.obj testlarc.obj
#FILES_H   = main.h   c1.h   euclid.h   fileplot.h   lab2015larc.h   larc.h   r1.h   r2.h   r3.h   r4.h   r6.h   r1op.h   r2op.h   r3op.h   r4op.h   r6op.h   testlarc.h  

#--------------------------------------
# directories
#--------------------------------------
DIR_BIN = c:\p\bcc\bin      # binaries (execuatables) directory
DIR_INC = c:\p\bcc\Include  # include files directory
DIR_LIB = c:\p\bcc\Lib      # libraries directory

#--------------------------------------
# Programs
#--------------------------------------
PRG_COMPILE = $(DIR_BIN)\bcc32
PRG_LINK    = $(DIR_BIN)\bcc32

#--------------------------------------
# master build control
#--------------------------------------
$(TARGET).exe: $(FILES_OBJ) $(FILES_H) makefile
  time /T
  $(PRG_LINK)  -I$(DIR_INC) -I. -L$(DIR_LIB) -e$(TARGET).exe $(FILES_OBJ)
  dir $(TARGET).exe
#  @echo Done ... sto/graphics pdfs are ready!

#--------------------------------------
# implicit build control
#--------------------------------------
.cpp.obj:
	$(PRG_COMPILE) -I$(DIR_INC) -I. -L$(DIR_LIB) -c -w-ncf $&.cpp

#--------------------------------------
# explicit build control
#--------------------------------------
#lat5_m3_tall.pdf:  $(PATHLAT)/lat5_m3.tex shell_lat5_m3_tall.tex
#	$(PRG_TYPESET) shell_$&.tex
#	copy /V/Y shell_$&.pdf $&.pdf
#	del shell_$&.pdf

#--------------------------------------
# commands
#--------------------------------------
clean:
  del *.obj
  del *.tds
  del *.tex

scrub:
  make clean
  del *.exe

new:
  del *.obj
  make

zip:
  zip -o -9 -r $(TARGET).zip * -x *.obj *.tds *.exe *.dat *.sty *.zip *.bat *.pdf tmp\*
  unzip -l $(TARGET).zip
  dir *.zip
