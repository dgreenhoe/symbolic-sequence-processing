//===============================================================================
// Daniel J. Greenhoe
//===============================================================================
//=======================================
// headers
//=======================================
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "main.h"
#include "symseq.h"
#include "r1.h"
#include "r2.h"
#include "r3.h"
#include "r4.h"
#include "r6.h"
#include "c1.h"
#include "c4.h"
#include "c6.h"
#include "r1op.h"
#include "r2op.h"
#include "r3op.h"
#include "r4op.h"
#include "r6op.h"
#include "elliptic.h"
#include "larc.h"
#include "mca.h"
#include "euclid.h"
#include "die.h"
#include "realdie.h"
#include "fairdie.h"
#include "spinner.h"
#include "dna.h"
#include "dnan.h"
#include "dft.h"
#include "fileplot.h"
#include "lab2015ssp.h"

//-----------------------------------------------------------------------------
//! \brief generate plot files for fair die sequence
//-----------------------------------------------------------------------------
int lab_fdie_ocs(const unsigned seed, const long N, const char *basefilename)
{
  time_t time1; time(&time1);   //starting time stamp (passed to plotting routine)
  fdieseq x(N);   //for die sequence
  seqR1 Rxx(2*N+1); //for auto-correlation of x

  char comment[2*1024];          //comment to be passed to plotting function
  char buf[2*1024];          //general purpose buffer
  char filename[1024];
  FILE *lptr; // pointer to log  file

  //----------------------------------------------
  //open log file
  //----------------------------------------------
  printf(         "----------------------------------------------------------------\n");
  sprintf(comment,"Experiment: calculate statistics of length %ld fair die sequence",N); printf("%s\n",comment);
  printf(         "----------------------------------------------------------------\n");
  sprintf(filename,"%s_ocs_%ld",basefilename,N);
  lptr=log_open (filename,time1,comment);

  //----------------------------------------------
  //sprintf(buf,"Generate fair die sequence...");printofe(lptr,buf,time1);
  //----------------------------------------------
  sprintf(buf,"length %ld uniform die sequence using seed value 0x%x:\n",N,seed);
  x.randomize(seed);
  x.list(0,299,buf ,"...\n",1,lptr);
  x.histogram(1,lptr);

  //----------------------------------------------
  //plot die sequence
  //----------------------------------------------
  sprintf(comment,"length %ld die sequence",N);
  sprintf(filename,"%s_%x_51_seq.tex",basefilename,seed);
  sprintf(buf,"generate plot file \"%s\"\n",filename); printof(lptr,buf);
  plot_ocs_seq((symseq *)&x, 0,50, time1, "die", filename, comment,lptr);

  //----------------------------------------------
  //plot histogram
  //----------------------------------------------
  sprintf(comment,"HISTOGRAM of length %ld real die sequence",N);
  sprintf(filename,"%s_%x_%ld_histo.tex",basefilename,seed,N);
  sprintf(buf,"generate plot file \"%s\"\n",filename); printof(lptr,buf);
  plot_ocs_histo((symseq *)&x, time1, "die", filename, comment,lptr);

  //----------------------------------------------
  sprintf(buf,"Perform auto-correlation operations...");printofe(lptr,buf,time1);
  //----------------------------------------------
  if(x.Rxxo(&Rxx,1)){fprintf(stderr,"ERROR computing auto-correlation.\n"); return -1;}
  // |      |   |____________switch to turn on counting display
  // |      |________________pointer to output correlation sequence
  // |_______________________input real die sequence
  sprintf(buf,"done");printofe(lptr,buf,time1);
  sprintf(buf,"auto-correlation sequence [r_n|n=N-50...N], where N=%ld:\n",N);
  Rxx.list(N-50,N,buf,"\n",1,lptr);
  sprintf(buf,"auto-correlation sequence [r_n|n=N...N+50], where N=%ld:\n",N);
  Rxx.list(N,N+50,buf,"\n",1,lptr);
  sprintf(comment,"AUTO-CORRELATION of length %ld real die sequence",N);
  sprintf(filename,"%s_%x_%ld_auto.tex",basefilename,seed,N);
  sprintf(buf,"generate plot file \"%s\"\n",filename); printof(lptr,buf);
  plot_ocs_auto(&Rxx, 1, time1, filename,comment,lptr);

  //----------------------------------------------
  // close log file
  //----------------------------------------------
  sprintf(buf,"lab_fdie_ocs experiment complete"); printofe(lptr,buf,time1);
  plot_close(lptr,time1);
  return 0;
}

//-----------------------------------------------------------------------------
//! \brief generate plot files for real die sequence
//-----------------------------------------------------------------------------
int lab_rdie_ocs(const unsigned seed, const long N, const char *basefilename){
  time_t time1; time(&time1);   //starting time stamp (passed to plotting routine)
  rdieseq x(N);   //for die sequence
  seqR1 Rxx(2*N+1); //for auto-correlation of x

  char comment[2*1024];          //comment to be passed to plotting function
  char buf[2*1024];          //general purpose buffer
  char filename[1024];
  FILE *lptr; // pointer to log  file

  //----------------------------------------------
  //open log file
  //----------------------------------------------
  printf(         "----------------------------------------------------------------\n");
  sprintf(comment,"Experiment: calculate statistics of length %ld real die sequence",N); printf("%s\n",comment);
  printf(         "----------------------------------------------------------------\n");
  sprintf(filename,"%s_ocs_%ld",basefilename,N);
  lptr=log_open (filename,time1,comment);

  //----------------------------------------------
  //sprintf(buf,"Generate real die sequence...");printofe(lptr,buf,time1);
  //----------------------------------------------
  sprintf(buf,"generate length %ld uniform die sequence generated using seed value 0x%x:\n",N,seed);
  printof(lptr,buf);
  x.randomize(seed);
  x.list(0,299,"","...\n",1,lptr);
  x.histogram(1,lptr);

  //----------------------------------------------
  //plot die sequence
  //----------------------------------------------
  sprintf(comment,"length %ld die sequence",N);
  sprintf(filename,"%s_%x_51_seq.tex",basefilename,seed);
  sprintf(buf,"generate plot file \"%s\"\n",filename); printof(lptr,buf);
  plot_ocs_seq((symseq *)&x, 0,50, time1, "die", filename, comment,lptr);

  //----------------------------------------------
  //plot histogram
  //----------------------------------------------
  sprintf(comment,"HISTOGRAM of length %ld die sequence",N);
  sprintf(filename,"%s_%x_%ld_histo.tex",basefilename,seed,N);
  sprintf(buf,"generate plot file \"%s\"\n",filename); printof(lptr,buf);
  plot_ocs_histo((symseq *)&x, time1, "die", filename, comment,lptr);

  //----------------------------------------------
  sprintf(buf,"Perform auto-correlation operations...\n");printof(lptr,buf);
  //----------------------------------------------
  if(x.Rxxo(&Rxx,1)){fprintf(stderr,"ERROR computing auto-correlation.\n"); return -1;}
  // |      |   |____________switch to turn on counting display
  // |      |________________pointer to output correlation sequence
  // |_______________________input real die sequence
  sprintf(buf,"done");printofe(lptr,buf,time1);
  sprintf(buf,"auto-correlation sequence [r_n|n=N-50...N], where N=%ld:\n",N);
  Rxx.list(N-50,N,buf,"\n",1,lptr);
  sprintf(buf,"auto-correlation sequence [r_n|n=N...N+50], where N=%ld:\n",N);
  Rxx.list(N,N+50,buf,"\n",1,lptr);
  sprintf(comment,"AUTO-CORRELATION of length %ld real die sequence",N);
  sprintf(filename,"%s_%x_%ld_auto.tex",basefilename,seed,N);
  sprintf(buf,"generate plot file \"%s\"\n",filename); printof(lptr,buf);
  plot_ocs_auto(&Rxx, 1, time1, filename,comment,lptr);

  //----------------------------------------------
  // close log file
  //----------------------------------------------
  sprintf(buf,"lab_rdie_ocs experiment complete"); printofe(lptr,buf,time1);
  plot_close(lptr,time1);
  return 0;
  }

//-----------------------------------------------------------------------------
//! \brief generate plot files for spinner sequence
//-----------------------------------------------------------------------------
int lab_spin_ocs(const unsigned seed, const long N, const char *basefilename){
  time_t time1; time(&time1);   //starting time stamp (passed to plotting routine)
  spinseq x(N);   //for die sequence
  seqR1 Rxx(2*N+1); //for auto-correlation of x

  char comment[2*1024];          //comment to be passed to plotting function
  char buf[2*1024];          //general purpose buffer
  char filename[1024];
  FILE *lptr; // pointer to log  file

  //----------------------------------------------
  //open log file
  //----------------------------------------------
  printf(         "---------------------------------------------------------------\n");
  sprintf(comment,"Experiment: calculate statistics of length %ld spinner sequence",N); printf("%s\n",comment);
  printf(         "---------------------------------------------------------------\n");
  sprintf(filename,"%s_ocs_%ld",basefilename,N);
  lptr=log_open (filename,time1,comment);

  //----------------------------------------------
  //sprintf(buf,"Generate spinner sequence...");printofe(lptr,buf,time1);
  //----------------------------------------------
  sprintf(buf,"generate length %ld uniform spinner sequence using seed value 0x%x:\n",N,seed);
  printof(lptr,buf);
  x.randomize(seed);
  x.list(0,299,"","...\n",1,lptr);
  x.histogram(1,lptr);

  //----------------------------------------------
  //plot spinner sequence
  //----------------------------------------------
  sprintf(comment,"length %ld spinner sequence",N);
  sprintf(filename,"%s_%x_51_seq.tex",basefilename,seed);
  sprintf(buf,"generate plot file \"%s\"\n",filename); printof(lptr,buf);
  plot_ocs_seq((symseq *)&x, 0,50, time1, "spin", filename, comment,lptr);

  //----------------------------------------------
  //plot histogram
  //----------------------------------------------
  sprintf(comment,"HISTOGRAM of length %ld spinner sequence",N);
  sprintf(filename,"%s_%x_%ld_histo.tex",basefilename,seed,N);
  sprintf(buf,"generate plot file \"%s\"\n",filename); printof(lptr,buf);
  plot_ocs_histo((symseq *)&x, time1, "spin", filename, comment,lptr);

  //----------------------------------------------
  sprintf(buf,"Perform auto-correlation operations...");printofe(lptr,buf,time1);
  //----------------------------------------------
  if(x.Rxxo(&Rxx,1)){fprintf(stderr,"ERROR computing auto-correlation.\n"); return -1;}
  // |      |   |____________switch to turn on counting display
  // |      |________________pointer to output correlation sequence
  // |_______________________input real die sequence
  sprintf(buf,"done");printofe(lptr,buf,time1);
  sprintf(buf,"auto-correlation sequence [r_n|n=N-50...N], where N=%ld:\n",N);
  Rxx.list(N-50,N,buf,"\n",1,lptr);
  sprintf(buf,"auto-correlation sequence [r_n|n=N...N+50], where N=%ld:\n",N);
  Rxx.list(N,N+50,buf,"\n",1,lptr);
  sprintf(comment,"AUTO-CORRELATION of length %ld spinner sequence",N);
  sprintf(filename,"%s_%x_%ld_auto.tex",basefilename,seed,N);
  sprintf(buf,"generate plot file \"%s\"\n",filename); printof(lptr,buf);
  plot_ocs_auto(&Rxx, 1, time1, filename,comment,lptr);

  //----------------------------------------------
  // close log file
  //----------------------------------------------
  sprintf(buf,"lab_spin_ocs experiment complete"); printofe(lptr,buf,time1);
  plot_close(lptr,time1);
  return 0;
  }

//-----------------------------------------------------------------------------
//! \brief generate plot files for weighted die sequence
//-----------------------------------------------------------------------------
int lab_wdie_ocs(const unsigned seed, const long N, const char *basefilename){
  time_t time1; time(&time1);   //starting time stamp (passed to plotting routine)
  fdieseq x(N);   //for die sequence
  seqR1 Rxx(2*N+1); //for auto-correlation of x

  char comment[2*1024];          //comment to be passed to plotting function
  char buf[2*1024];          //general purpose buffer
  char filename[1024];
  FILE *lptr; // pointer to log  file

  //----------------------------------------------
  //open log file
  //----------------------------------------------
  printf(         "--------------------------------------------------------------------\n");
  sprintf(comment,"Experiment: calculate statistics of length %ld weighted die sequence",N); printf("%s\n",comment);
  printf(         "--------------------------------------------------------------------\n");
  sprintf(filename,"%s_ocs_%ld",basefilename,N);
  lptr=log_open (filename,time1,comment);

  //----------------------------------------------
  sprintf(buf,"Generate weighted die sequence...");printofe(lptr,buf,time1);
  //----------------------------------------------
  x.randomize(seed,5,5,5,5,75,5);
  x.histogram(1,lptr);
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  //----------------------------------------------
  //plot die sequence
  //----------------------------------------------
  sprintf(comment,"length %ld weighted die sequence with weights 5,5,5,5,75,5",N);
  sprintf(filename,"%s_%x_51_seq.tex",basefilename,seed);
  if(plot_ocs_seq((symseq *)&x, 0,50, time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_seq(...)\n"); return -1;}

  //----------------------------------------------
  //plot histogram
  //----------------------------------------------
  sprintf(comment,"HISTOGRAM of length %ld weighted die sequence with weights 5,5,5,5,75,5",N);
  sprintf(filename,"%s_%x_%ld_histo.tex",basefilename,seed,N);
  if(plot_ocs_histo((symseq *)&x, time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_histo(...)\n"); return -1;}

  //----------------------------------------------
  sprintf(buf,"Perform auto-correlation operations...");printofe(lptr,buf,time1);
  //----------------------------------------------
  if(x.Rxxo(&Rxx,1)){fprintf(stderr,"ERROR computing auto-correlation.\n"); return -1;}
  // |       |   |____________switch to turn on counting display
  // |       |________________pointer to output correlation sequence
  // |________________________input fair die sequence
  sprintf(buf,"done.");printofe(lptr,buf,time1);
  sprintf(buf,"auto-correlation sequence [r_n|n=N-50...N], where N=%ld:\n",N);
  Rxx.list(N-50,N,buf,"\n",1,lptr);
  sprintf(buf,"auto-correlation sequence [r_n|n=N...N+50], where N=%ld:\n",N);
  Rxx.list(N,N+50,buf,"\n",1,lptr);
  sprintf(comment,"AUTO-CORRELATION of length %ld weighted die sequence with weights 5,5,5,5,75,5",N);
  sprintf(filename,"%s_%x_%ld_auto.tex",basefilename,seed,N);
  if(plot_ocs_auto(&Rxx, 1, time1, filename,comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_auto(...,%s,...)\n",filename); return -1;}

  //----------------------------------------------
  // close log file
  //----------------------------------------------
  sprintf(buf,"lab_wdie_ocs experiment complete"); printofe(lptr,buf,time1);
  plot_close(lptr,time1);
  return 0;
  }


//-----------------------------------------------------------------------------
//! \brief generate plot files for weighted real die sequence
//-----------------------------------------------------------------------------
int lab_wrdie_ocs(const unsigned seed, const long N, const char *basefilename){
  time_t time1; time(&time1);   //starting time stamp (passed to plotting routine)
  rdieseq x(N);   //for die sequence
  seqR1 Rxx(2*N+1); //for auto-correlation of x

  char comment[2*1024];          //comment to be passed to plotting function
  char buf[2*1024];          //general purpose buffer
  char filename[1024];
  FILE *lptr; // pointer to log  file

  //----------------------------------------------
  //open log file
  //----------------------------------------------
  printf(         "---------------------------------------------------------------------------------------------------------\n");
  sprintf(comment,"Experiment: calculate statistics of length %ld weighted real die sequence with weights*100=(5,5,5,5,75,5)",N); printf("%s\n",comment);
  printf(         "---------------------------------------------------------------------------------------------------------\n");
  sprintf(filename,"%s_ocs_%ld",basefilename,N);
  lptr=log_open (filename,time1,comment);

  //----------------------------------------------
  //sprintf(buf,"Generate weighted die sequence...");printofe(lptr,buf,time1);
  //----------------------------------------------
  sprintf(buf,"length %ld weighted die sequence using seed value 0x%x and weights*100=(5,5,5,5,75,5):\n",N,seed);
  x.randomize(seed,5,5,5,5,75,5);
  x.list(0,299,buf ,"...\n",1,lptr);
  x.histogram(1,lptr);

  //----------------------------------------------
  //plot die sequence
  //----------------------------------------------
  sprintf(comment,"length %ld weighted die sequence with weights 5,5,5,5,75,5",N);
  sprintf(filename,"%s_%x_51_seq.tex",basefilename,seed);
  if(plot_ocs_seq((symseq *)&x, 0,50, time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_seq(...)\n"); return -1;}

  //----------------------------------------------
  //plot histogram
  //----------------------------------------------
  sprintf(comment,"HISTOGRAM of length %ld weighted die sequence with weights 5,5,5,5,75,5",N);
  sprintf(filename,"%s_%x_%ld_histo.tex",basefilename,seed,N);
  if(plot_ocs_histo((symseq *)&x, time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_histo(...)\n"); return -1;}

  //----------------------------------------------
  sprintf(buf,"Perform auto-correlation operations...");printofe(lptr,buf,time1);
  //----------------------------------------------
  if(x.Rxxo(&Rxx,1)){fprintf(stderr,"ERROR computing auto-correlation.\n"); return -1;}
  // |       |   |____________switch to turn on counting display
  // |       |________________pointer to output correlation sequence
  // |________________________input fair die sequence
  sprintf(buf,"done.");printofe(lptr,buf,time1);
  sprintf(buf,"auto-correlation sequence [r_n|n=N-50...N], where N=%ld:\n",N);
  Rxx.list(N-50,N,buf,"\n",1,lptr);
  sprintf(buf,"auto-correlation sequence [r_n|n=N...N+50], where N=%ld:\n",N);
  Rxx.list(N,N+50,buf,"\n",1,lptr);
  sprintf(comment,"AUTO-CORRELATION of length %ld weighted die sequence with weights 5,5,5,5,75,5",N);
  sprintf(filename,"%s_%x_%ld_auto.tex",basefilename,seed,N);
  if(plot_ocs_auto(&Rxx, 1, time1, filename,comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_auto(...,%s,...)\n",filename); return -1;}

  //----------------------------------------------
  // close log file
  //----------------------------------------------
  sprintf(buf,"lab_wrdie_ocs experiment complete"); printofe(lptr,buf,time1);
  plot_close(lptr,time1);
  return 0;
  }

//-----------------------------------------------------------------------------
//! \brief generate plot files for spinner sequence
//-----------------------------------------------------------------------------
int lab_wspin_ocs(const unsigned seed, const long N, const char *basefilename){
  time_t time1; time(&time1);   //starting time stamp (passed to plotting routine)
  spinseq x(N);   //for spinner sequence
  seqR1 Rxx(2*N+1); //for auto-correlation of x

  char comment[2*1024];          //comment to be passed to plotting function
  char buf[2*1024];          //general purpose buffer
  char filename[1024];
  FILE *lptr; // pointer to log  file

  //----------------------------------------------
  //open log file
  //----------------------------------------------
  printf(         "--------------------------------------------------------------------------------------------------------\n");
  sprintf(comment,"Experiment: calculate statistics of length %ld weighted spinner sequence with weights*100=(5,5,5,5,75,5)",N); printf("%s\n",comment);
  printf(         "--------------------------------------------------------------------------------------------------------\n");
  sprintf(filename,"%s_ocs_%ld",basefilename,N);
  lptr=log_open (filename,time1,comment);

  //----------------------------------------------
  //sprintf(buf,"Generate weighted spinner sequence...");printof(lptr,buf);
  //----------------------------------------------
  sprintf(buf,"length %ld weighted spinner sequence using seed value 0x%x and weights*100=(5,5,5,5,75,5):\n",N,seed);
  x.randomize(seed,5,5,5,5,75,5);
  x.list(0,299,buf ,"...\n",1,lptr);
  x.histogram(1,lptr);

  //----------------------------------------------
  //plot spinner sequence
  //----------------------------------------------
  sprintf(comment,"length %ld weighted spinner sequence with weights 5,5,5,5,75,5",N);
  sprintf(filename,"%s_%x_51_seq.tex",basefilename,seed);
  plot_ocs_seq((symseq *)&x, 0,50, time1, "spinner", filename, comment,lptr);

  //----------------------------------------------
  //plot histogram
  //----------------------------------------------
  sprintf(comment,"HISTOGRAM of length %ld weighted spinner sequence with weights 5,5,5,5,75,5",N);
  sprintf(filename,"%s_%x_%ld_histo.tex",basefilename,seed,N);
  plot_ocs_histo((symseq *)&x, time1, "spinner", filename, comment,lptr);

  //----------------------------------------------
  sprintf(buf,"Perform auto-correlation operations...");printofe(lptr,buf,time1);
  //----------------------------------------------
  if(x.Rxxo(&Rxx,1)){fprintf(stderr,"ERROR computing auto-correlation.\n"); return -1;}
  // |       |   |____________switch to turn on counting display
  // |       |________________pointer to output correlation sequence
  // |________________________input fair spinner sequence
  sprintf(buf,"done");printofe(lptr,buf,time1);
  sprintf(buf,"auto-correlation sequence [r_n|n=N-50...N], where N=%ld:\n",N);
  Rxx.list(N-50,N,buf,"\n",1,lptr);
  sprintf(buf,"auto-correlation sequence [r_n|n=N...N+50], where N=%ld:\n",N);
  Rxx.list(N,N+50,buf,"\n",1,lptr);
  sprintf(comment,"AUTO-CORRELATION of length %ld weighted spinner sequence with weights 5,5,5,5,75,5",N);
  sprintf(filename,"%s_%x_%ld_auto.tex",basefilename,seed,N);
  plot_ocs_auto(&Rxx, 1, time1, filename,comment,lptr);

  //----------------------------------------------
  // close log file
  //----------------------------------------------
  sprintf(buf,"lab_wspin_ocs experiment complete"); printofe(lptr,buf,time1);
  plot_close(lptr,time1);
  return 0;
  }

//-----------------------------------------------------------------------------
//! \brief generate plot files for artificial dna sequence
//-----------------------------------------------------------------------------
int lab_dna_ocs(const unsigned seed, const long N, const char *basefilename){
  time_t time1; time(&time1);   //starting time stamp (passed to plotting routine)
  dnaseq x(N);   //for die sequence
  seqR1 Rxx(2*N+1); //for auto-correlation of x

  char comment[2*1024];          //comment to be passed to plotting function
  char buf[2*1024];          //general purpose buffer
  char filename[1024];
  FILE *lptr; // pointer to log  file

  //----------------------------------------------
  //open log file
  //----------------------------------------------
  printf(         "-----------------------------------------------------------\n");
  sprintf(comment,"Experiment: calculate statistics of length %ld DNA sequence",N); printf("%s\n",comment);
  printf(         "-----------------------------------------------------------\n");
  sprintf(filename,"%s_ocs_%ld",basefilename,N);
  lptr=log_open (filename,time1,comment);

  //----------------------------------------------
  //sprintf(buf,"generate length %ld DNA sequence using seed value 0x%x",N,seed);printofe(lptr,buf,time1);
  //----------------------------------------------
  sprintf(buf,"length %ld DNA sequence using seed value 0x%x:\n",N,seed);
  x.randomize(seed);
  x.list(0,299,buf ,"...\n",1,lptr);
  x.histogram(1,lptr);

  //----------------------------------------------
  //plot die sequence
  //----------------------------------------------
  sprintf(comment,"artificial length %ld dna sequence",N);
  sprintf(filename,"%s_%x_51_seq.tex",basefilename,seed);
  if(plot_ocs_seq((symseq *)&x, 0,50, time1, "dna", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_seq(...)\n"); return -1;}

  //----------------------------------------------
  //plot histogram
  //----------------------------------------------
  sprintf(comment,"HISTOGRAM of artificial length %ld dna sequence",N);
  sprintf(filename,"%s_%x_%ld_histo.tex",basefilename,seed,N);
  if(plot_ocs_histo((symseq *)&x, time1, "dna", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_histo(...)\n"); return -1;}

  //----------------------------------------------
  sprintf(buf,"Perform auto-correlation operations...");printofe(lptr,buf,time1);
  //----------------------------------------------
  if(x.Rxxo(&Rxx,1)){fprintf(stderr,"\nERROR computing Rxx.\n"); return -1;} //auto-correlation seqR1 of xR3hlMN
  // |       |   |____________switch to turn on counting display
  // |       |________________pointer to output correlation sequence
  // |________________________input real die sequence
  sprintf(buf,"done.");printofe(lptr,buf,time1);
  sprintf(buf,"auto-correlation sequence [r_n|n=N-50...N], where N=%ld:\n",N);
  Rxx.list(N-50,N,buf,"\n",1,lptr);
  sprintf(buf,"auto-correlation sequence [r_n|n=N...N+50], where N=%ld:\n",N);
  Rxx.list(N,N+50,buf,"\n",1,lptr);
  sprintf(comment,"AUTO-CORRELATION of artificial length %ld dna sequence",N);
  sprintf(filename,"%s_%x_%ld_auto.tex",basefilename,seed,N);
  if(plot_ocs_auto(&Rxx, 1, time1, filename,comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_auto(...,%s,...)\n",filename); return -1;}

  //----------------------------------------------
  // close log file
  //----------------------------------------------
  sprintf(buf,"lab_dna_ocs experiment complete"); printofe(lptr,buf,time1);
  plot_close(lptr,time1);
  return 0;
  }

//-----------------------------------------------------------------------------
//! \brief generate plot files for dna sequence read from FASTA file
//-----------------------------------------------------------------------------
int lab_dna_ocs(const char *datafilename, const char *basefilename){
  time_t time1; time(&time1);   //starting time stamp (passed to plotting routine)
  const long N=numsym_fasta_file(datafilename);
  dnaseq x(N);   //for die sequence
  seqR1 Rxx(2*N+1); //for auto-correlation of x
  char header[1024];
  if(N==-1)return -1;
  FILE *lptr; // pointer to log  file
  char comment[2*1024];          //comment to be passed to plotting function
  char buf[2*1024];          //general purpose buffer
  char filename[1024];

  //----------------------------------------------
  //open log file
  //----------------------------------------------
  printf(         "-----------------------------------------------------------\n");
  sprintf(comment,"calculate statistics of length %ld dna sequence from FASTA file %s\n%s",N,datafilename,header); printf("%s\n",comment);
  printf(         "-----------------------------------------------------------\n");
  sprintf(buf,"%s",basefilename);
  lptr=log_open (buf,time1,comment);

  //----------------------------------------------
  // read data from FASTA file
  //----------------------------------------------
  sprintf(buf,"length %ld DNA sequence from FASTA file \"%s\":\n",N,datafilename);
  read_fasta_file(datafilename,header,&x);
  x.list(0,299,buf ,"...\n",1,lptr);
  x.histogram(1,lptr);

  //----------------------------------------------
  //plot dna sequence
  //----------------------------------------------
  sprintf(comment,"length %ld dna sequence from %s. %s",N,datafilename,header);
  sprintf(filename,"%s_51_seq.tex",basefilename);
  if(plot_ocs_seq((symseq *)&x, 0,50, time1, "dna", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_seq(...)\n"); return -1;}

  //----------------------------------------------
  //plot histogram
  //----------------------------------------------
  sprintf(comment,"HISTOGRAM of length %ld dna sequence from %s. %s",N,datafilename,header);
  sprintf(filename,"%s_histo.tex",basefilename);
  if(plot_ocs_histo((symseq *)&x, time1, "dna", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_histo(...)\n"); return -1;}

  //----------------------------------------------
  sprintf(buf,"Perform auto-correlation operations...");printofe(lptr,buf,time1);
  //----------------------------------------------
  if(x.Rxxo(&Rxx,1)){fprintf(stderr,"\nERROR computing Rxx.\n"); return -1;} //auto-correlation seqR1 of xR3hlMN
  // |       |   |____________switch to turn on counting display
  // |       |________________pointer to output correlation sequence
  // |________________________input real die sequence
  sprintf(buf,"done.");printofe(lptr,buf,time1);
  sprintf(buf,"auto-correlation sequence [r_n|n=N-50...N], where N=%ld:\n",N);
  Rxx.list(N-50,N,buf,"\n",1,lptr);
  sprintf(buf,"auto-correlation sequence [r_n|n=N...N+50], where N=%ld:\n",N);
  Rxx.list(N,N+50,buf,"\n",1,lptr);
  sprintf(comment,"AUTO-CORRELATION of length %ld dna sequence from %s. %s",N,datafilename,header);
  sprintf(filename,"%s_auto.tex",basefilename);
  if(plot_ocs_auto(&Rxx, 1, time1, filename,comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_auto(...,%s,...)\n",filename); return -1;}

  //----------------------------------------------
  // close log file
  //----------------------------------------------
  sprintf(buf,"lab_dna_ocs experiment complete"); printofe(lptr,buf,time1);
  plot_close(lptr,time1);
  return 0;
  }

//-----------------------------------------------------------------------------
//! \brief generate plot files for scaffold dna sequence read from FASTA file
//-----------------------------------------------------------------------------
int lab_dnan_ocs(const char *datafilename, const char *basefilename){
  time_t time1; time(&time1);   //starting time stamp (passed to plotting routine)
  const long N=numsym_fasta_file(datafilename);
  dnanseq x(N);   //for die sequence
  seqR1 Rxx(2*N+1); //for auto-correlation of x
  char header[1024];
  if(N==-1)return -1;
  FILE *lptr; // pointer to log  file
  char comment[2*1024];          //comment to be passed to plotting function
  char buf[2*1024];          //general purpose buffer
  char filename[1024];

  //----------------------------------------------
  //open log file
  //----------------------------------------------
  printf(         "----------------------------------------------------------------------------\n");
  sprintf(comment,"Experiment: calculate statistics of length %ld dnan sequence from FASTA file %s\n%s\n",N,datafilename,header);printf("%s",comment);
  printf(         "----------------------------------------------------------------------------\n");
  sprintf(buf,"%s_ocs",basefilename);
  lptr=log_open (buf,time1,comment);

  //----------------------------------------------
  // read data from FASTA file
  //----------------------------------------------
  sprintf(buf,"length %ld DNA sequence from FASTA file \"%s\":\n",N,datafilename);
  read_fasta_file(datafilename,header,&x);
  x.list(0,299,buf ,"...\n",1,lptr);
  x.histogram(1,lptr);

  //----------------------------------------------
  //plot die sequence
  //----------------------------------------------
  sprintf(comment,"length %ld scaffold dna sequence from %s. %s",N,datafilename,header);
  //printf("%s",comment);
  sprintf(filename,"%s_51_seq.tex",basefilename);
  if(plot_ocs_seq((symseq *)&x, 0,50, time1, "dnan", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_seq(...)\n"); return -1;}

  //----------------------------------------------
  //plot histogram
  //----------------------------------------------
  sprintf(comment,"HISTOGRAM of length %ld scaffold dna sequence from %s. %s",N,datafilename,header);
  printf("%s",comment);
  sprintf(filename,"%s_histo.tex",basefilename);
  if(plot_ocs_histo((symseq *)&x, time1, "dnan", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_histo(...)\n"); return -1;}

  //----------------------------------------------
  sprintf(buf,"Perform auto-correlation operations...");printofe(lptr,buf,time1);
  //----------------------------------------------
  if(x.Rxxo(&Rxx,1)){fprintf(stderr,"\nERROR computing Rxx.\n"); return -1;} //auto-correlation seqR1 of xR3hlMN
  // |       |   |____________switch to turn on counting display
  // |       |________________pointer to output correlation sequence
  // |________________________input real die sequence
  sprintf(buf,"done.");printofe(lptr,buf,time1);
  sprintf(buf,"auto-correlation sequence [r_n|n=N-50...N], where N=%ld:\n",N);
  Rxx.list(N-50,N,buf,"\n",1,lptr);
  sprintf(buf,"auto-correlation sequence [r_n|n=N...N+50], where N=%ld:\n",N);
  Rxx.list(N,N+50,buf,"\n",1,lptr);
  sprintf(comment,"AUTO-CORRELATION of length %ld scaffold dna sequence from %s. %s",N,datafilename,header);
  sprintf(filename,"%s_auto.tex",basefilename);
  if(plot_ocs_auto(&Rxx, 1, time1, filename,comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_auto(...,%s,...)\n",filename); return -1;}

  //----------------------------------------------
  // close log file
  //----------------------------------------------
  sprintf(buf,"lab_dnan_ocs experiment complete"); printofe(lptr,buf,time1);
  plot_close(lptr,time1);
  return 0;
  }

//-----------------------------------------------------------------------------
//! \brief generate TeX plot file for the low pass filtering of real die sequence
//!        filtered by M tap rectangular filter,
//!        and generate TeX plot file
//-----------------------------------------------------------------------------
int lab_rdie_lp(const unsigned seed, const long N, const long M, const char *basefilename){
  time_t time1; time(&time1);   //starting time stamp (passed to plotting routine)

  rdieseq x(N+M);        //for die sequence
  seqR1   xR1(N+M);      //for R^1 seqR1 mapped from x sequence
  seqR3   xR3(N+M);      //for R^3 seqR1 mapped from x sequence
  seqR1   r(M);          //for M-tap rectangular filter 
  seqR1   h(M);          //for M-tap Hanning filter
  seqR1   xR1r(N+M+M-1); //for x mapped to R1 and filtered by r
  seqR1   xR1h(N+M+M-1); //for x mapped to R1 and filtered by h
  seqR3   xR3r(N+M+M-1); //for x mapped to R3 and filtered by r
  seqR3   xR3h(N+M+M-1); //for x mapped to R3 and filtered by h
               //     |  |_______________M-1: convolution with h appends M-1 elements
               //     |__________________ Nm: length of xR1 sequence
  rdieseq xR1re(N+M+M-1);   //for xR1r mapped back using Euclidean metric
  rdieseq xR1he(N+M+M-1);   //for xR1h mapped back using Euclidean metric
  rdieseq xR3re(N+M+M-1);   //for xR3r mapped back using Euclidean metric
  rdieseq xR3rl(N+M+M-1);   //for xR3r mapped back using Lagrange arc distance
  rdieseq xR3he(N+M+M-1);   //for xR3h mapped back using Euclidean metric
  rdieseq xR3hl(N+M+M-1);   //for xR3h mapped back using Lagrange arc distance
                           
  rdieseq xR1reN(N);       //for xR1re with beginning M-1 and ending M-1 elements removed
  rdieseq xR1heN(N);       //for xR1he with beginning M-1 and ending M-1 elements removed
  rdieseq xR3reN(N);       //for xR3re with beginning M-1 and ending M-1 elements removed
  rdieseq xR3rlN(N);       //for xR3rl with beginning M-1 and ending M-1 elements removed
  rdieseq xR3heN(N);       //for xR3he with beginning M-1 and ending M-1 elements removed
  rdieseq xR3hlN(N);       //for xR3hl with beginning M-1 and ending M-1 elements removed
                           
  seqR1 RxR1re(2*N+1);  //for auto-correlation seqR1 of xR1re 
  seqR1 RxR1he(2*N+1);  //for auto-correlation seqR1 of xR1he 
  seqR1 RxR3re(2*N+1);  //for auto-correlation seqR1 of xR3re 
  seqR1 RxR3rl(2*N+1);  //for auto-correlation seqR1 of xR3rl 
  seqR1 RxR3he(2*N+1);  //for auto-correlation seqR1 of xR3he 
  seqR1 RxR3hl(2*N+1);  //for auto-correlation seqR1 of xR3hl 

  char comment[2*1024];          //comment to be passed to plotting function
  char buf[2*1024];          //comment to be passed to plotting function
  long n;
  char filename[1024];
  FILE *lptr; // pointer to log  file

  //----------------------------------------------
  //open log file
  //----------------------------------------------
  printf(         "---------------------------------------------------------------\n");
  sprintf(comment,"Experiment: low pass filter length %6ld real die sequence",N);printf("%s\n",comment);
  printf(         "---------------------------------------------------------------\n");
  sprintf(buf,"%s_%ldm%ld",basefilename,N,M);
  lptr=log_open (buf,time1,comment);

  //----------------------------------------------
  // generate die sequence
  //----------------------------------------------
  sprintf(buf,"length %ld uniform die sequence generated using seed value 0x%x:\n",N,seed);
  printof(lptr,buf);
  x.randomize(seed);
  x.list(0,299,buf,"...\n",1,lptr);
  x.histogram(1,lptr);

  //----------------------------------------------
  // map die sequence to R^n sequences
  //----------------------------------------------
  xR1=x.dietoR1();               //xR1 = x seqR1 mapped to R^1 sequence
  xR3=x.dietoR3();               //xR3 = x die seqR1 mapped to R^3 sequence
  xR1.list(0,24,"die sequence mapped to R^1 sequence:\n","...\n",1,lptr);
  xR3.list(0,24,"die sequence mapped to R^3 sequence:\n","...\n",1,lptr);

  //----------------------------------------------
  // generate filter coefficients
  //----------------------------------------------
  r=1.0;       
  h.hanning(); 
  sprintf(buf,"length %ld rectangular sequence before normalization:",M); r.list(buf,"\n",1,lptr);
  sprintf(buf,"length %ld Hanning sequence before normalization:",M);    h.list(buf,"\n",1,lptr);
  r.normalize(); 
  h.normalize(); 
  sprintf(buf,"length %ld rectangular sequence after normalization:",M); r.list(buf,"\n",1,lptr);
  sprintf(buf,"length %ld Hanning sequence after normalization:",M);    h.list(buf,"\n",1,lptr);
  //r=1.0;       r.normalize();
  //h.hanning(); h.normalize();
  //r.list("rectangular filter sequence:\n","\n",1,lptr);
  //h.list("Hanning filter sequence:\n","\n",1,lptr);

  //----------------------------------------------
  // filter mapped sequences
  //----------------------------------------------
  printf("xR1*r..."); convolve(&xR1, &r, &xR1r); //xR1r=xR1*r; 
  printf("xR1*h..."); convolve(&xR1, &h, &xR1h); //xR1h=xR1*h; 
  printf("xR3*r..."); convolve(&xR3, &r, &xR3r); //xR3r=xR3*r; 
  printf("xR3*h..."); convolve(&xR3, &h, &xR3h); //xR3h=xR3*h; 
  xR1r.list(0,9,"\nR^1 seq. mapped from die seq. and filtered by rectangular filter:\n","...\n",1,lptr);
  xR1h.list(0,9,"\nR^1 seq. mapped from die seq. and filtered by Hanning filter:\n","...\n",1,lptr);
  xR3r.list(0,58,"\nR^3 seq. mapped from die seq. and filtered by rectangular filter:\n","...\n",1,lptr);
  xR3h.list(0,8,"\nR^3 seq. mapped from die seq. and filtered by Hanning filter:\n","...\n",1,lptr);

  //----------------------------------------------
  // mapping operations from R^n back to die sequence space
  //----------------------------------------------
  xR1re=rdie_R1todie_euclid(xR1r); xR1re.list(0,249,"die seq. filtered in R^1 using rect. seq. and mapped back using Euclidean\n","...\n",1,lptr);
  xR1he=rdie_R1todie_euclid(xR1h); xR1he.list(0,249,"die seq. filtered in R^1 using Hann. seq. and mapped back using Euclidean\n","...\n",1,lptr);
  xR3re=rdie_R3todie_euclid(xR3r); xR3re.list(0,249,"die seq. filtered in R^3 using rect. seq. and mapped back using Euclidean\n","...\n",1,lptr);
  xR3rl=rdie_R3todie_larc  (xR3r); xR3rl.list(0,249,"die seq. filtered in R^3 using rect. seq. and mapped back using Lagrange\n","...\n",1,lptr);
  xR3he=rdie_R3todie_euclid(xR3h); xR3he.list(0,249,"die seq. filtered in R^3 using Hann. seq. and mapped back using Euclidean\n","...\n",1,lptr);
  xR3hl=rdie_R3todie_larc  (xR3h); xR3hl.list(0,249,"die seq. filtered in R^3 using Hann. seq. and mapped back using Lagrange\n","...\n",1,lptr);

  //----------------------------------------------
  //Compare Euclidean and Lagrange mappings
  //----------------------------------------------
  n=cmp(&xR3re,&xR3rl,0,lptr);
  if(n==0l){sprintf(buf,"Euclidean and Lagrange sequences after using %ld-tap rectangular filter are identical.\n",M);printofe(lptr,buf,time1);}
  else     {sprintf(buf,"Euclidean and Lagrange sequences after using %ld-tap rectangular filter differ at %ld locations.\n",M,n);printofe(lptr,buf,time1);}
  n=cmp(&xR3he,&xR3hl,0,lptr);
  if(n==0l){sprintf(buf,"Euclidean and Lagrange sequences after using %ld-tap Hanning filter are identical.\n",M);printofe(lptr,buf,time1);}
  else     {sprintf(buf,"Euclidean and Lagrange sequences after using %ld-tap Hanning filter differ at %ld locations.\n",M,n);printofe(lptr,buf,time1);}

  //----------------------------------------------
  sprintf(buf,"Remove beginning %ld element and ending %ld element transition regions introduced by filtering ... ",M-1,M);printof(lptr,buf);
  //----------------------------------------------
  copy(M-1, M-1+N-1, &xR1re, &xR1reN);
  copy(M-1, M-1+N-1, &xR1he, &xR1heN);
  copy(M-1, M-1+N-1, &xR3re, &xR3reN);
  copy(M-1, M-1+N-1, &xR3rl, &xR3rlN);
  copy(M-1, M-1+N-1, &xR3he, &xR3heN);
  copy(M-1, M-1+N-1, &xR3hl, &xR3hlN);

  //for(n=0;n<N;n++){
  //  xR1reN.put(n,xR1re.get(n+M-1));
  //  xR1heN.put(n,xR1he.get(n+M-1));
  //  xR3reN.put(n,xR3re.get(n+M-1));
  //  xR3rlN.put(n,xR3rl.get(n+M-1));
  //  xR3heN.put(n,xR3he.get(n+M-1));
  //  xR3hlN.put(n,xR3hl.get(n+M-1));
  //  }
  printof(lptr,"done.\n");

  //----------------------------------------------
  printf("Perform auto-correlation operations...\n");
  //----------------------------------------------
  if(xR1reN.Rxxo(&RxR1re,1)){fprintf(stderr,"ERROR computing Rxx for x=xR1re.\n"); return -1;} //auto-correlation seqR1 of truncated xR1re
  if(xR1heN.Rxxo(&RxR1he,1)){fprintf(stderr,"ERROR computing Rxx for x=xR1he.\n"); return -1;} //auto-correlation seqR1 of truncated xR1he
  if(xR3reN.Rxxo(&RxR3re,1)){fprintf(stderr,"ERROR computing Rxx for x=xR3re.\n"); return -1;} //auto-correlation seqR1 of truncated xR3re
  if(xR3rlN.Rxxo(&RxR3rl,1)){fprintf(stderr,"ERROR computing Rxx for x=xR3rl.\n"); return -1;} //auto-correlation seqR1 of truncated xR3rl
  if(xR3heN.Rxxo(&RxR3he,1)){fprintf(stderr,"ERROR computing Rxx for x=xR3he.\n"); return -1;} //auto-correlation seqR1 of truncated xR3he
  if(xR3hlN.Rxxo(&RxR3hl,1)){fprintf(stderr,"ERROR computing Rxx for x=xR3hl.\n"); return -1;} //auto-correlation seqR1 of truncated xR3hl
  // |            |      |____________switch to turn on counting display
  // |            |___________________pointer to output correlation sequence
  // |________________________________input real die sequence
  printf("done.\n");

  sprintf(buf,"R^1, rect., Euclidean auto-correlation sequence about N=%ld: ...\n",N);
  RxR1re.list(N-50,N,buf,"\n",1,lptr);
  RxR1re.list(N,N+50,"","...\n",1,lptr);

  sprintf(buf,"R^1, Hann., Euclidean auto-correlation sequence about N=%ld: ...\n",N);
  RxR1he.list(N-50,N,buf,"\n",1,lptr);
  RxR1he.list(N,N+50,"","...\n",1,lptr);

  sprintf(buf,"R^3, rect., Euclidean auto-correlation sequence about N=%ld: ...\n",N);
  RxR3re.list(N-50,N,buf,"\n",1,lptr);
  RxR3re.list(N,N+50,"","...\n",1,lptr);

  sprintf(buf,"R^3, rect., Lagrange auto-correlation sequence about N=%ld: ...\n",N);
  RxR3rl.list(N-50,N,buf,"\n",1,lptr);
  RxR3rl.list(N,N+50,"","...\n",1,lptr);

  sprintf(buf,"R^3, Hann., Euclidean auto-correlation sequence about N=%ld: ...\n",N);
  RxR3he.list(N-50,N,buf,"\n",1,lptr);
  RxR3he.list(N,N+50,"","...\n",1,lptr);

  sprintf(buf,"R^3, Hann., Lagrange auto-correlation sequence about N=%ld:...\n",N);
  RxR3hl.list(N-50,N,buf,"\n",1,lptr);
  RxR3hl.list(N,N+50,"","...\n",1,lptr);

  //----------------------------------------------
  //plot seqR1 tex files
  //----------------------------------------------
  sprintf(comment,"length %ld real die SEQUENCE filtered in R^1 using %ld tap Rectangular filter and mapped backc to real die seqR1 using Euclidean metric",N,M);
  sprintf(filename,"%s_%ld_R1_rect%ld_euclid_seq.tex",basefilename,N,M);
  printf("Plot file %50s ... ",filename);
  if(plot_ocs_seq((symseq *)&xR1re, M-1, M-1+48,time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_seq(...)\n"); return -1;}
  printf("done.\n");

  sprintf(comment,"length %ld real die SEQUENCE filtered in R^1 using %ld tap Hanning filter and mapped backc to real die seqR1 using Euclidean metric",N,M);
  sprintf(filename,"%s_%ld_R1_hann%ld_euclid_seq.tex",basefilename,N,M);
  printf("Plot file %50s ... ",filename);
  if(plot_ocs_seq((symseq *)&xR1he, M-1, M-1+48,time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_seq(...)\n"); return -1;}
  printf("done.\n");

  sprintf(comment,"length %ld real die SEQUENCE filtered in R^3 using %ld tap rectangular filter and mapped back to real die seqR1 using Euclidean metric",N,M);
  sprintf(filename,"%s_%ld_R3_rect%ld_euclid_seq.tex",basefilename,N,M);
  printf("Plot file %50s ... ",filename);
  if(plot_ocs_seq((symseq *)&xR3re, M-1, M-1+48,time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_seq(...)\n"); return -1;}
  printf("done.\n");

  sprintf(comment,"length %ld real die SEQUENCE filtered in R^3 using %ld tap rectangular filter and mapped back to real die seqR1 using Lagrange arc distance",N,M);
  sprintf(filename,"%s_%ld_R3_rect%ld_larc_seq.tex",basefilename,N,M);
  printf("Plot file %50s ... ",filename);
  if(plot_ocs_seq((symseq *)&xR3rl, M-1, M-1+48,time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_seq(...)\n"); return -1;}
  printf("done.\n");

  sprintf(comment,"length %ld real die SEQUENCE filtered in R^3 using %ld tap Hanning filter and mapped back to real die seqR1 using Euclidean metric",N,M);
  sprintf(filename,"%s_%ld_R3_hann%ld_euclid_seq.tex",basefilename,N,M);
  printf("Plot file %50s ... ",filename);
  if(plot_ocs_seq((symseq *)&xR3he, M-1, M-1+48,time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_seq(...)\n"); return -1;}
  printf("done.\n");

  sprintf(comment,"length %ld real die SEQUENCE filtered in R^3 using %ld tap Hanning filter and mapped back to real die seqR1 using Lagrange arc distance",N,M);
  sprintf(filename,"%s_%ld_R3_hann%ld_larc_seq.tex",basefilename,N,M);
  printf("Plot file %50s ... ",filename);
  if(plot_ocs_seq((symseq *)&xR3hl, M-1, M-1+48,time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_seq(...)\n"); return -1;}
  printf("done.\n");

  //if(xR3re==xR3rl)fprintf(lptr,"Euclidean and Lagrange sequences are identical.\n");
  //else            fprintf(lptr,"Euclidean and Lagrange sequences are different.\n");
  //if(xR3re==xR3rl)printf("Euclidean and Lagrange sequences are identical.\n");
  //else            printf("Euclidean and Lagrange sequences are different.\n");

  //----------------------------------------------
  //plot histogram tex files
  //----------------------------------------------
  sprintf(comment,"HISTOGRAM of length %ld real die seqR1 filtered in R^1 using %ld tap Rectangular filter",N,M);
  sprintf(filename,"%s_%ld_R1_rect%ld_euclid_histo.tex",basefilename,N,M);
  printf("Plot file %50s ... ",filename);
  if(plot_ocs_histo((symseq *)&xR1re, M-1, N+M-2, time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_histo(...)\n"); return -1;}
  printf("done.\n");

  sprintf(comment,"HISTOGRAM of length %ld real die seqR1 filtered in R^1 using %ld tap Hanning filter",N,M);
  sprintf(filename,"%s_%ld_R1_hann%ld_euclid_histo.tex",basefilename,N,M);
  printf("Plot file %50s ... ",filename);
  if(plot_ocs_histo((symseq *)&xR1he, M-1, N+M-2, time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_histo(...)\n"); return -1;}
  printf("done.\n");

  sprintf(comment,"HISTOGRAM of length %ld real die seqR1 filtered in R^3 using %ld tap rectangular filter and mapped back to real die seqR1 using Euclidean metric",N,M);
  sprintf(filename,"%s_%ld_R3_rect%ld_euclid_histo.tex",basefilename,N,M);
  printf("Plot file %50s ... ",filename);
  if(plot_ocs_histo((symseq *)&xR3re, M-1, N+M-2, time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_histo(...)\n"); return -1;}
  printf("done.\n");

  sprintf(comment,"HISTOGRAM of length %ld real die seqR1 filtered in R^3 using %ld tap rectangular filter and mapped back to real die seqR1 using Lagrange arc distance",N,M);
  sprintf(filename,"%s_%ld_R3_rect%ld_larc_histo.tex",basefilename,N,M);
  printf("Plot file %50s ... ",filename);
  if(plot_ocs_histo((symseq *)&xR3rl, M-1, N+M-2, time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_histo(...)\n"); return -1;}
  printf("done.\n");

  sprintf(comment,"HISTOGRAM of length %ld real die seqR1 filtered in R^3 using %ld tap Hanning filter and mapped back to real die seqR1 using Euclidean metric",N,M);
  sprintf(filename,"%s_%ld_R3_hann%ld_euclid_histo.tex",basefilename,N,M);
  printf("Plot file %50s ... ",filename);
  if(plot_ocs_histo((symseq *)&xR3he, M-1, N+M-2, time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_histo(...)\n"); return -1;}
  printf("done.\n");

  sprintf(comment,"HISTOGRAM of length %ld real die seqR1 filtered in R^3 using %ld tap Hanning filter and mapped back to real die seqR1 using Lagrange arc distance",N,M);
  sprintf(filename,"%s_%ld_R3_hann%ld_larc_histo.tex",basefilename,N,M);
  printf("Plot file %50s ... ",filename);
  if(plot_ocs_histo((symseq *)&xR3hl, M-1, N+M-2, time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_histo(...)\n"); return -1;}
  printf("done.\n");

  //----------------------------------------------
  //plot autocorrelation tex files
  //----------------------------------------------
  sprintf(comment,"AUTO-CORRELATION of length %ld real die seqR1 filtered in R^1 using %ld tap Rectangular filter",N,M);
  sprintf(filename,"%s_%ld_R1_rect%ld_euclid_auto.tex",basefilename,N,M);
  printf("Plot file %50s ... ",filename);
  if(plot_ocs_auto(&RxR1re, M, time1, filename,comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_auto(...,%s,...)\n",filename); return -1;}
  printf("done.\n");

  sprintf(comment,"AUTO-CORRELATION of length %ld real die seqR1 filtered in R^1 using %ld tap Hanning filter",N,M);
  sprintf(filename,"%s_%ld_R1_hann%ld_euclid_auto.tex",basefilename,N,M);
  printf("Plot file %50s ... ",filename);
  if(plot_ocs_auto(&RxR1he, M, time1, filename,comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_auto(...,%s,...)\n",filename); return -1;}
  printf("done.\n");

  sprintf(comment,"AUTO-CORRELATION of length %ld real die seqR1 filtered in R^3 using %ld tap rectangular filter and mapped back to real die seqR1 using Euclidean metric",N,M);
  sprintf(filename,"%s_%ld_R3_rect%ld_euclid_auto.tex",basefilename,N,M);
  printf("Plot file %50s ... ",filename);
  if(plot_ocs_auto(&RxR3re, M, time1, filename,comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_auto(...,%s,...)\n",filename); return -1;}
  printf("done.\n");

  sprintf(comment,"AUTO-CORRELATION of length %ld real die seqR1 filtered in R^e using %ld tap rectangular filter and mapped back to real die seqR1 using Lagrange arc distance",N,M);
  sprintf(filename,"%s_%ld_R3_rect%ld_larc_auto.tex",basefilename,N,M);
  printf("Plot file %50s ... ",filename);
  if(plot_ocs_auto(&RxR3rl, M, time1, filename,comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_auto(...,%s,...)\n",filename); return -1;}
  printf("done.\n");

  sprintf(comment,"AUTO-CORRELATION of length %ld real die seqR1 filtered in R^3 using %ld tap Hanning filter and mapped back to real die seqR1 using Euclidean metric",N,M);
  sprintf(filename,"%s_%ld_R3_hann%ld_euclid_auto.tex",basefilename,N,M);
  printf("Plot file %50s ... ",filename);
  if(plot_ocs_auto(&RxR3he, M, time1, filename,comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_auto(...,%s,...)\n",filename); return -1;}
  printf("done.\n");

  sprintf(comment,"AUTO-CORRELATION of length %ld real die seqR1 filtered in R^3 using %ld tap Hanning filter and mapped back to real die seqR1 using Lagrange arc distance",N,M);
  sprintf(filename,"%s_%ld_R3_hann%ld_larc_auto.tex",basefilename,N,M);
  printf("Plot file %50s ... ",filename);
  if(plot_ocs_auto(&RxR3hl, M, time1, filename,comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_auto(...,%s,...)\n",filename); return -1;}
  printf("done.\n");

  //----------------------------------------------
  // close log file
  //----------------------------------------------
  sprintf(buf,"lab_rdie_lp experiment complete"); printofe(lptr,buf,time1);
  plot_close(lptr,time1);
  return 0;
  }

//-----------------------------------------------------------------------------
//! \brief generate TeX plot file for the high pass filtering of weighted real die sequence
//!        filtered by M tap rectangular filter,
//!        and generate TeX plot file
//-----------------------------------------------------------------------------
int lab_wrdie_hp(const unsigned seed, const long N, const long M, const char *basefilename){
  time_t time1; time(&time1);   //starting time stamp (passed to plotting routine)
  rdieseq x(M*(N+2)-(M-1));//for real die sequence
  seqR1 r(M);              //for M-tap rectangular filter
  seqR1 h(M);              //for M-tap Hanning filter

  seqR1   xR1(M*(N+2)-(M-1));//for x mapped to R^1
  seqR3   xR3(M*(N+2)-(M-1));//for x mapped to R^3

  seqR1   xR1r(M*(N+2)); //for xR1 filtered by r
  seqR1   xR1h(M*(N+2)); //for xR1 filtered by h
  seqR3   xR3r(M*(N+2)); //for xR3 filtered by r
  seqR3   xR3h(M*(N+2)); //for xR3 filtered by h

  rdieseq xR1re(M*(N+2));//for xR1r mapped back to die seqR1 using Euclidean metric
  rdieseq xR1he(M*(N+2));//for xR1h mapped back to die seqR1 using Euclidean metric
  rdieseq xR3re(M*(N+2));//for xR3r mapped back to die seqR1 using Euclidean metric
  rdieseq xR3he(M*(N+2));//for xR3h mapped back to die seqR1 using Euclidean metric
  rdieseq xR3rl(M*(N+2));//for xR3r mapped back to die seqR1 using Lagrange  distance
  rdieseq xR3hl(M*(N+2));//for xR3h mapped back to die seqR1 using Lagrange  distance
          
  rdieseq xR1reM(N+2);   //for xR1re downsampled by factor of M
  rdieseq xR1heM(N+2);   //for xR1he downsampled by factor of M
  rdieseq xR3reM(N+2);   //for xR3re downsampled by factor of M
  rdieseq xR3heM(N+2);   //for xR3he downsampled by factor of M
  rdieseq xR3rlM(N+2);   //for xR3rl downsampled by factor of M
  rdieseq xR3hlM(N+2);   //for xR3hl downsampled by factor of M
          
  rdieseq xR1reMN(N);    //for xR1reM with first and last element removed
  rdieseq xR1heMN(N);    //for xR1heM with first and last element removed
  rdieseq xR3reMN(N);    //for xR3reM with first and last element removed
  rdieseq xR3heMN(N);    //for xR3heM with first and last element removed
  rdieseq xR3rlMN(N);    //for xR3rlM with first and last element removed
  rdieseq xR3hlMN(N);    //for xR3hlM with first and last element removed

  seqR1   RxR1re(2*N+1); //for auto-correlation of xR1reMN
  seqR1   RxR1he(2*N+1); //for auto-correlation of xR1heMN
  seqR1   RxR3re(2*N+1); //for auto-correlation of xR3reMN
  seqR1   RxR3he(2*N+1); //for auto-correlation of xR3heMN
  seqR1   RxR3rl(2*N+1); //for auto-correlation of xR3rlMN
  seqR1   RxR3hl(2*N+1); //for auto-correlation of xR3hlMN

  char comment[2*1024];          //comment to be passed to plotting function
  char buf[2*1024];          //general purpose buffer
  long n;
  char filename[1024];
  FILE *lptr; // pointer to log  file

  //----------------------------------------------
  //open log file
  //----------------------------------------------
  printf(         "-----------------------------------------------------------\n");
  sprintf(comment,"Experiment: map nonstationary length %ld die sequence\n%%  with length %ld periodic distribution into R^1 using PAM mapping\n%%  and use DFT:R^1-->C^1 to analyze",N,M); printf("%s",comment);
  printf(         "-----------------------------------------------------------\n");
  sprintf(buf,"%s_%ldm%ld",basefilename,N,M);
  lptr=log_open (buf,time1,comment);

  //----------------------------------------------
  sprintf(buf,"Generate weighted die sequence...");printofe(lptr,buf,time1);
  //----------------------------------------------
  x.randomize(seed,5,5,5,5,75,5);
  x.histogram(1,lptr);
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  //----------------------------------------------
  sprintf(buf,"Map from die space ");printofe(lptr,buf,time1);
  //----------------------------------------------
  printf("to R^1..."); xR1=x.dietoR1();  //xR1 = x seqR1 mapped to R^1 sequence
  printf("to R^3..."); xR3=x.dietoR3();  //xR3 = x die seqR1 mapped to R^3 sequence
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  //----------------------------------------------
  sprintf(buf,"Generate filter coefficients...");printofe(lptr,buf,time1);
  //----------------------------------------------
  r=1.0;       r.lptohp(); printf("\n  r_n=");r.list();
  h.hanning(); h.lptohp(); printf("\n  h_n=");h.list();putchar('\n'); 
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  //----------------------------------------------
  sprintf(buf,"Perform filtering operations...");printofe(lptr,buf,time1);
  //----------------------------------------------
  printf("xR1*r..."); convolve(&xR1, &r, &xR1r); //xR1r=xR1*r; 
  printf("xR1*h..."); convolve(&xR1, &h, &xR1h); //xR1h=xR1*h; 
  printf("xR3*r..."); convolve(&xR3, &r, &xR3r); //xR3r=xR3*r; 
  printf("xR3*h..."); convolve(&xR3, &h, &xR3h); //xR3h=xR3*h; 
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  //----------------------------------------------
  sprintf(buf,"Perform mapping operations from R^n back to die seqR1 space...");printofe(lptr,buf,time1);
  //----------------------------------------------
  xR1re=rdie_R1todie_euclid(xR1r);//xR1r mapped back to real die seq. using Euclid
  xR1he=rdie_R1todie_euclid(xR1h);//xR1h mapped back to real die seq. using Euclid
  xR3re=rdie_R3todie_euclid(xR3r);//xR3r mapped back to real die seq. using Euclid
  xR3rl=rdie_R3todie_larc  (xR3r);//xR3r mapped back to real die seq. using Lagrange
  xR3he=rdie_R3todie_euclid(xR3h);//xR3h mapped back to real die seq. using Euclid
  xR3hl=rdie_R3todie_larc  (xR3h);//xR3h mapped back to real die seq. using Lagrange
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  //----------------------------------------------
  sprintf(buf,"Compare Euclidean and Lagrange mappings after %ld-tap rectangular filtering...\n  ",M);printofe(lptr,buf,time1);
  //----------------------------------------------
  n=cmp(&xR3re,&xR3rl,0,lptr);
  if(n==0l){sprintf(buf,"Euclidean and Lagrange sequences after using %ld-tap rectangular filter are identical.\n",M);printofe(lptr,buf,time1);}
  else     {sprintf(buf,"Euclidean and Lagrange sequences after using %ld-tap rectangular filter differ at %ld locations.\n",M,n);printofe(lptr,buf,time1);}
  //----------------------------------------------
  sprintf(buf,"Compare Euclidean and Lagrange mappings after %ld-tap Hanning filtering...\n  ",M);printofe(lptr,buf,time1);
  //----------------------------------------------
  n=cmp(&xR3he,&xR3hl,0,lptr);
  if(n==0l){sprintf(buf,"Euclidean and Lagrange sequences after using %ld-tap Hanning filter are identical.\n",M);printofe(lptr,buf,time1);}
  else     {sprintf(buf,"Euclidean and Lagrange sequences after using %ld-tap Hanning filter differ at %ld locations.\n",M,n);printofe(lptr,buf,time1);}

  //----------------------------------------------
  sprintf(buf,"Down sample by factor of %ld...\n  ",M);printofe(lptr,buf,time1);
  //----------------------------------------------
  printf("xR1re|%ld...",M); downsample(M, &xR1re, &xR1reM);//xR1reM = xR1re downsampled by M
  printf("xR1he|%ld...",M); downsample(M, &xR1he, &xR1heM);//xR1heM = xR1he downsampled by M
  printf("xR3re|%ld...",M); downsample(M, &xR3re, &xR3reM);//xR3reM = xR3re downsampled by M
  printf("xR3rl|%ld...",M); downsample(M, &xR3rl, &xR3rlM);//xR3rlM = xR3rl downsampled by M
  printf("xR3he|%ld...",M); downsample(M, &xR3he, &xR3heM);//xR3heM = xR3he downsampled by M
  printf("xR3hl|%ld...",M); downsample(M, &xR3hl, &xR3hlM);//xR3hlM = xR3hl downsampled by M
  sprintf(buf,"\ndone.");printofe(lptr,buf,time1);

  //----------------------------------------------
  sprintf(buf,"Remove beginning %ld element and ending %ld element transition regions... ",M-1,M);printofe(lptr,buf,time1);
  //----------------------------------------------
  for(n=0;n<N;n++){
    xR1reMN.put(n,xR1reM.get(n+1));
    xR1heMN.put(n,xR1heM.get(n+1));
    xR3reMN.put(n,xR3reM.get(n+1));
    xR3rlMN.put(n,xR3rlM.get(n+1));
    xR3heMN.put(n,xR3heM.get(n+1));
    xR3hlMN.put(n,xR3hlM.get(n+1));
    }
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  //----------------------------------------------
  sprintf(buf,"Perform auto-correlation operations...");printofe(lptr,buf,time1);
  //----------------------------------------------
  if(xR1reMN.Rxxo(&RxR1re,1)){fprintf(stderr,"ERROR computing Rxx for x=xR1re.\n"); return -1;} //auto-correlation seqR1 of xR1reMN
  if(xR1heMN.Rxxo(&RxR1he,1)){fprintf(stderr,"ERROR computing Rxx for x=xR1he.\n"); return -1;} //auto-correlation seqR1 of xR1heMN
  if(xR3reMN.Rxxo(&RxR3re,1)){fprintf(stderr,"ERROR computing Rxx for x=xR3re.\n"); return -1;} //auto-correlation seqR1 of xR3reMN
  if(xR3rlMN.Rxxo(&RxR3rl,1)){fprintf(stderr,"ERROR computing Rxx for x=xR3rl.\n"); return -1;} //auto-correlation seqR1 of xR3rlMN
  if(xR3heMN.Rxxo(&RxR3he,1)){fprintf(stderr,"ERROR computing Rxx for x=xR3he.\n"); return -1;} //auto-correlation seqR1 of xR3heMN
  if(xR3hlMN.Rxxo(&RxR3hl,1)){fprintf(stderr,"ERROR computing Rxx for x=xR3hl.\n"); return -1;} //auto-correlation seqR1 of xR3hlMN
  // |             |      |____________switch to turn on counting display
  // |             |___________________pointer to output correlation sequence
  // |_________________________________input real die sequence
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  //----------------------------------------------
  //plot seqR1 tex files
  //----------------------------------------------
  sprintf(comment,"length %ld real die SEQUENCE filtered in R^1 using %ld tap Rectangular filter and mapped backc to real die seqR1 using Euclidean metric",N,M);
  sprintf(filename,"%s_%ld_R1_rect%ld_euclid_seq.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_seq((symseq *)&xR1reM, 0,50, time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_seq(...)\n"); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  sprintf(comment,"length %ld real die SEQUENCE filtered in R^1 using %ld tap Hanning filter and mapped backc to real die seqR1 using Euclidean metric",N,M);
  sprintf(filename,"%s_%ld_R1_hann%ld_euclid_seq.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_seq((symseq *)&xR1heM, 0,50, time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_seq(...)\n"); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  sprintf(comment,"length %ld real die SEQUENCE filtered in R^3 using %ld tap rectangular filter and mapped back to real die seqR1 using Euclidean metric",N,M);
  sprintf(filename,"%s_%ld_R3_rect%ld_euclid_seq.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_seq((symseq *)&xR3reM, 0,50, time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_seq(...)\n"); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  sprintf(comment,"length %ld real die SEQUENCE filtered in R^3 using %ld tap rectangular filter and mapped back to real die seqR1 using Lagrange arc distance",N,M);
  sprintf(filename,"%s_%ld_R3_rect%ld_larc_seq.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_seq((symseq *)&xR3rlM, 0,50, time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_seq(...)\n"); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  sprintf(comment,"length %ld real die SEQUENCE filtered in R^3 using %ld tap Hanning filter and mapped back to real die seqR1 using Euclidean metric",N,M);
  sprintf(filename,"%s_%ld_R3_hann%ld_euclid_seq.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_seq((symseq *)&xR3heM, 0,50, time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_seq(...)\n"); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  sprintf(comment,"length %ld real die SEQUENCE filtered in R^3 using %ld tap Hanning filter and mapped back to real die seqR1 using Lagrange arc distance",N,M);
  sprintf(filename,"%s_%ld_R3_hann%ld_larc_seq.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_seq((symseq *)&xR3hlM, 0,50, time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_seq(...)\n"); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);


  //----------------------------------------------
  //plot histogram tex files
  //----------------------------------------------
  sprintf(comment,"HISTOGRAM of length %ld real die seqR1 filtered in R^1 using %ld tap Rectangular filter with each tap set to %.6lf",N,M,1.0);
  sprintf(filename,"%s_%ld_R1_rect%ld_euclid_histo.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_histo((symseq *)&xR1reMN, time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_histo(...)\n"); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  sprintf(comment,"HISTOGRAM of length %ld real die seqR1 filtered in R^1 using %ld tap Hanning filter",N,M);
  sprintf(filename,"%s_%ld_R1_hann%ld_euclid_histo.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_histo((symseq *)&xR1heMN, time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_histo(...)\n"); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  sprintf(comment,"HISTOGRAM of length %ld real die seqR1 filtered in R^3 using %ld tap rectangular filter and mapped back to real die seqR1 using Euclidean metric",N,M);
  sprintf(filename,"%s_%ld_R3_rect%ld_euclid_histo.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_histo((symseq *)&xR3reMN, time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_histo(...)\n"); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  sprintf(comment,"HISTOGRAM of length %ld real die seqR1 filtered in R^3 using %ld tap rectangular filter and mapped back to real die seqR1 using Lagrange arc distance",N,M);
  sprintf(filename,"%s_%ld_R3_rect%ld_larc_histo.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_histo((symseq *)&xR3rlMN, time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_histo(...)\n"); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  sprintf(comment,"HISTOGRAM of length %ld real die seqR1 filtered in R^3 using %ld tap Hanning filter and mapped back to real die seqR1 using Euclidean metric",N,M);
  sprintf(filename,"%s_%ld_R3_hann%ld_euclid_histo.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_histo((symseq *)&xR3heMN, time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_histo(...)\n"); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  sprintf(comment,"HISTOGRAM of length %ld real die seqR1 filtered in R^3 using %ld tap Hanning filter and mapped back to real die seqR1 using Lagrange arc distance",N,M);
  sprintf(filename,"%s_%ld_R3_hann%ld_larc_histo.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_histo((symseq *)&xR3hlMN, time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_histo(...)\n"); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  //----------------------------------------------
  //plot autocorrelation tex files
  //----------------------------------------------
  sprintf(comment,"AUTO-CORRELATION of length %ld real die seqR1 filtered in R^1 using %ld tap Rectangular filter with each tap set to %.6lf",N,M,1.0);
  sprintf(filename,"%s_%ld_R1_rect%ld_euclid_auto.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_auto(&RxR1re, M, time1, filename,comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_auto(...,%s,...)\n",filename); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  sprintf(comment,"AUTO-CORRELATION of length %ld real die seqR1 filtered in R^1 using %ld tap Hanning filter",N,M);
  sprintf(filename,"%s_%ld_R1_hann%ld_euclid_auto.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_auto(&RxR1he, M, time1, filename,comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_auto(...,%s,...)\n",filename); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  sprintf(comment,"AUTO-CORRELATION of length %ld real die seqR1 filtered in R^3 using %ld tap rectangular filter and mapped back to real die seqR1 using Euclidean metric",N,M);
  sprintf(filename,"%s_%ld_R3_rect%ld_euclid_auto.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_auto(&RxR3re, M, time1, filename,comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_auto(...,%s,...)\n",filename); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  sprintf(comment,"AUTO-CORRELATION of length %ld real die seqR1 filtered in R^e using %ld tap rectangular filter and mapped back to real die seqR1 using Lagrange arc distance",N,M);
  sprintf(filename,"%s_%ld_R3_rect%ld_larc_auto.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_auto(&RxR3rl, M, time1, filename,comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_auto(...,%s,...)\n",filename); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  sprintf(comment,"AUTO-CORRELATION of length %ld real die seqR1 filtered in R^3 using %ld tap Hanning filter and mapped back to real die seqR1 using Euclidean metric",N,M);
  sprintf(filename,"%s_%ld_R3_hann%ld_euclid_auto.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_auto(&RxR3he, M, time1, filename,comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_auto(...,%s,...)\n",filename); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  sprintf(comment,"AUTO-CORRELATION of length %ld real die seqR1 filtered in R^3 using %ld tap Hanning filter and mapped back to real die seqR1 using Lagrange arc distance",N,M);
  sprintf(filename,"%s_%ld_R3_hann%ld_larc_auto.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_auto(&RxR3hl, M, time1, filename,comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_auto(...,%s,...)\n",filename); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  //----------------------------------------------
  // close log file
  //----------------------------------------------
  plot_close(lptr,time1);
  return 0;
  }



//-----------------------------------------------------------------------------
//! \brief generate TeX plot file for the low pass filtering of spinner sequence
//!        filtered by M tap rectangular filter,
//!        and generate TeX plot file
//-----------------------------------------------------------------------------
int lab_spin_lp(const unsigned seed, const long N, const long M, const char *basefilename){
  time_t time1; time(&time1);   //starting time stamp (passed to plotting routine)
  spinseq    x(N+M);        //for die sequence
  seqR1   xR1(N+M);      //for R^1 seqR1 mapped from x sequence
  seqR2 xR2(N+M);      //for R^2 seqR1 mapped from x sequence
  seqR1   r(M);          //for M-tap rectangular filter
  seqR1   h(M);          //for M-tap Hanning filter
  seqR1   xR1r(N+M+M-1); //for x mapped to R1 and filtered by r
  seqR1   xR1h(N+M+M-1); //for x mapped to R1 and filtered by h
  seqR2 xR2r(N+M+M-1); //for x mapped to R2 and filtered by r
  seqR2 xR2h(N+M+M-1); //for x mapped to R2 and filtered by h
               //     |  |_______________M-1: convolution with h appends M-1 elements
               //     |__________________ Nm: length of xR1 sequence
  spinseq xR1re(N+M+M-1);   //for xR1r mapped back using Euclidean metric
  spinseq xR1he(N+M+M-1);   //for xR1h mapped back using Euclidean metric
  spinseq xR2re(N+M+M-1);   //for xR2r mapped back using Euclidean metric
  spinseq xR2rl(N+M+M-1);   //for xR2r mapped back using Lagrange arc distance
  spinseq xR2he(N+M+M-1);   //for xR2h mapped back using Euclidean metric
  spinseq xR2hl(N+M+M-1);   //for xR2h mapped back using Lagrange arc distance
                           
  spinseq xR1reN(N);       //for xR1re with beginning M-1 and ending M-1 elements removed
  spinseq xR1heN(N);       //for xR1he with beginning M-1 and ending M-1 elements removed
  spinseq xR2reN(N);       //for xR2re with beginning M-1 and ending M-1 elements removed
  spinseq xR2rlN(N);       //for xR2rl with beginning M-1 and ending M-1 elements removed
  spinseq xR2heN(N);       //for xR2he with beginning M-1 and ending M-1 elements removed
  spinseq xR2hlN(N);       //for xR2hl with beginning M-1 and ending M-1 elements removed
                           
  seqR1 RxR1re(2*N+1);  //for auto-correlation seqR1 of xR1re 
  seqR1 RxR1he(2*N+1);  //for auto-correlation seqR1 of xR1he 
  seqR1 RxR2re(2*N+1);  //for auto-correlation seqR1 of xR2re 
  seqR1 RxR2rl(2*N+1);  //for auto-correlation seqR1 of xR2rl 
  seqR1 RxR2he(2*N+1);  //for auto-correlation seqR1 of xR2he 
  seqR1 RxR2hl(2*N+1);  //for auto-correlation seqR1 of xR2hl 

  char comment[2*1024];          //comment to be passed to plotting function
  char buf[2*1024];          //comment to be passed to plotting function
  long n;
  char filename[1024];
  FILE *lptr; // pointer to log  file

  //----------------------------------------------
  //open log file
  //----------------------------------------------
  printf(         "-----------------------------------------------------------\n");
  sprintf(comment,"map nonstationary length %ld die sequence\n%%  with length %ld periodic distribution into R^1 using PAM mapping\n%%  and use DFT:R^1-->C^1 to analyze",N,M); printf("%s",comment);
  printf(         "-----------------------------------------------------------\n");
  sprintf(buf,"%s_%ldm%ld",basefilename,N,M);
  lptr=log_open (buf,time1,comment);

  //----------------------------------------------
  printf("Generate random seqR1 using seed value 0x%4x...",seed);
  //----------------------------------------------
  x.randomize(seed);
  printf("done.\n");

  //----------------------------------------------
  printf("Perform mapping operations from spinner space to R^n...");
  //----------------------------------------------
  xR1=x.spintoR1();               //xR1 = x seqR1 mapped to R^1 sequence
  xR2=x.spintoR2();               //xR2 = x die seqR1 mapped to R^2 sequence
  printf("done.\n");

  //----------------------------------------------
  printf("Perform filtering operations...");
  //----------------------------------------------
  r=1.0;       r.normalize(); putchar('\n');r.list();putchar('\n');
  h.hanning(); h.normalize(); putchar('\n');h.list();putchar('\n'); 
  printf("xR1*r..."); convolve(&xR1, &r, &xR1r); //xR1r=xR1*r; 
  printf("xR1*h..."); convolve(&xR1, &h, &xR1h); //xR1h=xR1*h; 
  printf("xR3*r..."); convolve(&xR2, &r, &xR2r); //xR2r=xR2*r; 
  printf("xR3*h..."); convolve(&xR2, &h, &xR2h); //xR2h=xR2*h; 
  printf("done.\n");
  //----------------------------------------------
  printf("Perform mapping operations from R^n back to die seqR1 space...");
  //----------------------------------------------
  xR1re=spin_R1tospin_euclid(xR1r);//xR1r mapped back to spinner seq. using Euclid
  xR1he=spin_R1tospin_euclid(xR1h);//xR1h mapped back to spinner seq. using Euclid
  xR2re=spin_R2tospin_euclid(xR2r);//xR2r mapped back to spinner seq. using Euclid
  xR2rl=spin_R2tospin_larc  (xR2r);//xR2r mapped back to spinner seq. using Lagrange
  xR2he=spin_R2tospin_euclid(xR2h);//xR2h mapped back to spinner seq. using Euclid
  xR2hl=spin_R2tospin_larc  (xR2h);//xR2h mapped back to spinner seq. using Lagrange
  printf("done.\n");

  //----------------------------------------------
  sprintf(buf,"Compare Euclidean and Lagrange mappings after %ld-tap rectangular filtering...\n  ",M);printofe(lptr,buf,time1);
  //----------------------------------------------
  n=cmp(&xR2re,&xR2rl,0,lptr);
  if(n==0l){sprintf(buf,"Euclidean and Lagrange sequences following %ld-tap rectangular filter are identical.\n",M);printofe(lptr,buf,time1);}
  else     {sprintf(buf,"Euclidean and Lagrange sequences following %ld-tap rectangular filter differ at %ld locations.\n",M,n);printofe(lptr,buf,time1);}
  //----------------------------------------------
  sprintf(buf,"Compare Euclidean and Lagrange mappings after %ld-tap Hanning filtering...\n  ",M);printofe(lptr,buf,time1);
  //----------------------------------------------
  n=cmp(&xR2he,&xR2hl,0,lptr);
  if(n==0l){sprintf(buf,"Euclidean and Lagrange sequences following %ld-tap Hanning filter are identical.\n",M);printofe(lptr,buf,time1);}
  else     {sprintf(buf,"Euclidean and Lagrange sequences following %ld-tap Hanning filter differ at %ld locations.\n",M,n);printofe(lptr,buf,time1);}

  //----------------------------------------------
  printf("Remove beginning %ld element and ending %ld element transition regions introduced by filtering ... ",M-1,M);
  //----------------------------------------------
  for(n=0;n<N;n++){
    xR1reN.put(n,xR1re.get(n+M-1));
    xR1heN.put(n,xR1he.get(n+M-1));
    xR2reN.put(n,xR2re.get(n+M-1));
    xR2rlN.put(n,xR2rl.get(n+M-1));
    xR2heN.put(n,xR2he.get(n+M-1));
    xR2hlN.put(n,xR2hl.get(n+M-1));
    }
  printf("done.\n");

  //----------------------------------------------
  printf("Perform auto-correlation operations...\n");
  //----------------------------------------------
  if(xR1reN.Rxxo(&RxR1re,1)){fprintf(stderr,"ERROR computing Rxx for x=xR1re.\n"); return -1;} //auto-correlation seqR1 of truncated xR1re
  if(xR1heN.Rxxo(&RxR1he,1)){fprintf(stderr,"ERROR computing Rxx for x=xR1he.\n"); return -1;} //auto-correlation seqR1 of truncated xR1he
  if(xR2reN.Rxxo(&RxR2re,1)){fprintf(stderr,"ERROR computing Rxx for x=xR2re.\n"); return -1;} //auto-correlation seqR1 of truncated xR2re
  if(xR2rlN.Rxxo(&RxR2rl,1)){fprintf(stderr,"ERROR computing Rxx for x=xR2rl.\n"); return -1;} //auto-correlation seqR1 of truncated xR2rl
  if(xR2heN.Rxxo(&RxR2he,1)){fprintf(stderr,"ERROR computing Rxx for x=xR2he.\n"); return -1;} //auto-correlation seqR1 of truncated xR2he
  if(xR2hlN.Rxxo(&RxR2hl,1)){fprintf(stderr,"ERROR computing Rxx for x=xR2hl.\n"); return -1;} //auto-correlation seqR1 of truncated xR2hl
  // |          |      |____________switch to turn on counting display
  // |          |___________________pointer to output correlation sequence
  // |______________________________input spinner sequence
  printf("done.\n");

  //----------------------------------------------
  //plot seqR1 tex files
  //----------------------------------------------
  sprintf(comment,"length %ld spinner SEQUENCE filtered in R^1 using %ld tap Rectangular filter and mapped backc to spinner seqR1 using Euclidean metric",N,M);
  sprintf(filename,"%s_%ld_R1_rect%ld_euclid_seq.tex",basefilename,N,M);
  printf("Plot file %50s ... ",filename);
  if(plot_ocs_seq((symseq *)&xR1re, M-1, M-1+48,time1, "spin", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_seq(...)\n"); return -1;}
  printf("done.\n");

  sprintf(comment,"length %ld spinner SEQUENCE filtered in R^1 using %ld tap Hanning filter and mapped backc to spinner seqR1 using Euclidean metric",N,M);
  sprintf(filename,"%s_%ld_R1_hann%ld_euclid_seq.tex",basefilename,N,M);
  printf("Plot file %50s ... ",filename);
  if(plot_ocs_seq((symseq *)&xR1he, M-1, M-1+48,time1, "spin", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_seq(...)\n"); return -1;}
  printf("done.\n");

  sprintf(comment,"length %ld spinner SEQUENCE filtered in R^2 using %ld tap rectangular filter and mapped back to spinner seqR1 using Euclidean metric",N,M);
  sprintf(filename,"%s_%ld_R2_rect%ld_euclid_seq.tex",basefilename,N,M);
  printf("Plot file %50s ... ",filename);
  if(plot_ocs_seq((symseq *)&xR2re, M-1, M-1+48,time1, "spin", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_seq(...)\n"); return -1;}
  printf("done.\n");

  sprintf(comment,"length %ld spinner SEQUENCE filtered in R^2 using %ld tap rectangular filter and mapped back to spinner seqR1 using Lagrange arc distance",N,M);
  sprintf(filename,"%s_%ld_R2_rect%ld_larc_seq.tex",basefilename,N,M);
  printf("Plot file %50s ... ",filename);
  if(plot_ocs_seq((symseq *)&xR2rl, M-1, M-1+48,time1, "spin", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_seq(...)\n"); return -1;}
  printf("done.\n");

  sprintf(comment,"length %ld spinner SEQUENCE filtered in R^2 using %ld tap Hanning filter and mapped back to spinner seqR1 using Euclidean metric",N,M);
  sprintf(filename,"%s_%ld_R2_hann%ld_euclid_seq.tex",basefilename,N,M);
  printf("Plot file %50s ... ",filename);
  if(plot_ocs_seq((symseq *)&xR2he, M-1, M-1+48,time1, "spin", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_seq(...)\n"); return -1;}
  printf("done.\n");

  sprintf(comment,"length %ld spinner SEQUENCE filtered in R^2 using %ld tap Hanning filter and mapped back to spinner seqR1 using Lagrange arc distance",N,M);
  sprintf(filename,"%s_%ld_R2_hann%ld_larc_seq.tex",basefilename,N,M);
  printf("Plot file %50s ... ",filename);
  if(plot_ocs_seq((symseq *)&xR2hl, M-1, M-1+48,time1, "spin", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_seq(...)\n"); return -1;}
  printf("done.\n");

  //----------------------------------------------
  //plot histogram tex files
  //----------------------------------------------
  sprintf(comment,"HISTOGRAM of length %ld spinner seqR1 filtered in R^1 using %ld tap Rectangular filter ",N,M);
  sprintf(filename,"%s_%ld_R1_rect%ld_euclid_histo.tex",basefilename,N,M);
  printf("Plot file %50s ... ",filename);
  if(plot_ocs_histo((symseq *)&xR1re, M-1, N+M-2, time1, "spin", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_histo(...)\n"); return -1;}
  printf("done.\n");

  sprintf(comment,"HISTOGRAM of length %ld spinner seqR1 filtered in R^1 using %ld tap Hanning filter",N,M);
  sprintf(filename,"%s_%ld_R1_hann%ld_euclid_histo.tex",basefilename,N,M);
  printf("Plot file %50s ... ",filename);
  if(plot_ocs_histo((symseq *)&xR1he, M-1, N+M-2, time1, "spin", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_histo(...)\n"); return -1;}
  printf("done.\n");

  sprintf(comment,"HISTOGRAM of length %ld spinner seqR1 filtered in R^2 using %ld tap rectangular filter and mapped back to spinner seqR1 using Euclidean metric",N,M);
  sprintf(filename,"%s_%ld_R2_rect%ld_euclid_histo.tex",basefilename,N,M);
  printf("Plot file %50s ... ",filename);
  if(plot_ocs_histo((symseq *)&xR2re, M-1, N+M-2, time1, "spin", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_histo(...)\n"); return -1;}
  printf("done.\n");

  sprintf(comment,"HISTOGRAM of length %ld spinner seqR1 filtered in R^2 using %ld tap rectangular filter and mapped back to spinner seqR1 using Lagrange arc distance",N,M);
  sprintf(filename,"%s_%ld_R2_rect%ld_larc_histo.tex",basefilename,N,M);
  printf("Plot file %50s ... ",filename);
  if(plot_ocs_histo((symseq *)&xR2rl, M-1, N+M-2, time1, "spin", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_histo(...)\n"); return -1;}
  printf("done.\n");

  sprintf(comment,"HISTOGRAM of length %ld spinner seqR1 filtered in R^2 using %ld tap Hanning filter and mapped back to spinner seqR1 using Euclidean metric",N,M);
  sprintf(filename,"%s_%ld_R2_hann%ld_euclid_histo.tex",basefilename,N,M);
  printf("Plot file %50s ... ",filename);
  if(plot_ocs_histo((symseq *)&xR2he, M-1, N+M-2, time1, "spin", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_histo(...)\n"); return -1;}
  printf("done.\n");

  sprintf(comment,"HISTOGRAM of length %ld spinner seqR1 filtered in R^2 using %ld tap Hanning filter and mapped back to spinner seqR1 using Lagrange arc distance",N,M);
  sprintf(filename,"%s_%ld_R2_hann%ld_larc_histo.tex",basefilename,N,M);
  printf("Plot file %50s ... ",filename);
  if(plot_ocs_histo((symseq *)&xR2hl, M-1, N+M-2, time1, "spin", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_histo(...)\n"); return -1;}
  printf("done.\n");

  //----------------------------------------------
  //plot autocorrelation tex files
  //----------------------------------------------
  sprintf(comment,"AUTO-CORRELATION of length %ld spinner seqR1 filtered in R^1 using %ld tap Rectangular filter ",N,M);
  sprintf(filename,"%s_%ld_R1_rect%ld_euclid_auto.tex",basefilename,N,M);
  printf("Plot file %50s ... ",filename);
  if(plot_ocs_auto(&RxR1re, M, time1, filename,comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_auto(...,%s,...)\n",filename); return -1;}
  printf("done.\n");

  sprintf(comment,"AUTO-CORRELATION of length %ld spinner seqR1 filtered in R^1 using %ld tap Hanning filter",N,M);
  sprintf(filename,"%s_%ld_R1_hann%ld_euclid_auto.tex",basefilename,N,M);
  printf("Plot file %50s ... ",filename);
  if(plot_ocs_auto(&RxR1he, M, time1, filename,comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_auto(...,%s,...)\n",filename); return -1;}
  printf("done.\n");

  sprintf(comment,"AUTO-CORRELATION of length %ld spinner seqR1 filtered in R^2 using %ld tap rectangular filter and mapped back to spinner seqR1 using Euclidean metric",N,M);
  sprintf(filename,"%s_%ld_R2_rect%ld_euclid_auto.tex",basefilename,N,M);
  printf("Plot file %50s ... ",filename);
  if(plot_ocs_auto(&RxR2re, M, time1, filename,comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_auto(...,%s,...)\n",filename); return -1;}
  printf("done.\n");

  sprintf(comment,"AUTO-CORRELATION of length %ld spinner seqR1 filtered in R^e using %ld tap rectangular filter and mapped back to spinner seqR1 using Lagrange arc distance",N,M);
  sprintf(filename,"%s_%ld_R2_rect%ld_larc_auto.tex",basefilename,N,M);
  printf("Plot file %50s ... ",filename);
  if(plot_ocs_auto(&RxR2rl, M, time1, filename,comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_auto(...,%s,...)\n",filename); return -1;}
  printf("done.\n");

  sprintf(comment,"AUTO-CORRELATION of length %ld spinner seqR1 filtered in R^2 using %ld tap Hanning filter and mapped back to spinner seqR1 using Euclidean metric",N,M);
  sprintf(filename,"%s_%ld_R2_hann%ld_euclid_auto.tex",basefilename,N,M);
  printf("Plot file %50s ... ",filename);
  if(plot_ocs_auto(&RxR2he, M, time1, filename,comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_auto(...,%s,...)\n",filename); return -1;}
  printf("done.\n");

  sprintf(comment,"AUTO-CORRELATION of length %ld spinner seqR1 filtered in R^2 using %ld tap Hanning filter and mapped back to spinner seqR1 using Lagrange arc distance",N,M);
  sprintf(filename,"%s_%ld_R2_hann%ld_larc_auto.tex",basefilename,N,M);
  printf("Plot file %50s ... ",filename);
  if(plot_ocs_auto(&RxR2hl, M, time1, filename,comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_auto(...,%s,...)\n",filename); return -1;}
  printf("done.\n");

  //----------------------------------------------
  // close log file
  //----------------------------------------------
  plot_close(lptr,time1);
  return 0;
  }

//-----------------------------------------------------------------------------
//! \brief generate TeX plot file for the high pass filtering of weighted real die sequence
//!        filtered by M tap rectangular filter scaled to normalized values,
//!        and generate TeX plot file
//-----------------------------------------------------------------------------
int lab_wspin_hp(const unsigned seed, const long N, const long M, const char *basefilename){
  time_t time1; time(&time1);   //starting time stamp (passed to plotting routine)
  spinseq x(M*(N+2)-(M-1));   //for spinner sequence
  seqR1 r(M);              //for M-tap rectangular filter
  seqR1 h(M);              //for M-tap Hanning filter

  seqR1    xR1(M*(N+2)-(M-1)); //for x mapped to R^1
  seqR2    xR2(M*(N+2)-(M-1)); //for x mapped to R^2

  seqR1    xR1r(M*(N+2)); //for xR1 filtered by r
  seqR1    xR1h(M*(N+2)); //for xR1 filtered by h
  seqR2    xR2r(M*(N+2)); //for xR2 filtered by r
  seqR2    xR2h(M*(N+2)); //for xR2 filtered by h

  spinseq  xR1re(M*(N+2));//for xR1r mapped back to die seqR1 using Euclidean metric
  spinseq  xR1he(M*(N+2));//for xR1h mapped back to die seqR1 using Euclidean metric
  spinseq  xR2re(M*(N+2));//for xR2r mapped back to die seqR1 using Euclidean metric
  spinseq  xR2he(M*(N+2));//for xR2h mapped back to die seqR1 using Euclidean metric
  spinseq  xR2rl(M*(N+2));//for xR2r mapped back to die seqR1 using Lagrange  distance
  spinseq  xR2hl(M*(N+2));//for xR2h mapped back to die seqR1 using Lagrange  distance
           
  spinseq  xR1reM(N+2);   //for xR1re downsampled by factor of M
  spinseq  xR1heM(N+2);   //for xR1he downsampled by factor of M
  spinseq  xR2reM(N+2);   //for xR2re downsampled by factor of M
  spinseq  xR2heM(N+2);   //for xR2he downsampled by factor of M
  spinseq  xR2rlM(N+2);   //for xR2rl downsampled by factor of M
  spinseq  xR2hlM(N+2);   //for xR2hl downsampled by factor of M
           
  spinseq  xR1reMN(N);    //for xR1reM with first and last element removed
  spinseq  xR1heMN(N);    //for xR1heM with first and last element removed
  spinseq  xR2reMN(N);    //for xR2reM with first and last element removed
  spinseq  xR2heMN(N);    //for xR2heM with first and last element removed
  spinseq  xR2rlMN(N);    //for xR2rlM with first and last element removed
  spinseq  xR2hlMN(N);    //for xR2hlM with first and last element removed

  seqR1    RxR1re(2*N+1); //for auto-correlation of xR1reMN
  seqR1    RxR1he(2*N+1); //for auto-correlation of xR1heMN
  seqR1    RxR2re(2*N+1); //for auto-correlation of xR2reMN
  seqR1    RxR2he(2*N+1); //for auto-correlation of xR2heMN
  seqR1    RxR2rl(2*N+1); //for auto-correlation of xR2rlMN
  seqR1    RxR2hl(2*N+1); //for auto-correlation of xR2hlMN

  char comment[2*1024];          //comment to be passed to plotting function
  char buf[2*1024];          //general purpose buffer
  long n;
  char filename[1024];

  FILE *lptr; // pointer to log  file

  //----------------------------------------------
  //open log file
  //----------------------------------------------
  printf(         "-----------------------------------------------------------\n");
  sprintf(comment,"Experiment: high pass filter length %ld weighted spinner sequence\n",N); printf("%s",comment);
  printf(         "-----------------------------------------------------------\n");
  sprintf(buf,"%s_%ldm%ld",basefilename,N,M);
  lptr=log_open (buf,time1,comment);

  //----------------------------------------------
  sprintf(buf,"Generate weighted die sequence...");printofe(lptr,buf,time1);
  //----------------------------------------------
  x.randomize(seed,5,5,5,5,75,5);
  x.histogram(1,lptr);
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  //----------------------------------------------
  sprintf(buf,"Map from die space ");printofe(lptr,buf,time1);
  //----------------------------------------------
  printf("to R^1..."); xR1=x.spintoR1();  //xR1 = x seqR1 mapped to R^1 sequence
  printf("to R^2..."); xR2=x.spintoR2();  //xR2 = x die seqR1 mapped to R^2 sequence
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  //----------------------------------------------
  sprintf(buf,"Generate filter coefficients...");printofe(lptr,buf,time1);
  //----------------------------------------------
  r=1.0;       r.lptohp(); printf("\n  r_n=");r.list();
  h.hanning(); h.lptohp(); printf("\n  h_n=");h.list();putchar('\n'); 
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  //----------------------------------------------
  sprintf(buf,"Perform filtering operations...");printofe(lptr,buf,time1);
  //----------------------------------------------
  printf("xR1*r..."); convolve(&xR1, &r, &xR1r); //xR1r=xR1*r; 
  printf("xR1*h..."); convolve(&xR1, &h, &xR1h); //xR1h=xR1*h; 
  printf("xR3*r..."); convolve(&xR2, &r, &xR2r); //xR2r=xR2*r; 
  printf("xR3*h..."); convolve(&xR2, &h, &xR2h); //xR2h=xR2*h; 
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  //----------------------------------------------
  sprintf(buf,"Perform mapping operations from R^n back to die seqR1 space...");printofe(lptr,buf,time1);
  //----------------------------------------------
  xR1re=spin_R1tospin_euclid(xR1r);//xR1r mapped back to spinner seq. using Euclid
  xR1he=spin_R1tospin_euclid(xR1h);//xR1h mapped back to spinner seq. using Euclid
  xR2re=spin_R2tospin_euclid(xR2r);//xR2r mapped back to spinner seq. using Euclid
  xR2rl=spin_R2tospin_larc  (xR2r);//xR2r mapped back to spinner seq. using Lagrange
  xR2he=spin_R2tospin_euclid(xR2h);//xR2h mapped back to spinner seq. using Euclid
  xR2hl=spin_R2tospin_larc  (xR2h);//xR2h mapped back to spinner seq. using Lagrange
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  //----------------------------------------------
  sprintf(buf,"Compare Euclidean and Lagrange mappings after %ld-tap rectangular filtering...\n  ",M);printofe(lptr,buf,time1);
  //----------------------------------------------
  n=cmp(&xR2re,&xR2rl,0,lptr);
  if(n==0l){sprintf(buf,"Euclidean and Lagrange sequences following %ld-tap rectangular filter are identical.\n",M);printofe(lptr,buf,time1);}
  else     {sprintf(buf,"Euclidean and Lagrange sequences following %ld-tap rectangular filter differ at %ld locations.\n",M,n);printofe(lptr,buf,time1);}
  //----------------------------------------------
  sprintf(buf,"Compare Euclidean and Lagrange mappings after %ld-tap Hanning filtering...\n  ",M);printofe(lptr,buf,time1);
  //----------------------------------------------
  n=cmp(&xR2he,&xR2hl,0,lptr);
  if(n==0l){sprintf(buf,"Euclidean and Lagrange sequences following %ld-tap Hanning filter are identical.\n",M);printofe(lptr,buf,time1);}
  else     {sprintf(buf,"Euclidean and Lagrange sequences following %ld-tap Hanning filter differ at %ld locations.\n",M,n);printofe(lptr,buf,time1);}

  //----------------------------------------------
  sprintf(buf,"Perform down sampling...\n  ");printofe(lptr,buf,time1);
  //----------------------------------------------
  printf("xR1re|%ld...",M); downsample(M, &xR1re, &xR1reM);//xR1reM = xR1re downsampled by M
  printf("xR1he|%ld...",M); downsample(M, &xR1he, &xR1heM);//xR1heM = xR1he downsampled by M
  printf("xR2re|%ld...",M); downsample(M, &xR2re, &xR2reM);//xR2reM = xR2re downsampled by M
  printf("xR2rl|%ld...",M); downsample(M, &xR2rl, &xR2rlM);//xR2rlM = xR2rl downsampled by M
  printf("xR2he|%ld...",M); downsample(M, &xR2he, &xR2heM);//xR2heM = xR2he downsampled by M
  printf("xR2hl|%ld...",M); downsample(M, &xR2hl, &xR2hlM);//xR2hlM = xR2hl downsampled by M
  sprintf(buf,"\ndone.");printofe(lptr,buf,time1);

  //----------------------------------------------
  sprintf(buf,"Remove beginning %ld element and ending %ld element transition regions... ",M-1,M);printofe(lptr,buf,time1);
  //----------------------------------------------
  for(n=0;n<N;n++){
    xR1reMN.put(n,xR1reM.get(n+1));
    xR1heMN.put(n,xR1heM.get(n+1));
    xR2reMN.put(n,xR2reM.get(n+1));
    xR2rlMN.put(n,xR2rlM.get(n+1));
    xR2heMN.put(n,xR2heM.get(n+1));
    xR2hlMN.put(n,xR2hlM.get(n+1));
    }
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  //----------------------------------------------
  sprintf(buf,"Perform auto-correlation operations...");printofe(lptr,buf,time1);
  //----------------------------------------------
  if(xR1reMN.Rxxo(&RxR1re,1)){fprintf(stderr,"ERROR computing Rxx for x=xR1reMN.\n"); return -1;} //auto-correlation seqR1 of xR1reMN
  if(xR1heMN.Rxxo(&RxR1he,1)){fprintf(stderr,"ERROR computing Rxx for x=xR1heMN.\n"); return -1;} //auto-correlation seqR1 of xR1heMN
  if(xR2reMN.Rxxo(&RxR2re,1)){fprintf(stderr,"ERROR computing Rxx for x=xR2reMN.\n"); return -1;} //auto-correlation seqR1 of xR2reMN
  if(xR2rlMN.Rxxo(&RxR2rl,1)){fprintf(stderr,"ERROR computing Rxx for x=xR2rlMN.\n"); return -1;} //auto-correlation seqR1 of xR2rlMN
  if(xR2heMN.Rxxo(&RxR2he,1)){fprintf(stderr,"ERROR computing Rxx for x=xR2heMN.\n"); return -1;} //auto-correlation seqR1 of xR2heMN
  if(xR2hlMN.Rxxo(&RxR2hl,1)){fprintf(stderr,"ERROR computing Rxx for x=xR2hlMN.\n"); return -1;} //auto-correlation seqR1 of xR2hlMN
  // |            |      |____________switch to turn on counting display
  // |            |___________________pointer to output correlation sequence
  // |________________________________input spinner sequence
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  //----------------------------------------------
  //plot seqR1 tex files
  //----------------------------------------------
  sprintf(comment,"length %ld spinner SEQUENCE filtered in R^1 using %ld tap Rectangular filter and mapped backc to spinner seqR1 using Euclidean metric",N,M);
  sprintf(filename,"%s_%ld_R1_rect%ld_euclid_seq.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_seq((symseq *)&xR1reM, 0,50, time1, "spin", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_seq(...)\n"); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  sprintf(comment,"length %ld spinner SEQUENCE filtered in R^1 using %ld tap Hanning filter and mapped backc to spinner seqR1 using Euclidean metric",N,M);
  sprintf(filename,"%s_%ld_R1_hann%ld_euclid_seq.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_seq((symseq *)&xR1heM, 0,50, time1, "spin", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_seq(...)\n"); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  sprintf(comment,"length %ld spinner SEQUENCE filtered in R^2 using %ld tap rectangular filter and mapped back to spinner seqR1 using Euclidean metric",N,M);
  sprintf(filename,"%s_%ld_R2_rect%ld_euclid_seq.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_seq((symseq *)&xR2reM, 0,50, time1, "spin", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_seq(...)\n"); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  sprintf(comment,"length %ld spinner SEQUENCE filtered in R^2 using %ld tap rectangular filter and mapped back to spinner seqR1 using Lagrange arc distance",N,M);
  sprintf(filename,"%s_%ld_R2_rect%ld_larc_seq.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_seq((symseq *)&xR2rlM, 0,50, time1, "spin", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_seq(...)\n"); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  sprintf(comment,"length %ld spinner SEQUENCE filtered in R^2 using %ld tap Hanning filter and mapped back to spinner seqR1 using Euclidean metric",N,M);
  sprintf(filename,"%s_%ld_R2_hann%ld_euclid_seq.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_seq((symseq *)&xR2heM, 0,50, time1, "spin", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_seq(...)\n"); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  sprintf(comment,"length %ld spinner SEQUENCE filtered in R^2 using %ld tap Hanning filter and mapped back to spinner seqR1 using Lagrange arc distance",N,M);
  sprintf(filename,"%s_%ld_R2_hann%ld_larc_seq.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_seq((symseq *)&xR2hlM, 0,50, time1, "spin", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_seq(...)\n"); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);


  //----------------------------------------------
  //plot histogram tex files
  //----------------------------------------------
  sprintf(comment,"HISTOGRAM of length %ld spinner seqR1 filtered in R^1 using %ld tap Rectangular filter ",N,M);
  sprintf(filename,"%s_%ld_R1_rect%ld_euclid_histo.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_histo((symseq *)&xR1reMN, time1, "spin", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_histo(...)\n"); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  sprintf(comment,"HISTOGRAM of length %ld spinner seqR1 filtered in R^1 using %ld tap Hanning filter",N,M);
  sprintf(filename,"%s_%ld_R1_hann%ld_euclid_histo.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_histo((symseq *)&xR1heMN, time1, "spin", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_histo(...)\n"); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  sprintf(comment,"HISTOGRAM of length %ld spinner seqR1 filtered in R^2 using %ld tap rectangular filter and mapped back to spinner seqR1 using Euclidean metric",N,M);
  sprintf(filename,"%s_%ld_R2_rect%ld_euclid_histo.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_histo((symseq *)&xR2reMN, time1, "spin", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_histo(...)\n"); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  sprintf(comment,"HISTOGRAM of length %ld spinner seqR1 filtered in R^2 using %ld tap rectangular filter and mapped back to spinner seqR1 using Lagrange arc distance",N,M);
  sprintf(filename,"%s_%ld_R2_rect%ld_larc_histo.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_histo((symseq *)&xR2rlMN, time1, "spin", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_histo(...)\n"); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  sprintf(comment,"HISTOGRAM of length %ld spinner seqR1 filtered in R^2 using %ld tap Hanning filter and mapped back to spinner seqR1 using Euclidean metric",N,M);
  sprintf(filename,"%s_%ld_R2_hann%ld_euclid_histo.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_histo((symseq *)&xR2heMN, time1, "spin", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_histo(...)\n"); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  sprintf(comment,"HISTOGRAM of length %ld spinner seqR1 filtered in R^2 using %ld tap Hanning filter and mapped back to spinner seqR1 using Lagrange arc distance",N,M);
  sprintf(filename,"%s_%ld_R2_hann%ld_larc_histo.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_histo((symseq *)&xR2hlMN, time1, "spin", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_histo(...)\n"); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  //----------------------------------------------
  //plot autocorrelation tex files
  //----------------------------------------------
  sprintf(comment,"AUTO-CORRELATION of length %ld spinner seqR1 filtered in R^1 using %ld tap Rectangular filter",N,M);
  sprintf(filename,"%s_%ld_R1_rect%ld_euclid_auto.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_auto(&RxR1re, M, time1, filename,comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_auto(...,%s,...)\n",filename); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  sprintf(comment,"AUTO-CORRELATION of length %ld spinner seqR1 filtered in R^1 using %ld tap Hanning filter",N,M);
  sprintf(filename,"%s_%ld_R1_hann%ld_euclid_auto.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_auto(&RxR1he, M, time1, filename,comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_auto(...,%s,...)\n",filename); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  sprintf(comment,"AUTO-CORRELATION of length %ld spinner seqR1 filtered in R^2 using %ld tap rectangular filter and mapped back to spinner seqR1 using Euclidean metric",N,M);
  sprintf(filename,"%s_%ld_R2_rect%ld_euclid_auto.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_auto(&RxR2re, M, time1, filename,comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_auto(...,%s,...)\n",filename); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  sprintf(comment,"AUTO-CORRELATION of length %ld spinner seqR1 filtered in R^e using %ld tap rectangular filter and mapped back to spinner seqR1 using Lagrange arc distance",N,M);
  sprintf(filename,"%s_%ld_R2_rect%ld_larc_auto.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_auto(&RxR2rl, M, time1, filename,comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_auto(...,%s,...)\n",filename); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  sprintf(comment,"AUTO-CORRELATION of length %ld spinner seqR1 filtered in R^2 using %ld tap Hanning filter and mapped back to spinner seqR1 using Euclidean metric",N,M);
  sprintf(filename,"%s_%ld_R2_hann%ld_euclid_auto.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_auto(&RxR2he, M, time1, filename,comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_auto(...,%s,...)\n",filename); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  sprintf(comment,"AUTO-CORRELATION of length %ld spinner seqR1 filtered in R^2 using %ld tap Hanning filter and mapped back to spinner seqR1 using Lagrange arc distance",N,M);
  sprintf(filename,"%s_%ld_R2_hann%ld_larc_auto.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_auto(&RxR2hl, M, time1, filename,comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_auto(...,%s,...)\n",filename); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  //----------------------------------------------
  // close log file
  //----------------------------------------------
  plot_close(lptr,time1);
  return 0;
  }

//-----------------------------------------------------------------------------
//! \brief generate TeX plot file for the low pass filtering of fair die sequence
//!        filtered by M tap rectangular filter,
//!        and generate TeX plot file
//-----------------------------------------------------------------------------
int lab_fdie_lp(const unsigned seed, const long N, const long M, const char *basefilename){
  time_t time1; time(&time1);   //starting time stamp (passed to plotting routine)

  fdieseq x(N+M);        //for die sequence
  seqR1   xR1(N+M);      //for R^1 seqR1 mapped from x sequence
  seqR6   xR6(N+M);      //for R^6 seqR1 mapped from x sequence
  seqR1   r(M);          //for M-tap rectangular filter
  seqR1   h(M);          //for M-tap Hanning filter
  seqR1   xR1r(N+M+M-1); //for x mapped to R1 and filtered by r
  seqR1   xR1h(N+M+M-1); //for x mapped to R1 and filtered by h
  seqR6   xR6r(N+M+M-1); //for x mapped to R6 and filtered by r
  seqR6   xR6h(N+M+M-1); //for x mapped to R6 and filtered by h
          //    |  |_______________M-1: convolution with h appends M-1 elements
          //    |__________________ Nm: length of xR1 sequence
  fdieseq xR1re(N+M+M-1);   //for xR1r mapped back using Euclidean metric
  fdieseq xR1he(N+M+M-1);   //for xR1h mapped back using Euclidean metric
  fdieseq xR6re(N+M+M-1);   //for xR6r mapped back using Euclidean metric
  fdieseq xR6rl(N+M+M-1);   //for xR6r mapped back using Lagrange arc distance
  fdieseq xR6he(N+M+M-1);   //for xR6h mapped back using Euclidean metric
  fdieseq xR6hl(N+M+M-1);   //for xR6h mapped back using Lagrange arc distance
                           
  fdieseq xR1reN(N);       //for xR1re with beginning M-1 and ending M-1 elements removed
  fdieseq xR1heN(N);       //for xR1he with beginning M-1 and ending M-1 elements removed
  fdieseq xR6reN(N);       //for xR6re with beginning M-1 and ending M-1 elements removed
  fdieseq xR6rlN(N);       //for xR6rl with beginning M-1 and ending M-1 elements removed
  fdieseq xR6heN(N);       //for xR6he with beginning M-1 and ending M-1 elements removed
  fdieseq xR6hlN(N);       //for xR6hl with beginning M-1 and ending M-1 elements removed
                           
  seqR1   RxR1re(2*N+1);   //for auto-correlation seqR1 of xR1re 
  seqR1   RxR1he(2*N+1);   //for auto-correlation seqR1 of xR1he 
  seqR1   RxR6re(2*N+1);   //for auto-correlation seqR1 of xR6re 
  seqR1   RxR6rl(2*N+1);   //for auto-correlation seqR1 of xR6rl 
  seqR1   RxR6he(2*N+1);   //for auto-correlation seqR1 of xR6he 
  seqR1   RxR6hl(2*N+1);   //for auto-correlation seqR1 of xR6hl 

  char comment[2*1024];          //comment to be passed to plotting function
  char buf[2*1024];          //comment to be passed to plotting function
  long n;
  char filename[1024];
  FILE *lptr; // pointer to log  file

  //----------------------------------------------
  //open log file
  //----------------------------------------------
  printf(         "-----------------------------------------------------------\n");
  sprintf(comment,"Experiment: low pass filter length %ld fair die sequence",N); printf("%s\n",comment);
  printf(         "-----------------------------------------------------------\n");
  sprintf(buf,"%s_%ldm%ld",basefilename,N,M);
  lptr=log_open (buf,time1,comment);

  //----------------------------------------------
  sprintf(buf,"Generate length %ld die sequence using seed value 0x%4x...",N,seed);printofe(lptr,buf,time1);
  //----------------------------------------------
  x.randomize(seed);
  x.list(0,99,"","...\n",1,lptr);
  x.histogram(1,lptr);

  //----------------------------------------------
  sprintf(buf,"Perform mapping operations from die space to R^n...");printofe(lptr,buf,time1);
  //----------------------------------------------
  xR1=x.dietoR1();               //xR1 = x seqR1 mapped to R^1 sequence
  xR6=x.dietoR6();               //xR6 = x die seqR1 mapped to R^6 sequence
  xR1.list(0,24,"sequence mapped to R^1:\n","...\n",1,lptr);  
  xR6.list(0,6,"sequence mapped to R^6:\n","...\n",1,lptr);  

  //----------------------------------------------
  printofe(lptr,"Perform filtering operations...",time1);
  //----------------------------------------------
  r=1.0;       
  h.hanning(); 
  sprintf(buf,"length %ld rectangular sequence before normalization:",M); r.list(buf,"\n",1,lptr);
  sprintf(buf,"length %ld Hanning sequence before normalization:",M);    h.list(buf,"\n",1,lptr);
  r.normalize(); 
  h.normalize(); 
  sprintf(buf,"length %ld rectangular sequence after normalization:",M); r.list(buf,"\n",1,lptr);
  sprintf(buf,"length %ld Hanning sequence after normalization:",M);    h.list(buf,"\n",1,lptr);
  printof(lptr,"xR1*r..."); convolve(&xR1, &r, &xR1r); //yR1 = xR1 filtered by r (convolved with r)
  printof(lptr,"xR1*h..."); convolve(&xR1, &h, &xR1h); //yR1 = xR1 filtered by h (convolved with h)
  printof(lptr,"xR6*r..."); convolve(&xR6, &r, &xR6r); //yR6 = xR6 filtered by r (convolved with r)
  printof(lptr,"xR6*h..."); convolve(&xR6, &h, &xR6h); //yR6 = xR6 filtered by h (convolved with h)
  printof(lptr,"done.\n");
  //----------------------------------------------
  printofe(lptr,"Perform mapping operations from R^n back to die seqR1 space...",time1);
  //----------------------------------------------
  printf("xR1re\n"); xR1re=rdie_R1todie_euclid(xR1r);//xR1r mapped back to fair die seq. using Euclid
  printf("xR1he\n"); xR1he=rdie_R1todie_euclid(xR1h);//xR1h mapped back to fair die seq. using Euclid
  printf("xR6re\n"); xR6re=fdie_R6todie_euclid(xR6r);//xR6r mapped back to fair die seq. using Euclid
  printf("xR6rl\n"); xR6rl=fdie_R6todie_larc  (xR6r);//xR6r mapped back to fair die seq. using Lagrange
  printf("xR6he\n"); xR6he=fdie_R6todie_euclid(xR6h);//xR6h mapped back to fair die seq. using Euclid
  printf("xR6hl\n"); xR6hl=fdie_R6todie_larc  (xR6h);//xR6h mapped back to fair die seq. using Lagrange
  printf("done.\n");
  //----------------------------------------------
  //sprintf(buf,"Compare Euclidean and Lagrange mappings after %ld-tap rectangular filtering...\n  ",M);printofe(lptr,buf,time1);
  //----------------------------------------------
  n=cmp(&xR6re,&xR6rl,0,lptr);
  if(n==0l){sprintf(buf,"Euclidean and Lagrange sequences after using %ld-tap rectangular filter are identical.",M);printofe(lptr,buf,time1);}
  else     {sprintf(buf,"Euclidean and Lagrange sequences after using %ld-tap rectangular filter differ at %ld locations.",M,n);printofe(lptr,buf,time1);}
  //----------------------------------------------
  //sprintf(buf,"Compare Euclidean and Lagrange mappings after %ld-tap Hanning filtering...\n  ",M);printofe(lptr,buf,time1);
  //----------------------------------------------
  n=cmp(&xR6he,&xR6hl,0,lptr);
  if(n==0l){sprintf(buf,"Euclidean and Lagrange sequences after using %ld-tap Hanning filter are identical.",M);printofe(lptr,buf,time1);}
  else     {sprintf(buf,"Euclidean and Lagrange sequences after using %ld-tap Hanning filter differ at %ld locations.",M,n);printofe(lptr,buf,time1);}
  //----------------------------------------------
  sprintf(buf,"Remove beginning %ld element and ending %ld element transition regions introduced by filtering ... ",M-1,M);
  printofe(lptr,buf,time1);
  //----------------------------------------------
  for(n=0;n<N;n++){
    xR1reN.put(n,xR1re.get(n+M-1));
    xR1heN.put(n,xR1he.get(n+M-1));
    xR6reN.put(n,xR6re.get(n+M-1));
    xR6rlN.put(n,xR6rl.get(n+M-1));
    xR6heN.put(n,xR6he.get(n+M-1));
    xR6hlN.put(n,xR6hl.get(n+M-1));
    }
  //----------------------------------------------
  sprintf(buf,"Perform auto-correlation operations...\n");printofe(lptr,buf,time1);
  //----------------------------------------------
  if(xR1reN.Rxxo(&RxR1re,1)){fprintf(stderr,"ERROR computing Rxx for x=xR1re.\n"); return -1;} //auto-correlation seqR1 of truncated xR1re
  if(xR1heN.Rxxo(&RxR1he,1)){fprintf(stderr,"ERROR computing Rxx for x=xR1he.\n"); return -1;} //auto-correlation seqR1 of truncated xR1he
  if(xR6reN.Rxxo(&RxR6re,1)){fprintf(stderr,"ERROR computing Rxx for x=xR6re.\n"); return -1;} //auto-correlation seqR1 of truncated xR6re
  if(xR6rlN.Rxxo(&RxR6rl,1)){fprintf(stderr,"ERROR computing Rxx for x=xR6rl.\n"); return -1;} //auto-correlation seqR1 of truncated xR6rl
  if(xR6heN.Rxxo(&RxR6he,1)){fprintf(stderr,"ERROR computing Rxx for x=xR6he.\n"); return -1;} //auto-correlation seqR1 of truncated xR6he
  if(xR6hlN.Rxxo(&RxR6hl,1)){fprintf(stderr,"ERROR computing Rxx for x=xR6hl.\n"); return -1;} //auto-correlation seqR1 of truncated xR6hl
  // |          |      |____________switch to turn on counting display
  // |          |___________________pointer to output correlation sequence
  // |______________________________input fair die sequence
  printofe(lptr,"done.\n",time1);

  //----------------------------------------------
  //plot seqR1 tex files
  //----------------------------------------------
  sprintf(comment,"length %ld fair die SEQUENCE filtered in R^1 using %ld tap Rectangular filter and mapped backc to fair die seqR1 using Euclidean metric",N,M);
  sprintf(filename,"%s_%ld_R1_rect%ld_euclid_seq.tex",basefilename,N,M);
  printf("Plot file %59s ... ",filename);
  if(plot_ocs_seq((symseq *)&xR1re, M-1, M-1+48,time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_seq(...)\n"); return -1;}
  printf("done.\n");

  sprintf(comment,"length %ld fair die SEQUENCE filtered in R^1 using %ld tap Hanning filter and mapped backc to fair die seqR1 using Euclidean metric",N,M);
  sprintf(filename,"%s_%ld_R1_hann%ld_euclid_seq.tex",basefilename,N,M);
  printf("Plot file %59s ... ",filename);
  if(plot_ocs_seq((symseq *)&xR1he, M-1, M-1+48,time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_seq(...)\n"); return -1;}
  printf("done.\n");

  sprintf(comment,"length %ld fair die SEQUENCE filtered in R^6 using %ld tap rectangular filter and mapped back to fair die seqR1 using Euclidean metric",N,M);
  sprintf(filename,"%s_%ld_R6_rect%ld_euclid_seq.tex",basefilename,N,M);
  printf("Plot file %59s ... ",filename);
  if(plot_ocs_seq((symseq *)&xR6re, M-1, M-1+48,time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_seq(...)\n"); return -1;}
  printf("done.\n");

  sprintf(comment,"length %ld fair die SEQUENCE filtered in R^6 using %ld tap rectangular filter and mapped back to fair die seqR1 using Lagrange arc distance",N,M);
  sprintf(filename,"%s_%ld_R6_rect%ld_larc_seq.tex",basefilename,N,M);
  printf("Plot file %59s ... ",filename);
  if(plot_ocs_seq((symseq *)&xR6rl, M-1, M-1+48,time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_seq(...)\n"); return -1;}
  printf("done.\n");

  sprintf(comment,"length %ld fair die SEQUENCE filtered in R^6 using %ld tap Hanning filter and mapped back to fair die seqR1 using Euclidean metric",N,M);
  sprintf(filename,"%s_%ld_R6_hann%ld_euclid_seq.tex",basefilename,N,M);
  printf("Plot file %59s ... ",filename);
  if(plot_ocs_seq((symseq *)&xR6he, M-1, M-1+48,time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_seq(...)\n"); return -1;}
  printf("done.\n");

  sprintf(comment,"length %ld fair die SEQUENCE filtered in R^6 using %ld tap Hanning filter and mapped back to fair die seqR1 using Lagrange arc distance",N,M);
  sprintf(filename,"%s_%ld_R6_hann%ld_larc_seq.tex",basefilename,N,M);
  printf("Plot file %59s ... ",filename);
  if(plot_ocs_seq((symseq *)&xR6hl, M-1, M-1+48,time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_seq(...)\n"); return -1;}
  printf("done.\n");

  //----------------------------------------------
  //plot histogram tex files
  //----------------------------------------------
  sprintf(comment,"HISTOGRAM of length %ld fair die seqR1 filtered in R^1 using %ld tap Rectangular filter with each tap set to %.6lf",N,M,1.0);
  sprintf(filename,"%s_%ld_R1_rect%ld_euclid_histo.tex",basefilename,N,M);
  printf("Plot file %59s ... ",filename);
  if(plot_ocs_histo((symseq *)&xR1re, M-1, N+M-2, time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_histo(...)\n"); return -1;}
  printf("done.\n");

  sprintf(comment,"HISTOGRAM of length %ld fair die seqR1 filtered in R^1 using %ld tap Hanning filter",N,M);
  sprintf(filename,"%s_%ld_R1_hann%ld_euclid_histo.tex",basefilename,N,M);
  printf("Plot file %59s ... ",filename);
  if(plot_ocs_histo((symseq *)&xR1he, M-1, N+M-2, time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_histo(...)\n"); return -1;}
  printf("done.\n");

  sprintf(comment,"HISTOGRAM of length %ld fair die seqR1 filtered in R^6 using %ld tap rectangular filter and mapped back to fair die seqR1 using Euclidean metric",N,M);
  sprintf(filename,"%s_%ld_R6_rect%ld_euclid_histo.tex",basefilename,N,M);
  printf("Plot file %59s ... ",filename);
  if(plot_ocs_histo((symseq *)&xR6re, M-1, N+M-2, time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_histo(...)\n"); return -1;}
  printf("done.\n");

  sprintf(comment,"HISTOGRAM of length %ld fair die seqR1 filtered in R^6 using %ld tap rectangular filter and mapped back to fair die seqR1 using Lagrange arc distance",N,M);
  sprintf(filename,"%s_%ld_R6_rect%ld_larc_histo.tex",basefilename,N,M);
  printf("Plot file %59s ... ",filename);
  if(plot_ocs_histo((symseq *)&xR6rl, M-1, N+M-2, time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_histo(...)\n"); return -1;}
  printf("done.\n");

  sprintf(comment,"HISTOGRAM of length %ld fair die seqR1 filtered in R^6 using %ld tap Hanning filter and mapped back to fair die seqR1 using Euclidean metric",N,M);
  sprintf(filename,"%s_%ld_R6_hann%ld_euclid_histo.tex",basefilename,N,M);
  printf("Plot file %59s ... ",filename);
  if(plot_ocs_histo((symseq *)&xR6he, M-1, N+M-2, time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_histo(...)\n"); return -1;}
  printf("done.\n");

  sprintf(comment,"HISTOGRAM of length %ld fair die seqR1 filtered in R^6 using %ld tap Hanning filter and mapped back to fair die seqR1 using Lagrange arc distance",N,M);
  sprintf(filename,"%s_%ld_R6_hann%ld_larc_histo.tex",basefilename,N,M);
  printf("Plot file %59s ... ",filename);
  if(plot_ocs_histo((symseq *)&xR6hl, M-1, N+M-2, time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_histo(...)\n"); return -1;}
  printf("done.\n");

  //----------------------------------------------
  //plot autocorrelation tex files
  //----------------------------------------------
  sprintf(comment,"AUTO-CORRELATION of length %ld fair die seqR1 filtered in R^1 using %ld tap Rectangular filter with each tap set to %.6lf",N,M,1.0);
  sprintf(filename,"%s_%ld_R1_rect%ld_euclid_auto.tex",basefilename,N,M);
  printf("Plot file %59s ... ",filename);
  if(plot_ocs_auto(&RxR1re, M, time1, filename,comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_auto(...,%s,...)\n",filename); return -1;}
  printf("done.\n");

  sprintf(comment,"AUTO-CORRELATION of length %ld fair die seqR1 filtered in R^1 using %ld tap Hanning filter",N,M);
  sprintf(filename,"%s_%ld_R1_hann%ld_euclid_auto.tex",basefilename,N,M);
  printf("Plot file %59s ... ",filename);
  if(plot_ocs_auto(&RxR1he, M, time1, filename,comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_auto(...,%s,...)\n",filename); return -1;}
  printf("done.\n");

  sprintf(comment,"AUTO-CORRELATION of length %ld fair die seqR1 filtered in R^6 using %ld tap rectangular filter and mapped back to fair die seqR1 using Euclidean metric",N,M);
  sprintf(filename,"%s_%ld_R6_rect%ld_euclid_auto.tex",basefilename,N,M);
  printf("Plot file %59s ... ",filename);
  if(plot_ocs_auto(&RxR6re, M, time1, filename,comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_auto(...,%s,...)\n",filename); return -1;}
  printf("done.\n");

  sprintf(comment,"AUTO-CORRELATION of length %ld fair die seqR1 filtered in R^e using %ld tap rectangular filter and mapped back to fair die seqR1 using Lagrange arc distance",N,M);
  sprintf(filename,"%s_%ld_R6_rect%ld_larc_auto.tex",basefilename,N,M);
  printf("Plot file %59s ... ",filename);
  if(plot_ocs_auto(&RxR6rl, M, time1, filename,comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_auto(...,%s,...)\n",filename); return -1;}
  printf("done.\n");

  sprintf(comment,"AUTO-CORRELATION of length %ld fair die seqR1 filtered in R^6 using %ld tap Hanning filter and mapped back to fair die seqR1 using Euclidean metric",N,M);
  sprintf(filename,"%s_%ld_R6_hann%ld_euclid_auto.tex",basefilename,N,M);
  printf("Plot file %59s ... ",filename);
  if(plot_ocs_auto(&RxR6he, M, time1, filename,comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_auto(...,%s,...)\n",filename); return -1;}
  printf("done.\n");

  sprintf(comment,"AUTO-CORRELATION of length %ld fair die seqR1 filtered in R^6 using %ld tap Hanning filter and mapped back to fair die seqR1 using Lagrange arc distance",N,M);
  sprintf(filename,"%s_%ld_R6_hann%ld_larc_auto.tex",basefilename,N,M);
  printf("Plot file %59s ... ",filename);
  if(plot_ocs_auto(&RxR6hl, M, time1, filename,comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_auto(...,%s,...)\n",filename); return -1;}
  printf("done.\n");

  //----------------------------------------------
  // close log file
  //----------------------------------------------
  plot_close(lptr,time1);
  return 0;
  }


//-----------------------------------------------------------------------------
//! \brief generate TeX plot file for the high pass filtering of weighted die sequence
//!        filtered by M tap rectangular filter
//!        and generate TeX plot file
//-----------------------------------------------------------------------------
int lab_wdie_hp(const unsigned seed, const long N, const long M, const char *basefilename){
  time_t time1; time(&time1);   //starting time stamp (passed to plotting routine)
  fdieseq x(M*(N+2)-(M-1));   //for fair die sequence
  seqR1 r(M);              //for M-tap rectangular filter
  seqR1 h(M);              //for M-tap Hanning filter

  seqR1   xR1(M*(N+2)-(M-1)); //for x mapped to R^1
  seqR6 xR6(M*(N+2)-(M-1)); //for x mapped to R^6

  seqR1   xR1r(M*(N+2)); //for xR1 filtered by r
  seqR1   xR1h(M*(N+2)); //for xR1 filtered by h
  seqR6 xR6r(M*(N+2)); //for xR6 filtered by r
  seqR6 xR6h(M*(N+2)); //for xR6 filtered by h

  fdieseq  xR1re(M*(N+2)); //for xR1r mapped back to die seqR1 using Euclidean metric
  fdieseq  xR1he(M*(N+2)); //for xR1h mapped back to die seqR1 using Euclidean metric
  fdieseq  xR6re(M*(N+2)); //for xR6r mapped back to die seqR1 using Euclidean metric
  fdieseq  xR6he(M*(N+2)); //for xR6h mapped back to die seqR1 using Euclidean metric
  fdieseq  xR6rl(M*(N+2)); //for xR6r mapped back to die seqR1 using Lagrange  distance
  fdieseq  xR6hl(M*(N+2)); //for xR6h mapped back to die seqR1 using Lagrange  distance
           
  fdieseq  xR1reM(N+2); //for xR1re downsampled by factor of M
  fdieseq  xR1heM(N+2); //for xR1he downsampled by factor of M
  fdieseq  xR6reM(N+2); //for xR6re downsampled by factor of M
  fdieseq  xR6heM(N+2); //for xR6he downsampled by factor of M
  fdieseq  xR6rlM(N+2); //for xR6rl downsampled by factor of M
  fdieseq  xR6hlM(N+2); //for xR6hl downsampled by factor of M
           
  fdieseq  xR1reMN(N); //for xR1reM with first and last element removed
  fdieseq  xR1heMN(N); //for xR1heM with first and last element removed
  fdieseq  xR6reMN(N); //for xR6reM with first and last element removed
  fdieseq  xR6heMN(N); //for xR6heM with first and last element removed
  fdieseq  xR6rlMN(N); //for xR6rlM with first and last element removed
  fdieseq  xR6hlMN(N); //for xR6hlM with first and last element removed

  seqR1 RxR1re(2*N+1); //for auto-correlation of xR1reMN
  seqR1 RxR1he(2*N+1); //for auto-correlation of xR1heMN
  seqR1 RxR6re(2*N+1); //for auto-correlation of xR6reMN
  seqR1 RxR6he(2*N+1); //for auto-correlation of xR6heMN
  seqR1 RxR6rl(2*N+1); //for auto-correlation of xR6rlMN
  seqR1 RxR6hl(2*N+1); //for auto-correlation of xR6hlMN

  char comment[2*1024];          //comment to be passed to plotting function
  char buf[2*1024];          //general purpose buffer
  long n;
  char filename[1024];
  FILE *lptr; // pointer to log  file

  //----------------------------------------------
  //open log file
  //----------------------------------------------
  printf(         "-----------------------------------------------------------\n");
  sprintf(comment,"Experiment: high pass filter length %ld weighted die sequence",N); printf("%s\n",comment);
  printf(         "-----------------------------------------------------------\n");
  sprintf(buf,"%s_%ldm%ld",basefilename,N,M);
  lptr=log_open (buf,time1,comment);

  //----------------------------------------------
  sprintf(buf,"Generate weighted die sequence...");printofe(lptr,buf,time1);
  //----------------------------------------------
  x.randomize(seed,5,5,5,5,75,5);
  x.histogram(1,lptr);
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  //----------------------------------------------
  sprintf(buf,"Map from die space ");printofe(lptr,buf,time1);
  //----------------------------------------------
  printf("to R^1..."); xR1=x.dietoR1();  //xR1 = x seqR1 mapped to R^1 sequence
  printf("to R^6..."); xR6=x.dietoR6();  //xR6 = x die seqR1 mapped to R^6 sequence
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  //----------------------------------------------
  sprintf(buf,"Generate filter coefficients...");printofe(lptr,buf,time1);
  //----------------------------------------------
  r=1.0;       r.lptohp(); printf("\n  r_n=");r.list();
  h.hanning(); h.lptohp(); printf("\n  h_n=");h.list();putchar('\n'); 
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  //----------------------------------------------
  sprintf(buf,"Perform filtering operations...");printofe(lptr,buf,time1);
  //----------------------------------------------
  printf("xR1*r..."); convolve(&xR1, &r, &xR1r); //xR1r = xR1 filtered by r (convolved with r)
  printf("xR1*h..."); convolve(&xR1, &h, &xR1h); //xR1h = xR1 filtered by h (convolved with h)
  printf("xR6*r..."); convolve(&xR6, &r, &xR6r); //xR6r = xR6 filtered by r (convolved with r)
  printf("xR6*h..."); convolve(&xR6, &h, &xR6h); //xR6h = xR6 filtered by h (convolved with h)
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  //----------------------------------------------
  sprintf(buf,"Perform mapping operations from R^n back to die seqR1 space...");printofe(lptr,buf,time1);
  //----------------------------------------------
  xR1re=rdie_R1todie_euclid(xR1r);//xR1r mapped back to fair die seq. using Euclid
  xR1he=rdie_R1todie_euclid(xR1h);//xR1h mapped back to fair die seq. using Euclid
  xR6re=fdie_R6todie_euclid(xR6r);//xR6r mapped back to fair die seq. using Euclid
  xR6rl=fdie_R6todie_larc  (xR6r);//xR6r mapped back to fair die seq. using Lagrange
  xR6he=fdie_R6todie_euclid(xR6h);//xR6h mapped back to fair die seq. using Euclid
  xR6hl=fdie_R6todie_larc  (xR6h);//xR6h mapped back to fair die seq. using Lagrange
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  //----------------------------------------------
  sprintf(buf,"Compare Euclidean and Lagrange mappings after %ld-tap rectangular filtering...\n  ",M);printofe(lptr,buf,time1);
  //----------------------------------------------
  n=cmp(&xR6re,&xR6rl,0,lptr);
  if(n==0l){sprintf(buf,"Euclidean and Lagrange sequences after using %ld-tap rectangular filter are identical.\n",M);printofe(lptr,buf,time1);}
  else     {sprintf(buf,"Euclidean and Lagrange sequences after using %ld-tap rectangular filter differ at %ld locations.\n",M,n);printofe(lptr,buf,time1);}
  //----------------------------------------------
  sprintf(buf,"Compare Euclidean and Lagrange mappings after %ld-tap Hanning filtering...\n  ",M);printofe(lptr,buf,time1);
  //----------------------------------------------
  n=cmp(&xR6he,&xR6hl,0,lptr);
  if(n==0l){sprintf(buf,"Euclidean and Lagrange sequences after using %ld-tap Hanning filter are identical.\n",M);printofe(lptr,buf,time1);}
  else     {sprintf(buf,"Euclidean and Lagrange sequences after using %ld-tap Hanning filter differ at %ld locations.\n",M,n);printofe(lptr,buf,time1);}

  //----------------------------------------------
  sprintf(buf,"Perform down sampling...\n  "); printofe(lptr,buf,time1);
  //----------------------------------------------
  printf("xR1re|%ld...",M); downsample(M, &xR1re, &xR1reM);//xR1reM = xR1re downsampled by M
  printf("xR1he|%ld...",M); downsample(M, &xR1he, &xR1heM);//xR1heM = xR1he downsampled by M
  printf("xR6re|%ld...",M); downsample(M, &xR6re, &xR6reM);//xR6reM = xR6re downsampled by M
  printf("xR6rl|%ld...",M); downsample(M, &xR6rl, &xR6rlM);//xR6rlM = xR6rl downsampled by M
  printf("xR6he|%ld...",M); downsample(M, &xR6he, &xR6heM);//xR6heM = xR6he downsampled by M
  printf("xR6hl|%ld...",M); downsample(M, &xR6hl, &xR6hlM);//xR6hlM = xR6hl downsampled by M
  sprintf(buf,"\ndone.");printofe(lptr,buf,time1);

  //----------------------------------------------
  sprintf(buf,"Remove beginning %ld element and ending %ld element transition regions... ",M-1,M);printofe(lptr,buf,time1);
  //----------------------------------------------
  for(n=0;n<N;n++){
    xR1reMN.put(n,xR1reM.get(n+1));
    xR1heMN.put(n,xR1heM.get(n+1));
    xR6reMN.put(n,xR6reM.get(n+1));
    xR6rlMN.put(n,xR6rlM.get(n+1));
    xR6heMN.put(n,xR6heM.get(n+1));
    xR6hlMN.put(n,xR6hlM.get(n+1));
    }
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  //----------------------------------------------
  sprintf(buf,"Perform auto-correlation operations...");printofe(lptr,buf,time1);
  //----------------------------------------------
  if(xR1reMN.Rxxo(&RxR1re,1)){fprintf(stderr,"ERROR computing Rxx for x=xR1re.\n"); return -1;} //auto-correlation seqR1 of xR1reMN
  if(xR1heMN.Rxxo(&RxR1he,1)){fprintf(stderr,"ERROR computing Rxx for x=xR1he.\n"); return -1;} //auto-correlation seqR1 of xR1heMN
  if(xR6reMN.Rxxo(&RxR6re,1)){fprintf(stderr,"ERROR computing Rxx for x=xR6re.\n"); return -1;} //auto-correlation seqR1 of xR6reMN
  if(xR6rlMN.Rxxo(&RxR6rl,1)){fprintf(stderr,"ERROR computing Rxx for x=xR6rl.\n"); return -1;} //auto-correlation seqR1 of xR6rlMN
  if(xR6heMN.Rxxo(&RxR6he,1)){fprintf(stderr,"ERROR computing Rxx for x=xR6he.\n"); return -1;} //auto-correlation seqR1 of xR6heMN
  if(xR6hlMN.Rxxo(&RxR6hl,1)){fprintf(stderr,"ERROR computing Rxx for x=xR6hl.\n"); return -1;} //auto-correlation seqR1 of xR6hlMN
  // |            |      |____________switch to turn on counting display
  // |            |___________________pointer to output correlation sequence
  // |________________________________input fair die sequence
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  //----------------------------------------------
  //plot seqR1 tex files
  //----------------------------------------------
  sprintf(comment,"length %ld fair die SEQUENCE filtered in R^1 using %ld tap Rectangular filter and mapped backc to fair die seqR1 using Euclidean metric",N,M);
  sprintf(filename,"%s_%ld_R1_rect%ld_euclid_seq.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_seq((symseq *)&xR1reM, 0,50, time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_seq(...)\n"); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  sprintf(comment,"length %ld fair die SEQUENCE filtered in R^1 using %ld tap Hanning filter and mapped backc to fair die seqR1 using Euclidean metric",N,M);
  sprintf(filename,"%s_%ld_R1_hann%ld_euclid_seq.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_seq((symseq *)&xR1heM, 0,50, time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_seq(...)\n"); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  sprintf(comment,"length %ld fair die SEQUENCE filtered in R^6 using %ld tap rectangular filter and mapped back to fair die seqR1 using Euclidean metric",N,M);
  sprintf(filename,"%s_%ld_R6_rect%ld_euclid_seq.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_seq((symseq *)&xR6reM, 0,50, time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_seq(...)\n"); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  sprintf(comment,"length %ld fair die SEQUENCE filtered in R^6 using %ld tap rectangular filter and mapped back to fair die seqR1 using Lagrange arc distance",N,M);
  sprintf(filename,"%s_%ld_R6_rect%ld_larc_seq.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_seq((symseq *)&xR6rlM, 0,50, time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_seq(...)\n"); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  sprintf(comment,"length %ld fair die SEQUENCE filtered in R^6 using %ld tap Hanning filter and mapped back to fair die seqR1 using Euclidean metric",N,M);
  sprintf(filename,"%s_%ld_R6_hann%ld_euclid_seq.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_seq((symseq *)&xR6heM, 0,50, time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_seq(...)\n"); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  sprintf(comment,"length %ld fair die SEQUENCE filtered in R^6 using %ld tap Hanning filter and mapped back to fair die seqR1 using Lagrange arc distance",N,M);
  sprintf(filename,"%s_%ld_R6_hann%ld_larc_seq.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_seq((symseq *)&xR6hlM, 0,50, time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_seq(...)\n"); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);


  //----------------------------------------------
  //plot histogram tex files
  //----------------------------------------------
  sprintf(comment,"HISTOGRAM of length %ld fair die seqR1 filtered in R^1 using %ld tap Rectangular filter with each tap set to %.6lf",N,M,1.0);
  sprintf(filename,"%s_%ld_R1_rect%ld_euclid_histo.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_histo((symseq *)&xR1reMN, time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_histo(...)\n"); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  sprintf(comment,"HISTOGRAM of length %ld fair die seqR1 filtered in R^1 using %ld tap Hanning filter",N,M);
  sprintf(filename,"%s_%ld_R1_hann%ld_euclid_histo.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_histo((symseq *)&xR1heMN, time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_histo(...)\n"); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  sprintf(comment,"HISTOGRAM of length %ld fair die seqR1 filtered in R^6 using %ld tap rectangular filter and mapped back to fair die seqR1 using Euclidean metric",N,M);
  sprintf(filename,"%s_%ld_R6_rect%ld_euclid_histo.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_histo((symseq *)&xR6reMN, time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_histo(...)\n"); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  sprintf(comment,"HISTOGRAM of length %ld fair die seqR1 filtered in R^6 using %ld tap rectangular filter and mapped back to fair die seqR1 using Lagrange arc distance",N,M);
  sprintf(filename,"%s_%ld_R6_rect%ld_larc_histo.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_histo((symseq *)&xR6rlMN, time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_histo(...)\n"); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  sprintf(comment,"HISTOGRAM of length %ld fair die seqR1 filtered in R^6 using %ld tap Hanning filter and mapped back to fair die seqR1 using Euclidean metric",N,M);
  sprintf(filename,"%s_%ld_R6_hann%ld_euclid_histo.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_histo((symseq *)&xR6heMN, time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_histo(...)\n"); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  sprintf(comment,"HISTOGRAM of length %ld fair die seqR1 filtered in R^6 using %ld tap Hanning filter and mapped back to fair die seqR1 using Lagrange arc distance",N,M);
  sprintf(filename,"%s_%ld_R6_hann%ld_larc_histo.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_histo((symseq *)&xR6hlMN, time1, "die", filename, comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_histo(...)\n"); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  //----------------------------------------------
  //plot autocorrelation tex files
  //----------------------------------------------
  sprintf(comment,"AUTO-CORRELATION of length %ld fair die seqR1 filtered in R^1 using %ld tap Rectangular filter with each tap set to %.6lf",N,M,1.0);
  sprintf(filename,"%s_%ld_R1_rect%ld_euclid_auto.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_auto(&RxR1re, M, time1, filename,comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_auto(...,%s,...)\n",filename); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  sprintf(comment,"AUTO-CORRELATION of length %ld fair die seqR1 filtered in R^1 using %ld tap Hanning filter",N,M);
  sprintf(filename,"%s_%ld_R1_hann%ld_euclid_auto.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_auto(&RxR1he, M, time1, filename,comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_auto(...,%s,...)\n",filename); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  sprintf(comment,"AUTO-CORRELATION of length %ld fair die seqR1 filtered in R^6 using %ld tap rectangular filter and mapped back to fair die seqR1 using Euclidean metric",N,M);
  sprintf(filename,"%s_%ld_R6_rect%ld_euclid_auto.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_auto(&RxR6re, M, time1, filename,comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_auto(...,%s,...)\n",filename); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  sprintf(comment,"AUTO-CORRELATION of length %ld fair die seqR1 filtered in R^e using %ld tap rectangular filter and mapped back to fair die seqR1 using Lagrange arc distance",N,M);
  sprintf(filename,"%s_%ld_R6_rect%ld_larc_auto.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_auto(&RxR6rl, M, time1, filename,comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_auto(...,%s,...)\n",filename); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  sprintf(comment,"AUTO-CORRELATION of length %ld fair die seqR1 filtered in R^6 using %ld tap Hanning filter and mapped back to fair die seqR1 using Euclidean metric",N,M);
  sprintf(filename,"%s_%ld_R6_hann%ld_euclid_auto.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_auto(&RxR6he, M, time1, filename,comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_auto(...,%s,...)\n",filename); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  sprintf(comment,"AUTO-CORRELATION of length %ld fair die seqR1 filtered in R^6 using %ld tap Hanning filter and mapped back to fair die seqR1 using Lagrange arc distance",N,M);
  sprintf(filename,"%s_%ld_R6_hann%ld_larc_auto.tex",basefilename,N,M);
  sprintf(buf,"Plot file %50s ... ",filename);printofe(lptr,buf,time1);
  if(plot_ocs_auto(&RxR6hl, M, time1, filename,comment,lptr)){fprintf(stderr,"ERROR using plot_ocs_auto(...,%s,...)\n",filename); return -1;}
  sprintf(buf,"done.");printofe(lptr,buf,time1);

  //----------------------------------------------
  // close log file
  //----------------------------------------------
  plot_close(lptr,time1);
  return 0;
  }



//-----------------------------------------------------------------------------
//! \brief map nonstationary die sequence into R^1, C^1, and R^6 and 
//!        use DFT to analyze
//-----------------------------------------------------------------------------
int lab_die_nonstat34(const unsigned seed, const long N, const long M, const int vx, const char *basefilename){
  time_t time1; time(&time1);   //starting time stamp (passed to plotting routine)
  const int vA=vx,vB=vx,vD=vx,vE=vx,vF=vx, vC=100-5*vx;
  const int wA=vx,wB=vx,wC=vx,wE=vx,wF=vx, wD=100-5*vx;
  double threshold;
  double yR1n;
  double yC1n;
  double yR6n;
  long count;
  fdieseq x(N);    //pseudo-random die sequence
  seqR1   xR1(N);  //   xR1 = x mapped to R^1
  seqC1   xC1(N);  //   xC1 = x mapped to C^1
  seqR6   xR6(N);  //   xR6 = x mapped to R^6
  seqC1  DxR1(N);  //  DxR1 =  DFT(xR1) 
  seqC1  DxC1(N);  //  DxC1 =  DFT(xC1) 
  seqC6  DxR6(N);  //  DxR6 =  DFT(xR6) 
  seqR1 mDxR1(N);  // mDxR1 = |DFT(xR1)|
  seqR1 mDxC1(N);  // mDxC1 = |DFT(xC1)|
  seqR1 mDxR6(N);  // mDxR6 = |DFT(xR6)|

  long n;
  vectC6 y6n;
  char comment[1024];      //comment to be passed to plotting function
  char buf[4*1024];          //general purpose buffer
  char filename[1024];
  double pmargin;
  FILE *lptr; // pointer to log  file

  //----------------------------------------------
  //open log file
  //----------------------------------------------
  printf(         "-----------------------------------------------------------\n");
  sprintf(comment,"Experiment: perform Fourier analysis on length %ld die sequence with periodic distribution",N); printf("%s\n",comment);
  printf(         "-----------------------------------------------------------\n");
  sprintf(filename,"%s_%02d%02d_%ldm%ld",basefilename,vx,vC,N,M);
  lptr=log_open (filename,time1,comment);
  sprintf(buf,"open log file %s to %s",filename,comment);

  //----------------------------------------------
  sprintf(buf,"generate length %ld pseudo-random die sequence\n",N);printofe(lptr,buf,time1);
  //----------------------------------------------
  x.prngseed(seed);
  for(n=0;((n+2)*M/2)<=N;n+=2){
    x.randomize(    n*M/2, (n+1)*M/2-1, vA,vB,vC,vD,vE,vF); 
    x.randomize((n+1)*M/2, (n+2)*M/2-1, wA,wB,wC,wD,wE,wF);
    sprintf(buf,"die weight*100  for [x_{%4ld}...x_{%4ld}]: (%3d %3d %3d %3d %3d %3d)\n",n*M/2,    (n+1)*M/2-1, vA,vB,vC,vD,vE,vF);
    printof(lptr,buf);
    sprintf(buf,"die weight*100  for [x_{%4ld}...x_{%4ld}]: (%3d %3d %3d %3d %3d %3d)\n",(n+1)*M/2,(n+2)*M/2-1, wA,wB,wC,wD,wE,wF);
    printof(lptr,buf);
    }
  for(n=0;n<10;n++){
    sprintf(buf,"[x_{%4ld}...x_{%4ld}] = ",n*M/2,(n+1)*M/2-1);
    x.list(n*M/2,(n+1)*M/2-1,buf,"\n",1,lptr);
    } 

  x.histogram(1,lptr);
  x.histogram(0*M/2, 1*M/2-1, 1, lptr);
  x.histogram(1*M/2, 2*M/2-1, 1, lptr);
  x.histogram(2*M/2, 3*M/2-1, 1, lptr);
  x.histogram(3*M/2, 4*M/2-1, 1, lptr);

  //----------------------------------------------
  printofe(lptr,"Map from die space to R^1, C^1, and R^6",time1);
  //----------------------------------------------
  xR1=x.dietoR1pam();
  xC1=x.dietoC1();
  xR6=x.dietoR6();
  xR1.list (0,49,"x->R^1 values n=0..49:\n","\n",1,lptr);
  xC1.list (0,49,"x->C^1 values n=0..49:\n","\n",1,lptr);
  xR6.list1(0,49,"x->R^6 values n=0..49:\n","\n",1,lptr);

  //----------------------------------------------
  printofe(lptr,"Perform DFT operations...",time1);
  //----------------------------------------------
  printf("  DFT:R^1-->C^1 :"); dft(&xR1, &DxR1 ); printf("\n"); //  DxR1 =  DFT(xR1)
  printf("  DFT:C^1-->C^1 :"); dft(&xC1, &DxC1 ); printf("\n"); //  DxC1 =  DFT(xC1)
  printf("  DFT:R^6-->C^6 :"); dft(&xR6, &DxR6 ); printf("\n"); //  DxR6 =  DFT(xR6)
  mag(&DxR1, &mDxR1); // mDxR1 = |DFT(xR1)|
  mag(&DxC1, &mDxC1); // mDxC1 = |DFT(xC1)|
  mag(&DxR6, &mDxR6); // mDxR6 = |DFT(xR6)|
  mDxR1.list(0,49,"|DFT(R^1 sequence)| = \n","...\n",1,lptr);
  mDxC1.list(0,49,"|DFT(C^1 sequence)| = \n","...\n",1,lptr);
  mDxR6.list(0,49,"|DFT(R^6 sequence)| = \n","...\n",1,lptr);

//DxR1.list(N/2-10,N/2,"","\n",1,lptr);
//DxR1.list(N/2,N/2+10,"","\n",1,lptr);
//mDxR1.list(N/2-10,N/2,"","\n",1,lptr);
//mDxR1.list(N/2,N/2+10,"","\n",1,lptr);
//
//DxC1.list(N/2-10,N/2,"","\n",1,lptr);
//DxC1.list(N/2,N/2+10,"","\n",1,lptr);
//mDxC1.list(N/2-10,N/2,"","\n",1,lptr);
//mDxC1.list(N/2,N/2+10,"","\n",1,lptr);
//
//DxR6.list(N/2-5,N/2-1,"list DxR6:...\n","\n",1,lptr);
//DxR6.list(N/2,N/2,"","\n",1,lptr);
//DxR6.list(N/2+1,N/2+5,"","...\n",1,lptr);
//mDxR6.list(N/2-10,N/2,"","\n",1,lptr);
//mDxR6.list(N/2+1,N/2+10,"","\n",1,lptr);

  //----------------------------------------------
  printofe(lptr,"Perform amplitude analysis...",time1);
  //----------------------------------------------
  n=N/M; 
  yR1n=mDxR1.get(n);
  yC1n=mDxC1.get(n);
  yR6n=mDxR6.get(n);
  
  pmargin=1.00;

  threshold = pmargin*yR1n;
  count=mDxR1.gt(threshold,0L,N/2-1);
  sprintf(buf,"\nnumber of |DFT(R^1 sequence)| values > %lf = %.2lf*|DFT(x,%ld)| =\n  %.2lf*%lf in [x_n|n=0,...,%ld] is %ld",threshold,pmargin,n,pmargin,yC1n,N/2-1,count);
  printof(lptr,buf);
  mDxR1.gt(threshold,0L,N/2-1,":\n","\n",0,lptr);

  threshold = pmargin*yC1n;
  count=mDxC1.gt(threshold,0L,N-1);
  sprintf(buf,"\nnumber of |DFT(C^1 sequence)| values > %lf = %.2lf*|DFT(x,%ld)| =\n  %.2lf*%lf in [x_n|n=0,...,%ld] is %ld",threshold,pmargin,n,pmargin,yC1n,N-1,count);
  printof(lptr,buf);
  mDxC1.gt(threshold,0L,N-1,":\n","\n",0,lptr);

  threshold = pmargin*yR6n;
  count=mDxR6.gt(threshold,0L,N/2-1);
  sprintf(buf,"\nnumber of |DFT(R^6 sequence)| values > %lf = %.2lf*|DFT(x,%ld)| =\n  %.2lf*%lf in [x_n|n=0,...,%ld] is %ld",threshold,pmargin,n,pmargin,yC1n,N/2-1,count);
  printof(lptr,buf);
  mDxR6.gt(threshold,0L,N/2-1,":\n","\n",0,lptr);

  pmargin=0.90;
  threshold = pmargin*yR6n;
  count=mDxR6.gte(threshold,0L,N/2-1);
  sprintf(buf,"\nnumber of |DFT(R^6 sequence)| values >= %lf = %.2lf*|DFT(x,%ld)| =\n  %.2lf*%lf in [x_n|n=0,...,%ld] is %ld",threshold,pmargin,n,pmargin,yC1n,N/2-1,count);
  printof(lptr,buf);
  mDxR6.gte(threshold,0L,N/2-1,":\n","\n",0,lptr);

  pmargin=0.75;
  threshold = pmargin*yR6n;
  count=mDxR6.gte(threshold,0L,N/2-1);
  sprintf(buf,"\nnumber of |DFT(R^6 sequence)| values >= %lf = %.2lf*|DFT(x,%ld)| =\n  %.2lf*%lf in [x_n|n=0,...,%ld] is %ld",threshold,pmargin,n,pmargin,yC1n,N/2-1,count);
  printof(lptr,buf);
  mDxR6.gte(threshold,0L,N/2-1,":\n","\n",0,lptr);

  pmargin=0.65;
  threshold = pmargin*yR6n;
  count=mDxR6.gte(threshold,0L,N/2-1);
  sprintf(buf,"\nnumber of |DFT(R^6 sequence)| values >= %lf = %.2lf*|DFT(x,%ld)| =\n  %.2lf*%lf in [x_n|n=0,...,%ld] is %ld",threshold,pmargin,n,pmargin,yC1n,N/2-1,count);
  printof(lptr,buf);
  mDxR6.gte(threshold,0L,N/2-1,":\n","\n",0,lptr);

  printofe(lptr,"\n...done.",time1);

  //----------------------------------------------
  printof(lptr,"generate DFT magnitude TeX files ...\n");
  //----------------------------------------------
  sprintf(filename,"%sR1_%02d%02d_%ldm%ld.tex",basefilename,vx,vC,N,M);
  sprintf(buf,"   |DFT(R^1 seq.)|: \"%s\"\n",filename);printof(lptr,buf);
  plot_dft_seq( &mDxR1, 0, N/2, -3, 7, filename, comment);
  sprintf(filename,"%sC1_%02d%02d_%ldm%ld.tex",basefilename,vx,vC,N,M);
  sprintf(buf,"   |DFT(C^1 seq.)|: \"%s\"\n",filename);printof(lptr,buf);
  plot_dft_seq( &mDxC1, 0, N-1, -5, 5, filename, comment);
  sprintf(filename,"%sR6_%02d%02d_%ldm%ld.tex",basefilename,vx,vC,N,M);
  sprintf(buf,"   |DFT(R^6 seq.)|: \"%s\"\n",filename);printof(lptr,buf);
  if(N<=2000)plot_dft_seq( &mDxR6, 0, N/2, -3, 3, filename, comment);
  else       plot_dft_seq( &mDxR6, 0, 500, -4, 4, filename, comment);

  //----------------------------------------------
  // close log file
  //----------------------------------------------
  sprintf(buf,"lab_die_nonstat34 experiment complete"); printofe(lptr,buf,time1);
  plot_close(lptr,time1);
  return 0;
  }

//-----------------------------------------------------------------------------
//! \brief die sequence edge detection
//-----------------------------------------------------------------------------
int lab_die_edge(const unsigned seed, const long N, const long M, const long Mh, const int vx, const char *basefilename){
  time_t time1; time(&time1);   //starting time stamp (passed to plotting routine)
  const long D=10;
  const int vA=vx,vB=vx,vD=vx,vE=vx,vF=vx, vC=100-5*vx;

  fdieseq x(N);        //for die sequence
  seqR1   xR1(N);      //for R^1 seqR1 mapped from x sequence
  seqC1   xC1(N);      //for R^1 seqR1 mapped from x sequence
  seqR6   xR6(N);      //for R^6 seqR6 mapped from x sequence
  seqR1   h(Mh);        //for M-tap Hanning filter
  seqR1   xR1h( N+Mh-1); //for x mapped to R1 and filtered by h
  seqC1   xC1h( N+Mh-1); //for x mapped to R6 and filtered by h
  seqR6   xR6h( N+Mh-1); //for x mapped to R6 and filtered by h
  seqR1   xR1hm(N+Mh-1); //magnitude for x mapped to R1 and filtered by h
  seqR1   xC1hm(N+Mh-1); //magnitude for x mapped to R1 and filtered by h
  seqR1   xR6hm(N+Mh-1); //magnitude for x mapped to R6 and filtered by h

  double threshold,yn,pmargin,max,rmsval;
  long n,count;
  char comment[2*1024];      //comment to be passed to plotting function
  char buf[128];          //general purpose buffer
  FILE *lptr; // pointer to log  file

  //----------------------------------------------
  //open log file
  //----------------------------------------------
  printf(         "-----------------------------------------------------------\n");
  sprintf(comment,"Experiment: perform wavelet analysis on length %ld non-stationary die sequence",N); printf("%s\n",comment);
  printf(         "-----------------------------------------------------------\n");
  sprintf(buf,"%s_%ldm%ld_h%02ld_%02d%02d",basefilename,N,M,Mh,vx,100-5*vx);
  lptr=log_open (buf,time1,comment);

  //----------------------------------------------
  // generate die sequence
  //----------------------------------------------
  sprintf(buf,"weights*100 for [x_{%3ld}...x_{%3ld}]: (%2d %2d %2d %2d %2d %2d)\n",M, 2*M-1,vA,vB,vC,vD,vE,vF);printof(lptr,buf);
  x.randomize(seed);
  x.randomize(M,2*M-1, vA,vB,vC,vD,vE,vF); 
  x.histogram(1,lptr);
  x.histogram(  0,   M-1, 1, lptr);
  x.histogram(  M, 2*M-1, 1, lptr);
  x.histogram(2*M, N-1,   1, lptr);

  //----------------------------------------------
  //filter coefficients
  //----------------------------------------------
  for(n=0;n<Mh/2;n++)h.put(n,+1.0/(double)Mh);
  for(   ;n<Mh;  n++)h.put(n,-1.0/(double)Mh);
  h.list("Haar wavelet sequence:\n","\n\n",1,lptr);

  //----------------------------------------------
  sprintf(buf,"Perform filtering R^1-->R^1...");printofe(lptr,buf,time1);
  //----------------------------------------------
  xR1=x.dietoR1pam();
  convolve(&xR1, &h, &xR1h);// xR1h=h*xR1;
  mag(&xR1h,&xR1hm);
  xR1hm.truncate(N);
  n=M+Mh/2;
  sprintf(buf,"|W(R^1 sequence)| from n=%ld..%ld:\n",n-50,n+50);
  xR1hm.list(n-50,n+50,buf,"\n",1,lptr);
  pmargin=1.0;
  yn = xR1hm.get(n);
  threshold = pmargin*yn;
  count=xR1hm.gte(threshold,0,N-1);
  sprintf(buf,"\nnumber of |W(R^1 sequence)| values >= %lf = %.2lf*|W(x,%ld)| =\n  %.2lf*%lf in [x_n|n=0,...,%ld] is %ld\n",
                                                     threshold,pmargin,n,        pmargin,yn,N-1,count);
  printof(lptr,buf);
  sprintf(buf,"%sR1_%ldm%ld_h%02ld_%02d%02d_D%ld.tex",basefilename,N,M,Mh,vx,100-5*vx,D);
  plot_R1_seq(&xR1hm, 0, N-1, D, 0, 0.5, 160, 30, buf, comment);

  //----------------------------------------------
  sprintf(buf,"Perform filtering C^1-->C^1...");printofe(lptr,buf,time1);
  //----------------------------------------------
  xC1=x.dietoC1();
  convolve(&xC1, &h, &xC1h); //xC1h = xC1*h;
  mag(&xC1h,&xC1hm);
  n=M+Mh/2;
  sprintf(buf,"|W(C^1 sequence)| from n=%ld..%ld:\n",n-50,n+50);
  xC1hm.list(n-50,n+50,buf,"\n",1,lptr);
  max=xC1hm.max(n-50,n+50);
  sprintf(buf,"max|W(C^1 sequence)| from n=%ld..%ld is %lf:\n",n-50,n+50,max);printof(lptr,buf);
  rmsval=rms(&xC1hm,Mh,N-1);
  sprintf(buf,"rms|W(C^1 sequence)| from n=%ld..%ld is %lf:\n",Mh,N-1,rmsval);printof(lptr,buf);
  pmargin=1.0;
  n=M+Mh/2;
  yn = xC1hm.get(n);
  threshold = pmargin*yn;
  count=xC1hm.gte(threshold,0,N-1);
  sprintf(buf,"\nnumber of |W(C^1 sequence)| values >= %lf = %.2lf*|W(x,%ld)| =\n  %.2lf*%lf in [x_n|n=0,...,%ld] is %ld\n",
                                                     threshold,pmargin,n,        pmargin,yn,N-1,count);
  printof(lptr,buf);
  xC1hm.gte(threshold,0,N-1,buf,"\n",0,lptr);
  sprintf(buf,"%sC1_%ldm%ld_h%02ld_%02d%02d_D%ld.tex",basefilename,N,M,Mh,vx,100-5*vx,D);
  xC1hm.truncate(N);
  plot_R1_seq(&xC1hm, 0, N-1, D, 0, 0.30, 160, 30, buf, comment);

  //----------------------------------------------
  sprintf(buf,"Perform filtering R^1-->R^6...");printofe(lptr,buf,time1);
  //----------------------------------------------
  xR6=x.dietoR6();
  sprintf(buf,"[x'_%ld..x'_%ld] of length %ld die sequence mapped to R^6:\n",M-10,M+9,N); 
  convolve(&xR6, &h, &xR6h);//xR6h = xR6*h;
  mag(&xR6h,&xR6hm);
  n=M+Mh/2;
  sprintf(buf,"|W(R^6 sequence)| from n=%ld..%ld:\n",n-50,n+50);
  xR6hm.list(n-50,n+50,buf,"\n",1,lptr);
  max=xR6hm.max(n-50,n+50);
  sprintf(buf,"max|W(R^6 sequence)| from n=%ld..%ld is %lf:\n",n-50,n+50,max);printof(lptr,buf);
  rmsval=rms(&xR6hm,Mh,N-1);
  sprintf(buf,"rms|W(R^6 sequence)| from n=%ld..%ld is %lf:\n",Mh,N-1,rmsval);printof(lptr,buf);
  pmargin=1.0;
  n=M+Mh/2;
  yn = xR6hm.get(n);
  threshold = pmargin*yn;
  count=xR6hm.gte(threshold,0,N-1);
  sprintf(buf,"\nnumber of |W(R^4 sequence)| values >= %lf = %.2lf*|W(x,%ld)| =\n  %.2lf*%lf in [x_n|n=0,...,%ld] is %ld\n",
                                                     threshold,pmargin,n,        pmargin,yn,N-1,count);
  printof(lptr,buf);
  xR6hm.gte(threshold,0,N-1,buf,"\n",0,lptr);
  sprintf(buf,"%sR6_%ldm%ld_h%02ld_%02d%02d_D%ld.tex",basefilename,N,M,Mh,vx,100-5*vx,D);
  xR6hm.truncate(N);
  plot_R1_seq(&xR6hm, 0, N-1, D, 0, 0.30, 160, 30, buf, comment);
  //----------------------------------------------
  // Done.
  //----------------------------------------------
  //----------------------------------------------
  // close log file
  //----------------------------------------------
  sprintf(buf,"lab_die_edge experiment complete"); printofe(lptr,buf,time1);
  plot_close(lptr,time1);
  return 0;
  }


//-----------------------------------------------------------------------------
//! \brief DFT analysis of a DNA sequence
//-----------------------------------------------------------------------------
int lab_dna_dft(const char *datafilename, const char *basefilename){
  time_t time1; time(&time1);   //starting time stamp (passed to plotting routine)
  const long N=numsym_fasta_file(datafilename);
  const long M=8;
  const long hres1=0;
  const long hres2=N/3;
  const long hres3=N/2;
  const double plotmin=-1.0;
  double plotmax=10.0;
  dnaseq x(N);
  char header[1024];
  double rmsval,maxval;
  if(N==-1)return -1;

  seqR1   xR1(N);      //R^1 sequence mapped from x sequence
  seqC1   xC1(N);      //C^1 sequence mapped from x sequence
  seqR4   xR4(N);      //R^4 sequence mapped from x sequence
  seqC1   DxR1(N);     //DFT of R^1 sequence mapped from x sequence
  seqC1   DxC1(N);     //DFT of C^1 sequence mapped from x sequence
  seqC4   DxR4(N);     //DFT of R^4 sequence mapped from x sequence
  seqR1   mDxR1(N);    //magnitude of DFT of R^1 sequence mapped from x sequence
  seqR1   mDxC1(N);    //magnitude of DFT of C^1 sequence mapped from x sequence
  seqR1   mDxR4(N);    //magnitude of DFT of R^4 sequence mapped from x sequence

  char comment[2*1024];      //comment to be passed to plotting function
  char buf[128];          //general purpose buffer
  char filename[128];          //general purpose buffer
  FILE *lptr; // pointer to log  file

  //----------------------------------------------
  // read data from FASTA file
  //----------------------------------------------
  read_fasta_file(datafilename,header,&x);
  x.list(0,99);//putchar('\n');
  
  //----------------------------------------------
  //special test block: periodic data 
  //comment out this block after testing
  //----------------------------------------------
  //const long P=100;//period length
  //const long S=P/4;//sub-period
  //for(n=0;n<N;n+=P){//initialize period by period
  //  x.put(n+0*S,n+1*S-1,'C');
  //  x.put(n+1*S,n+2*S-1,'G');
  //  x.put(n+2*S,n+3*S-1,'T');
  //  x.put(n+3*S,n+4*S-1,'A');
  //  }
  //x.list();

  //----------------------------------------------
  //open log file
  //----------------------------------------------
  printf(         "-----------------------------------------------------------\n");
  sprintf(comment,"Experiment: perform Fourier analysis on length %ld DNA sequence",N); printf("%s\n",comment);
  printf(         "-----------------------------------------------------------\n");
  sprintf(buf,"%s",basefilename);
  lptr=log_open (buf,time1,comment);

  //----------------------------------------------
  sprintf(buf,"Perform DFT operation R^1-->C^1...\n");printof(lptr,buf);
  //----------------------------------------------
  xR1=x.dnatoR1pam();
  sprintf(buf,"[x'_%d..x'_%d] of length %ld dna sequence mapped to R^1 using PAM mapping:\n",0,50,N); 
  dft(&xR1,&DxR1);
  mag(&DxR1,&mDxR1);
  rmsval=rms(&mDxR1,1,N/2-1);
  maxval=mDxR1.max(1,N/2-1);
  sprintf(buf,"\nR^1 rms = %lf, max=%lf\n",rmsval,maxval);
  printof(lptr,buf);
  mDxR1.list(N/2-10,N/2,"N/2-10...\n","\n",lptr);
  mDxR1.list(N/2,N/2+10,"...N/2+10\n","\n",lptr);
  mDxR1.list(N/2-10,N/2,"N/2-10...\n","\n",lptr);
  mDxR1.list(N/2,N/2+10,"...N/2+10\n","\n",lptr);
  sprintf(filename,"%sR1.tex",basefilename);
  plotmax = ceil(10*log10(mDxR1.max(N/10,9*N/10)));
  plot_dft_seq( &mDxR1, 0, N/2, plotmin, plotmax, M,hres1,hres2,hres3, filename, comment);
  //----------------------------------------------
  sprintf(buf,"Perform DFT operation C^1-->C^1...");printof(lptr,buf);
  //----------------------------------------------
  xC1=x.dnatoC1();
  sprintf(buf,"[x'_%d..x'_%d] of length %ld dna sequence mapped to C^1 using QPSK mapping:\n",0,10,N); 
  xC1.list(0,10,buf,"\n",lptr);
  dft(&xC1,&DxC1);
  mag(&DxC1,&mDxC1);
  rmsval=rms(&mDxC1,1,N-1);
  maxval=mDxC1.max(1,N-1);
  sprintf(buf,"\nC^1 rms = %lf, max=%lf\n",rmsval,maxval);
  printof(lptr,buf);
  sprintf(filename,"%sC1.tex",basefilename);
  plotmax = ceil(10*log10(mDxC1.max(N/10,9*N/10)));
  plot_dft_seq( &mDxC1, 0, N-1, plotmin, plotmax, M,hres1,hres2,hres3, filename, comment);
  sprintf(buf,"done.");printof(lptr,buf);
 
  //----------------------------------------------
  sprintf(buf,"Perform DFT operation R^4-->C^4...\n");printof(lptr,buf);
  //----------------------------------------------
  xR4=x.dnatoR4();
  sprintf(buf,"[x'_%d..x'_%d] of length %ld dna sequence mapped to R^4 :\n",0,10,N);  
  xR4.list(0,10,buf,"\n",lptr);
  dft(&xR4,&DxR4);
  mag(&DxR4,&mDxR4);
  rmsval=rms(&mDxR4,1,N/2-1);
  maxval=mDxR4.max(1,N/2-1);
  sprintf(buf,"\nR^4 rms = %lf, max=%lf\n",rmsval, maxval);
  printof(lptr,buf);
  sprintf(filename,"%sR4.tex",basefilename);
  plotmax = ceil(10*log10(mDxR4.max(N/10,9*N/10)));
  plot_dft_seq( &mDxR4, 0, N/2, plotmin, plotmax, M,hres1,hres2,hres3, filename, comment);

  //----------------------------------------------
  // close log file
  //----------------------------------------------
  plot_close(lptr,time1);
  return 0;
  }


//---------------------------------------------------------------------------
// DFT analysis of a non-stationary artificial DNA sequence
//---------------------------------------------------------------------------
int lab_dna_nonstatCT(const unsigned seed, const long N, const long M, const int vx, const double plotmin, const double plotmax, const char *basefilename){
  time_t time1;time(&time1);   //starting time stamp (passed to plotting routine)
  const int vA=vx,vG=vx,vT=vx, vC=100-3*vx;
  const int wA=vx,wC=vx,wT=vx, wG=100-3*vx;
  double threshold;
  dnaseq x(N);
  seqR1   xR1(N);      //R^1 sequence mapped from x sequence
  seqC1   xC1(N);      //C^1 sequence mapped from x sequence
  seqR4   xR4(N);      //R^4 sequence mapped from x sequence
  seqC1   DxR1(N);     //DFT of R^1 sequence mapped from x sequence
  seqC1   DxC1(N);     //DFT of C^1 sequence mapped from x sequence
  seqC4   DxR4(N);     //DFT of R^4 sequence mapped from x sequence
  seqR1   ymag(N);     //magnitude of DFT sequences
  long n;
  complex ynR1;
  complex ynC1;
  vectC4  ynR4;
  const double pmargin=1.0;
  char filename[128];          //general purpose buffer
  char comment[2*1024];      //comment to be passed to plotting function
  char buf[2*1024];          //general purpose buffer
  FILE *lptr; // pointer to log  file

  //----------------------------------------------
  //open log file
  //----------------------------------------------
  printf(         "-----------------------------------------------------------\n");
  sprintf(comment,"Experiment: perform DFT analysis on nonstationary length %ld artificial dna sequence\n%%  with length %ld periodic distribution",N,M);printf("%s",comment);
  printf(         "-----------------------------------------------------------\n");
  sprintf(buf,"%s_%ldm%ld",basefilename,N,M);
  lptr=log_open (buf,time1,comment);

  //----------------------------------------------
  //sprintf(buf,"length %ld pseudo-random die sequence ",N); printof(lptr,buf);
  //----------------------------------------------
  sprintf(buf,"dna weight values for [x_{%3d}...x_{%3ld}]: (%d %d %d %d)\n",0,  M/2-1,vA,vC,vG,vT);printof(lptr,buf);
  sprintf(buf,"dna weight values for [x_{%3ld}...x_{%3ld}]: (%d %d %d %d)\n",M/2,M-1,  wA,wC,wG,wT);printof(lptr,buf);
  x.seed(seed);
  for(n=0;((n+2)*M/2)<=N;n+=2){
    x.randomize(    n*M/2, (n+1)*M/2-1, vA,vC,vG,vT); 
    x.randomize((n+1)*M/2, (n+2)*M/2-1, wA,wC,wG,wT);
    }
  sprintf(buf,"[x_0..x_499] of length %ld pseudo-random die sequence:\n",N); x.list(0,499,buf,"\n",lptr);
  x.histogram(                1,lptr);                  
  x.histogram(0*M/2, 1*M/2-1, 1,lptr);  
  x.histogram(1*M/2, 2*M/2-1, 1,lptr);  
  x.histogram(2*M/2, 3*M/2-1, 1,lptr);  
  x.histogram(3*M/2, 4*M/2-1, 1,lptr);  

  //----------------------------------------------
  sprintf(buf,"Perform DFT operation R^1-->C^1...\n");printof(lptr,buf);
  //----------------------------------------------
  xR1=x.dnatoR1pam();
  sprintf(buf,"[x'_%d..x'_%d] of length %ld dna sequence mapped to R^1 using PAM mapping:\n",0,50,N); 
  xR1.list(0,100,buf,"\n",lptr);
  dft(&xR1,&DxR1);
  mag(&DxR1,&ymag);
  n=N/M; 
  ynR1=dftn(&xR1,n);
  threshold = pmargin*ynR1.mag();
  sprintf(buf,"\n|DFT(x)| values > %lf = %.2lf*|DFT(x,%ld)| = %.2lf*%lf\nin domain [0,N/2-1]=[0,%ld]:\n",threshold,pmargin,n,pmargin,ynR1.mag(),N/2-1);
  ymag.gt(threshold,0,N/2-1,buf,"\n",lptr);
  sprintf(filename,"%sR1.tex",basefilename);
  plot_dft_seq( &ymag, 0, N/12, plotmin, plotmax, filename, comment);

  //----------------------------------------------
  sprintf(buf,"Perform DFT operation C^1-->C^1...");printof(lptr,buf);
  //----------------------------------------------
  xC1=x.dnatoC1();
  sprintf(buf,"[x'_%d..x'_%d] of length %ld dna sequence mapped to C^1 using QPSK mapping:\n",0,10,N); 
  xC1.list(0,10,buf,"\n",lptr);
  dft(&xC1,&DxC1);
  mag(&DxC1,&ymag);
  n=N/M; 
  ynC1=dftn(&xC1,n);
  threshold = pmargin*ynC1.mag();
  sprintf(buf,"\n|DFT(x)| values > %lf = %.2lf*|DFT(x,%ld)| = %.2lf*%lf\nin domain [0,N-1]=[0,%ld]:\n",threshold,pmargin,n,pmargin,ynC1.mag(),N-1);
  ymag.gt(threshold,0,N-1,buf,"\n",lptr);
  sprintf(filename,"%sC1.tex",basefilename);
  plot_dft_seq( &ymag, 0, N/12, plotmin, plotmax, filename, comment);
  sprintf(buf,"done.");printof(lptr,buf);
 
  //----------------------------------------------
  sprintf(buf,"Perform DFT operation R^4-->C^4...\n");printof(lptr,buf);
  //----------------------------------------------
  xR4=x.dnatoR4();
  sprintf(buf,"[x'_%d..x'_%d] of length %ld dna sequence mapped to R^4 :\n",0,10,N);  
  xR4.list(0,10,buf,"\n",lptr);
  dft(&xR4,&DxR4);
  mag(&DxR4,&ymag);
  n=N/M; 
  ynR4=dftn(&xR4,n);
  threshold = pmargin*ynR4.mag();
  sprintf(buf,"\n|DFT(x)| values >= %lf = %.2lf*|DFT(x,%ld)| = %.2lf*%lf\nin domain [0,N/2-1]=[0,%ld]:\n",threshold,pmargin,n,pmargin,ynR4.mag(),N/2-1);
  ymag.gte(threshold,0,N/2-1,buf,"\n",lptr);
  sprintf(filename,"%sR4.tex",basefilename);
  plot_dft_seq( &ymag, 0, N/12, plotmin, plotmax, filename, comment);

  //----------------------------------------------
  // close log file
  //----------------------------------------------
  plot_close(lptr,time1);
  return 0;
  }

//-----------------------------------------------------------------------------
//! \brief dna sequence edge detection
//-----------------------------------------------------------------------------
int lab_dna_edge(const unsigned seed, const long N, const long M, const long Mh, const int vx, const char *basefilename){
  time_t time1; time(&time1);   //starting time stamp (passed to plotting routine)
  const long D=5;
  const int vA=vx,vT=vx,vG=vx,vC=100-3*vx;

  dnaseq  x(N);          //for dna sequence
  seqR1   xR1(N);        //for R^1 seqR1 mapped from x sequence
  seqC1   xC1(N);        //for C^1 seqC1 mapped from x sequence
  seqR4   xR4(N);        //for R^4 seqR4 mapped from x sequence
  seqR1   h(Mh);         //for M-tap Hanning filter
  seqR1   xR1h( N+Mh-1); //for x mapped to R1 and filtered by h
  seqC1   xC1h( N+Mh-1); //for x mapped to R1 and filtered by h
  seqR4   xR4h( N+Mh-1); //for x mapped to R4 and filtered by h
  seqR1   xR1hm(N+Mh-1); //magnitude for x mapped to R1 and filtered by h
  seqR1   xC1hm(N+Mh-1); //magnitude for x mapped to R1 and filtered by h
  seqR1   xR4hm(N+Mh-1); //magnitude for x mapped to R4 and filtered by h

  double threshold,rmsval,yn,pmargin,max1,max2,max;
  long n,n2,count;
  char comment[2*1024];      //comment to be passed to plotting function
  char buf[128];          //general purpose buffer
  FILE *lptr; // pointer to log  file

  //----------------------------------------------
  //open log file
  //----------------------------------------------
  printf(         "-----------------------------------------------------------\n");
  sprintf(comment,"Experiment: Edge detection on nonstationary length %ld artificial dna sequence with length %ld periodic distribution",N,M);printf("%s",comment);
  printf(         "-----------------------------------------------------------\n");
  sprintf(buf,"%s_%ldm%ld_h%02ld_%02d%02d",basefilename,N,M,Mh,vx,100-3*vx);
  lptr=log_open (buf,time1,comment);

  //----------------------------------------------
  // generate dna sequence
  //----------------------------------------------
  sprintf(buf,"dna weight*100 for [x_{%3d}...x_{%3ld}]: (%d %d %d %d)\n",0, M-1,vA,vC,vT,vG);printof(lptr,buf);
  x.randomize(seed);
  x.randomize(M,2*M-1, vA,vC,vG,vT); 
  x.histogram(1,lptr);              
  x.histogram(  0,   M-1, 1, lptr);  
  x.histogram(  M, 2*M-1, 1, lptr);  
  x.histogram(2*M,   N-1, 1, lptr);  

  //----------------------------------------------
  //filter coefficients
  //----------------------------------------------
  for(n=0;n<Mh/2;n++)h.put(n,+1.0/(double)Mh);
  for(   ;n<Mh;  n++)h.put(n,-1.0/(double)Mh);
  h.list("fiter coefficients h_n = \n","\n",1,lptr);

  //----------------------------------------------
  sprintf(buf,"Perform filtering R^1-->R^1...");printofe(lptr,buf,time1);
  //----------------------------------------------
  xR1=x.dnatoR1pam();
  convolve(&xR1, &h, &xR1h);// xR1h=h*xR1;
  mag(&xR1h,&xR1hm);
  xR1hm.truncate(N);
  n=M+Mh/2;

  max1=xR1hm.max(n-50,n+50);
  sprintf(buf,"max|W(R^1 sequence)| from n=%ld..%ld is %lf:\n",n-50,n+50,max1);printof(lptr,buf);
  n2=2*M+Mh/2;
  max2=xR1hm.max(n2-50,n2+50);
  sprintf(buf,"max|W(R^1 sequence)| from n=%ld..%ld is %lf:\n",n2-50,n2+50,max2);printof(lptr,buf);
  rmsval=rms(&xR1hm,Mh,N-1);
  sprintf(buf,"rms|W(R^1 sequence)| from n=%ld..%ld is %lf:\n",Mh,N-1,rmsval);printof(lptr,buf);

  sprintf(buf,"|W(R^1 sequence)| from n=%ld..%ld:\n",n-50,n+50);
  xR1hm.list(n-50,n+50,buf,"\n",1,lptr);
  sprintf(buf,"|W(R^1 sequence)| from n=%ld..%ld:\n",n2-50,n2+50);
  xR1hm.list(n2-50,n2+50,buf,"\n",1,lptr);

  pmargin=1.0;
  //yn = xR1hm.get(n);
  yn = max1;
  threshold = pmargin*yn;
  count=xR1hm.gte(threshold,0,N-1);
  sprintf(buf,"\nnumber of |W(R^1 sequence)| values >= %lf = %.2lf*|W(x,n)| =\n  %.2lf*%lf in [x_n|n=0,...,%ld] is %ld\n",
                                                     threshold,pmargin,        pmargin,yn,N-1,count);
  printof(lptr,buf);
  xR1hm.gte(threshold,0,N-1,"\n","\n",1,lptr);
  sprintf(buf,"%sR1_%ldm%ld_h%02ld_%02d%02d_D%ld.tex",basefilename,N,M,Mh,vx,100-3*vx,D);
  plot_R1_seq(&xR1hm, 0, N-1, D, 0, 0.150, 160, 30, buf, comment);

  //----------------------------------------------
  sprintf(buf,"Perform filtering C^1-->C^1...");printofe(lptr,buf,time1);
  //----------------------------------------------
  xC1=x.dnatoC1();
  convolve(&xC1, &h, &xC1h); //xC1h = xC1*h;
  mag(&xC1h,&xC1hm);
  n=M+Mh/2;
  sprintf(buf,"|W(C^1 sequence)| from n=%ld..%ld:\n",n-50,n+50);
  xC1hm.list(n-50,n+50,buf,"\n",1,lptr);
  max=xC1hm.max(n-50,n+50);
  sprintf(buf,"max|W(C^1 sequence)| from n=%ld..%ld is %lf:\n",n-50,n+50,max);printof(lptr,buf);
  rmsval=rms(&xC1hm,Mh,N-1);
  sprintf(buf,"rms|W(C^1 sequence)| from n=%ld..%ld is %lf:\n",Mh,N-1,rmsval);printof(lptr,buf);
  pmargin=1.0;
  n=M+Mh/2;
  //yn = xC1hm.get(n);
  yn = max;
  threshold = pmargin*yn;
  count=xC1hm.gte(threshold,0,N-1);
  sprintf(buf,"\nnumber of |W(C^1 sequence)| values >= %lf = %.2lf*|W(x,%ld)| =\n  %.2lf*%lf in [x_n|n=0,...,%ld] is %ld\n",
                                                     threshold,pmargin,n,        pmargin,yn,N-1,count);
  xC1hm.gte(threshold,0,N-1,buf,"\n",0,lptr);
  sprintf(buf,"%sC1_%ldm%ld_h%02ld_%02d%02d_D%ld.tex",basefilename,N,M,Mh,vx,100-3*vx,D);
  xC1hm.truncate(N);
  plot_R1_seq(&xC1hm, 0, N-1, D, 0, 0.30, 160, 30, buf, comment);

  //----------------------------------------------
  sprintf(buf,"Perform filtering R^1-->R^4...");printofe(lptr,buf,time1);
  //----------------------------------------------
  xR4=x.dnatoR4();
  convolve(&xR4, &h, &xR4h);//xR4h = xR4*h;
  mag(&xR4h,&xR4hm);
  n=M+Mh/2;
  sprintf(buf,"|W(R^4 sequence)| from n=%ld..%ld:\n",n-50,n+50);
  xR4hm.list(n-50,n+50,buf,"\n",1,lptr);
  max=xR4hm.max(n-50,n+50);
  sprintf(buf,"max|W(R^4 sequence)| from n=%ld..%ld is %lf:\n",n-50,n+50,max);printof(lptr,buf);
  rmsval=rms(&xR4hm,Mh,N-1);
  sprintf(buf,"rms|W(R^4 sequence)| from n=%ld..%ld is %lf:\n",Mh,N-1,rmsval);printof(lptr,buf);
  pmargin=1.0;
  n=M+Mh/2;
  yn = xR4hm.get(n);
  threshold = pmargin*yn;
  count=xR4hm.gte(threshold,0,N-1);
  sprintf(buf,"\nnumber of |W(R^4 sequence)| values >= %lf = %.2lf*|W(x,%ld)| =\n  %.2lf*%lf in [x_n|n=0,...,%ld] is %ld\n",
                                                     threshold,pmargin,n,        pmargin,yn,N-1,count);
  xR4hm.gte(threshold,0,N-1,buf,"\n",0,lptr);
  sprintf(buf,"%sR4_%ldm%ld_h%02ld_%02d%02d_D%ld.tex",basefilename,N,M,Mh,vx,100-3*vx,D);
  xR4hm.truncate(N);
  plot_R1_seq(&xR4hm, 0, N-1, D, 0, 0.30, 160, 30, buf, comment);


  ////----------------------------------------------
  //sprintf(buf,"Perform filtering R^1-->R^4...\n");printof(lptr,buf);
  ////----------------------------------------------
  //xR4=x.dnatoR4();
  //sprintf(buf,"[x'_%ld..x'_%ld] of length %ld dna sequence mapped to R^4:\n",M-10,M+9,N); 
  //xR4.list(M-10,M+9,buf,"\n",lptr);
  //convolve(&xR4, &h, &xR4h); //xR4h = xR4*h;
  //mag(&xR4h,&xR4hm);
  ////min(&xR4h,&xR4hm);
  ////max(&xR4h,&xR4hm);
  ////downsample(D,&xR4hm,&xR4hmD);
  //threshold = 0.8*xR4hm.get(M+Mh/2);
  //sprintf(buf,"\nvalues >= %lf \n",threshold);
  //xR4hm.gte(threshold,0,N-1,buf,"\n",lptr);
  //sprintf(buf,"%sR4_%ldm%ld_h%02ld_%02d%02d_D%ld.tex",basefilename,N,M,Mh,vx,100-3*vx,D);
  //xR4hm.truncate(N);
  //plot_R1_seq(&xR4hm, 0, N-1, D, 0, 0.30, 160, 30, buf, comment);
  //----------------------------------------------
  // close log file
  //----------------------------------------------
  plot_close(lptr,time1);
  return 0;
  }

//-----------------------------------------------------------------------------
//! \brief dna sequence averaging
//-----------------------------------------------------------------------------
int lab_dna_averaging(const long Mh, const char *datafilename, const char *basefilename){
  time_t time1; time(&time1);   //starting time stamp (passed to plotting routine)
  const long D=20;
  const long N=numsym_fasta_file(datafilename);
  if(N==-1)return -1;

  dnaseq  x(N);         //for dna sequence
  seqR1   xR1(N);       //for R^1 seqR1 mapped from x sequence
  seqC1   xC1(N);       //for C^1 seqC1 mapped from x sequence
  seqR4   xR4(N);       //for R^4 seqR4 mapped from x sequence
  seqR1   h(Mh);        //for M-tap Hanning filter
  seqR1   xR1h(N+Mh-1); //for x mapped to R1 and filtered by h
  seqC1   xC1h(N+Mh-1); //for x mapped to R1 and filtered by h
  seqR4   xR4h(N+Mh-1); //for x mapped to R4 and filtered by h

  long n;
  char comment[2*1024];      //comment to be passed to plotting function
  char header[1024];
  char buf[128];          //general purpose buffer
  FILE *lptr; // pointer to log  file

  //----------------------------------------------
  // read data from FASTA file
  //----------------------------------------------
  read_fasta_file(datafilename,header,&x);
  x.list(0,99);//putchar('\n');
  
  //----------------------------------------------
  //open log file
  //----------------------------------------------
  printf(         "-----------------------------------------------------------\n");
  sprintf(comment,"Experiment: map length %ld dna sequence from %s\n%% and perform edge detection using Haar wavelet",N,datafilename);printf("%s",comment);
  printf(         "-----------------------------------------------------------\n");
  sprintf(buf,"%s_hs%ld",basefilename,Mh);
  lptr=log_open (buf,time1,comment);

  //----------------------------------------------
  //histograms
  //----------------------------------------------
  x.histogram( 9000,11000,1,lptr); 
  x.histogram(23000,25000,1,lptr); 
  x.histogram(34000,36000,1,lptr); 
  x.histogram(42000,44000,1,lptr); 

  //----------------------------------------------
  //filter coefficients
  //----------------------------------------------
  for(n=0;n<Mh;n++)h.put(n,1.0/(double)Mh);
  h.list("fiter coefficients h_n = ","\n",1,lptr);

  //----------------------------------------------
  sprintf(buf,"Perform filtering R^1-->R^1...\n");printof(lptr,buf);
  //----------------------------------------------
  xR1=x.dnatoR1bin();
  //xR1=x.dnatoR1pam();
  sprintf(buf,"[x'_%ld..x'_%ld] of length %ld dna sequence mapped to R^1 using binary mapping:\n",N/2-50,N/2+49,N); 
  xR1.list(N/2-50,N/2+49,buf,"\n",lptr);
  convolve(&xR1, &h, &xR1h);//xR1h = h*xR1;
  //mag(&xR1h,&xR1hm);
  //threshold = xR1hm.get(4079);
  //sprintf(buf,"\nvalues >= %lf \n",threshold);
  //xR1hm.gte(threshold,0,N-1,buf,"\n",lptr);
  sprintf(buf,"%sR1b_hs%02ld.tex",basefilename,Mh);
  xR1h.truncate(N);
  plot_R1_seq(&xR1h, 0, N-1, D, -0.50, 0.50, 160, 30, buf, comment);

//  //----------------------------------------------
//  sprintf(buf,"Perform filtering C^1-->C^1...\n");printof(lptr,buf);
//  //----------------------------------------------
//  xC1=x.dnatoC1();
//  sprintf(buf,"[x'_%ld..x'_%ld] of length %ld dna sequence mapped to C^1 using QPSK mapping:\n",N/2-10,N/2+9,N); 
//  xC1.list(N/2-10,N/2+9,buf,"\n",lptr);
//  xC1h = xC1*h;
//  mag(&xC1h,&xC1hm);
//  threshold = xC1hm.get(N/2+Mh/2);
//  sprintf(buf,"\nvalues >= %lf \n",threshold);
//  xC1hm.gte(threshold,0,N-1,buf,"\n",lptr);
//  sprintf(buf,"%sC1_h%02ld.tex",basefilename,Mh);
//  xC1hm.truncate(N);
//  plot_R1_seq(&xC1hm, 0, N-1, D, 0, 0.150, 160, 30, buf, comment);
//
//  //----------------------------------------------
//  sprintf(buf,"Perform filtering R^1-->R^4...\n");printof(lptr,buf);
//  //----------------------------------------------
//  xR4=x.dnatoR4();
//  sprintf(buf,"[x'_%ld..x'_%ld] of length %ld dna sequence mapped to R^4:\n",N/2-10,N/2+9,N); 
//  xR4.list(N/2-10,N/2+9,buf,"\n",lptr);
//  xR4h = xR4*h;
//  mag(&xR4h,&xR4hm);
//  threshold = 0.8*xR4hm.get(N/2+Mh/2);
//  sprintf(buf,"\nvalues >= %lf \n",threshold);
//  xR4hm.gte(threshold,0,N-1,buf,"\n",lptr);
//  sprintf(buf,"%sR4_h%02ld.tex",basefilename,Mh);
//  xR4hm.truncate(N);
//  plot_R1_seq(&xR4hm, 0, N-1, D, 0, 0.150, 160, 30, buf, comment);
  //----------------------------------------------
  // close log file
  //----------------------------------------------
  plot_close(lptr,time1);
  return 0;
  }


//-----------------------------------------------------------------------------
//! \brief dna sequence edge detection
//-----------------------------------------------------------------------------
int lab_dna_edge(const long Mh, const char *datafilename, const char *basefilename){
  time_t time1; time(&time1);   //starting time stamp (passed to plotting routine)
  const long D=20;
  const long N=numsym_fasta_file(datafilename);
  if(N==-1)return -1;

  dnaseq  x(N);          //for dna sequence
  seqR1   xR1(N);        //for R^1 seqR1 mapped from x sequence
  seqC1   xC1(N);        //for C^1 seqC1 mapped from x sequence
  seqR4   xR4(N);        //for R^4 seqR4 mapped from x sequence
  seqR1   h(Mh);         //for M-tap Hanning filter
  seqR1   xR1h( N+Mh-1); //for x mapped to R1 and filtered by h
  seqC1   xC1h( N+Mh-1); //for x mapped to R1 and filtered by h
  seqR4   xR4h( N+Mh-1); //for x mapped to R4 and filtered by h
  seqR1   xR1hm(N+Mh-1); //magnitude for x mapped to R1 and filtered by h
  seqR1   xC1hm(N+Mh-1); //magnitude for x mapped to R1 and filtered by h
  seqR1   xR4hm(N+Mh-1); //magnitude for x mapped to R4 and filtered by h

  long n;
  double threshold;
  char comment[2*1024];      //comment to be passed to plotting function
  char header[1024];
  char buf[128];          //general purpose buffer
  FILE *lptr; // pointer to log  file

  //----------------------------------------------
  // read data from FASTA file
  //----------------------------------------------
  read_fasta_file(datafilename,header,&x);
  x.list(0,99);//putchar('\n');
  
  //----------------------------------------------
  //open log file
  //----------------------------------------------
  printf(         "-----------------------------------------------------------\n");
  sprintf(comment,"Experiment: map length %ld dna sequence from %s\n%% and perform edge detection using Haar wavelet",N,datafilename);printf("%s",comment);
  printf(         "-----------------------------------------------------------\n");
  sprintf(buf,"%s_h%02ld",basefilename,Mh);
  lptr=log_open (buf,time1,comment);

  //----------------------------------------------
  //histograms
  //----------------------------------------------
  x.histogram( 9000,11000,1,lptr); 
  x.histogram(23000,25000,1,lptr); 
  x.histogram(34000,36000,1,lptr); 
  x.histogram(42000,44000,1,lptr); 

  //----------------------------------------------
  //filter coefficients
  //----------------------------------------------
  for(n=0;n<Mh/2;n++)h.put(n,-1.0/(double)Mh);
  for(   ;n<Mh;  n++)h.put(n,+1.0/(double)Mh);
  h.list("fiter coefficients h_n = ","\n",1,lptr);

  //----------------------------------------------
  sprintf(buf,"Perform filtering R^1-->R^1...");printofe(lptr,buf,time1);
  //----------------------------------------------
  xR1=x.dnatoR1pam();
  sprintf(buf,"[x'_%ld..x'_%ld] of length %ld dna sequence mapped to R^1 using PAM mapping:\n",N/2-50,N/2+49,N); 
  xR1.list(N/2-50,N/2+49,buf,"\n",lptr);
  convolve(&xR1, &h, &xR1h);//xR1h = h*xR1;
  mag(&xR1h,&xR1hm);
  threshold = xR1hm.get(4079);
  sprintf(buf,"\nvalues >= %lf \n",threshold);
  xR1hm.gte(threshold,0,N-1,buf,"\n",lptr);
  sprintf(buf,"%sR1_h%02ld.tex",basefilename,Mh);
  xR1hm.truncate(N);
  plot_R1_seq(&xR1hm, 0, N-1, D, 0, 0.150, 160, 30, buf, comment);

  //xR1=x.dnatoR1pam();
  //convolve(&xR1, &h, &xR1h);// xR1h=h*xR1;
  //mag(&xR1h,&xR1hm);
  //xR1hm.truncate(N);
  //n=M+Mh/2;
  //sprintf(buf,"|W(R^1 sequence)| from n=%ld..%ld:\n",n-50,n+50);
  //xR1hm.list(n-50,n+50,buf,"\n",1,lptr);
  //pmargin=1.0;
  //yn = xR1hm.get(n);
  //threshold = pmargin*yn;
  //count=xR1hm.gte(threshold,0,N-1);
  //sprintf(buf,"\nnumber of |W(R^1 sequence)| values >= %lf = %.2lf*|W(x,%ld)| =\n  %.2lf*%lf in [x_n|n=0,...,%ld] is %ld\n",
  //                                                   threshold,pmargin,n,        pmargin,yn,N-1,count);
  //sprintf(buf,"%sR1_%ldm%ld_h%02ld_%02d%02d_D%ld.tex",basefilename,N,M,Mh,vx,100-5*vx,D);
  //plot_R1_seq(&xR1hm, 0, N-1, D, 0, 0.150, 160, 30, buf, comment);

  //----------------------------------------------
  sprintf(buf,"Perform filtering C^1-->C^1...\n");printof(lptr,buf);
  //----------------------------------------------
  xC1=x.dnatoC1();
  sprintf(buf,"[x'_%ld..x'_%ld] of length %ld dna sequence mapped to C^1 using QPSK mapping:\n",N/2-10,N/2+9,N); 
  xC1.list(N/2-10,N/2+9,buf,"\n",lptr);
  convolve(&xC1, &h, &xC1h); //xC1h = xC1*h;
  mag(&xC1h,&xC1hm);
  threshold = xC1hm.get(N/2+Mh/2);
  sprintf(buf,"\nvalues >= %lf \n",threshold);
  xC1hm.gte(threshold,0,N-1,buf,"\n",lptr);
  sprintf(buf,"%sC1_h%02ld.tex",basefilename,Mh);
  xC1hm.truncate(N);
  plot_R1_seq(&xC1hm, 0, N-1, D, 0, 0.150, 160, 30, buf, comment);

  //----------------------------------------------
  sprintf(buf,"Perform filtering R^1-->R^4...\n");printof(lptr,buf);
  //----------------------------------------------
  xR4=x.dnatoR4();
  sprintf(buf,"[x'_%ld..x'_%ld] of length %ld dna sequence mapped to R^4:\n",N/2-10,N/2+9,N); 
  xR4.list(N/2-10,N/2+9,buf,"\n",lptr);
  convolve(&xR4, &h, &xR4h); //xR4h = xR4*h;
  mag(&xR4h,&xR4hm);
  threshold = 0.8*xR4hm.get(N/2+Mh/2);
  sprintf(buf,"\nvalues >= %lf \n",threshold);
  xR4hm.gte(threshold,0,N-1,buf,"\n",lptr);
  sprintf(buf,"%sR4_h%02ld.tex",basefilename,Mh);
  xR4hm.truncate(N);
  plot_R1_seq(&xR4hm, 0, N-1, D, 0, 0.150, 160, 30, buf, comment);
  //----------------------------------------------
  // close log file
  //----------------------------------------------
  plot_close(lptr,time1);
  return 0;
  }



